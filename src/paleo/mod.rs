//! Paleo-Genomic & Paleoproteomic Reconstruction (`paleo` feature).
//!
//! A mathematically rigorous Tri-Layer engine for reconstructing coding loci, resurrecting
//! inactivated pseudogenes, and planning de-extinction chassis edits for extinct taxa:
//! - **Layer 1 ([`asr`])**: Continuous-Time Markov Chain (JC69 / K80) Ancestral Sequence
//!   Reconstruction via Felsenstein's post-order pruning algorithm.
//! - **Layer 2 ([`biomolecule`])**: Bayesian paleoproteomic codon back-translation (`COL1A1`, `ENAM`)
//!   and Briggs post-mortem deamination-aware ancient DNA (`aDNA`) integration.
//! - **Layer 2.5 ([`resurrection`])**: Outgroup-bracketed pseudogene resurrection (repairing
//!   premature stop codons and frameshift defects in extant descendants).
//! - **Layer 3 ([`phenotype`])**: Commutative, group-discounted Bayesian phenotypic & isotope
//!   selection profiling.
//! - **Synthesis ([`engine`])**: [`PaleoReconstructionEngine`] producing [`PaleoGenomeDraft`]
//!   and prioritized [`ChassisGenomicEdit`] plans.

pub mod asr;
pub mod biomolecule;
pub mod engine;
pub mod phenotype;
pub mod resurrection;
pub mod types;

pub use asr::{ExtantSequence, FelsensteinAsr, SubstitutionModel};
pub use biomolecule::{
    AncientDnaAligner, AncientDnaRead, PaleoproteomeTranslator, PeptideFragment,
};
pub use engine::{
    ChassisGenomicEdit, EditEffectType, PaleoGenomeDraft, PaleoReconstructionEngine,
    ReconstructionMetrics,
};
pub use phenotype::{
    BoneHistologyType, CorrelationGroup, FunctionalPathway, IsotopeEvidence, OsteologyEvidence,
    PathwaySelectionProfile, PhenotypicEvidence,
};
pub use resurrection::{PseudogeneResurrector, RepairEvent, RepairType, ResurrectedGene};
pub use types::{
    base_to_index, index_to_base, translate_dna_codon, BasePosterior, LocusProvenance, PaleoError,
    ReconstructedLocus, SiteSource, BASES,
};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::phylo::PhyloNode;
    use crate::seq::Dna;

    fn sample_archosaur_tree() -> PhyloNode {
        // ((Chicken:68, Falcon:68):160, Alligator:228)
        let birds = PhyloNode::internal(
            vec![
                PhyloNode::leaf("Gallus_gallus", 68.0),
                PhyloNode::leaf("Falco_peregrinus", 68.0),
            ],
            160.0,
        );
        let alligator = PhyloNode::leaf("Alligator_mississippiensis", 228.0);
        PhyloNode::internal(vec![birds, alligator], 0.0)
    }

    #[test]
    fn test_felsenstein_asr_and_branch_decay() {
        let tree = sample_archosaur_tree();
        let seqs = vec![
            ExtantSequence::new("Gallus_gallus", Dna::new(b"ATGGGTCCCGCT").unwrap()),
            ExtantSequence::new("Falco_peregrinus", Dna::new(b"ATGGGTCCCGCT").unwrap()),
            ExtantSequence::new(
                "Alligator_mississippiensis",
                Dna::new(b"ATGGGTCCCGCA").unwrap(),
            ),
        ];
        let model = SubstitutionModel::default_vertebrate_coding();

        let locus_crown = FelsensteinAsr::reconstruct_locus(
            &tree,
            &seqs,
            &model,
            "Archosauria_anc",
            "COL1A1",
            0.0,
        )
        .unwrap();

        let locus_diverged = FelsensteinAsr::reconstruct_locus(
            &tree,
            &seqs,
            &model,
            "Tyrannosaurus_rex",
            "COL1A1",
            65.0,
        )
        .unwrap();

        assert_eq!(locus_crown.to_consensus_dna().sequence, b"ATGGGTCCCGCT");
        // Conserved sites (0..11) have higher confidence than divergent 3rd codon site (11)
        assert!(locus_crown.sites[0].max_prob() > locus_crown.sites[11].max_prob());
        // Longer unanchored extinct branch increases Shannon entropy monotonically
        assert!(locus_diverged.mean_entropy_bits > locus_crown.mean_entropy_bits);
    }

    #[test]
    fn test_paleoproteomic_back_translation_and_codon_bias() {
        let tree = sample_archosaur_tree();
        // Divergent third positions and uncertain middle codon
        let seqs = vec![
            ExtantSequence::new("Gallus_gallus", Dna::new(b"GGTCCCCCT").unwrap()),
            ExtantSequence::new("Falco_peregrinus", Dna::new(b"GGACCCCCT").unwrap()),
            ExtantSequence::new(
                "Alligator_mississippiensis",
                Dna::new(b"GGGCCTCCA").unwrap(),
            ),
        ];
        let model = SubstitutionModel::default_vertebrate_coding();
        let mut locus = FelsensteinAsr::reconstruct_locus(
            &tree,
            &seqs,
            &model,
            "Tyrannosaurus_rex",
            "COL1A1",
            80.0,
        )
        .unwrap();

        let pre_conf = locus.mean_posterior_confidence;

        // Fossil collagen peptide: Gly-Hyp-Ala ("GPA", written with 'O' for hydroxyproline!)
        // Wait: in extant sequences codon 2 was Pro ("CCT"/"CCA"), while fossil peptide proves T. rex had Ala ('A') at codon 2!
        let peptide = PeptideFragment::new(0, b"GOA", 0.98);
        // Reference genome strongly prefers "CCC" for Pro and "GCC" for Ala
        let codon_ref = Dna::new(b"CCCCCCCCCGCCGCCGCC").unwrap();

        let anchored =
            PaleoproteomeTranslator::anchor_with_peptides(&mut locus, &[peptide], Some(&codon_ref))
                .unwrap();
        assert_eq!(anchored, 3);
        assert!(locus.mean_posterior_confidence > pre_conf);

        let prot = locus.translate_consensus().unwrap();
        assert_eq!(prot.sequence, b"GPA");
        // Check that codon 2 (bases 6..9) is now Ala (G C _) with C at 3rd position due to codon_ref bias ("GCC")
        assert_eq!(&locus.to_consensus_dna().sequence[6..9], b"GCC");
        assert!((locus.provenance.proteomic_fraction - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_ancient_dna_briggs_deamination_correction() {
        // Prior favors 'C' (e.g. from Asian Elephant reference) with 0.82 confidence
        let sites = vec![
            BasePosterior::from_base(b'C', 0.82),
            BasePosterior::from_base(b'G', 0.90),
            BasePosterior::from_base(b'T', 0.90),
        ];
        let mut locus = ReconstructedLocus::new(
            "TRPM8",
            "Mammuthus_primigenius",
            sites,
            LocusProvenance::new_uniform(3, SiteSource::AncestralReconstruction),
        );

        // Read 1 has a 5'-terminal 'T' at offset 0 (classic post-mortem C->T deamination artifact)
        // Read 2 covers offset 0 internally (where deamination is near zero) and observes 'C'
        let read_deaminated = AncientDnaRead {
            base_offset: 0,
            sequence: b"TGT".to_vec(),
            base_accuracy: 0.99,
            deamination_5prime_max: 0.45,
            deamination_3prime_max: 0.25,
        };
        let read_internal = AncientDnaRead {
            base_offset: 0,
            sequence: b"CGT".to_vec(),
            base_accuracy: 0.99,
            deamination_5prime_max: 0.01,
            deamination_3prime_max: 0.01,
        };

        AncientDnaAligner::integrate_reads(&mut locus, &[read_deaminated, read_internal]).unwrap();
        // True ancestral base 'C' is preserved with high confidence
        assert_eq!(locus.sites[0].most_likely_base(), b'C');
        assert!(locus.sites[0].max_prob() > 0.95);
        assert!((locus.provenance.adna_fraction - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_pseudogene_resurrection_avian_enamel() {
        // Chicken ENAM pseudogene has a premature TGA stop at codon 1 and NNN defect at codon 2
        let chicken_pseudo = Dna::new(b"ATGTGANNNGCCTAA").unwrap();
        let alligator_func = ExtantSequence::new(
            "Alligator_mississippiensis",
            Dna::new(b"ATGCGACCCGCTTAA").unwrap(),
        );
        let turtle_func =
            ExtantSequence::new("Chelonia_mydas", Dna::new(b"ATGCGACCCGCTTAA").unwrap());

        let model = SubstitutionModel::default_vertebrate_coding();
        let res = PseudogeneResurrector::resurrect_inactivated_gene(
            "ENAM",
            "Tyrannosaurus_rex",
            &chicken_pseudo,
            68.0,
            &[(alligator_func, 228.0), (turtle_func, 250.0)],
            &model,
        )
        .unwrap();

        assert_eq!(res.repairs.len(), 2);
        assert_eq!(res.repairs[0].codon_index, 1);
        assert_eq!(res.repairs[0].repair_type, RepairType::PrematureStopReversal);
        assert_eq!(res.repairs[0].resurrected_codon, *b"CGA");
        assert_eq!(res.repairs[0].restored_amino_acid, b'R');

        assert_eq!(res.repairs[1].codon_index, 2);
        assert_eq!(res.repairs[1].repair_type, RepairType::FrameshiftDefectRepair);
        assert_eq!(res.repairs[1].resurrected_codon, *b"CCC");
        assert_eq!(res.repairs[1].restored_amino_acid, b'P');

        // Full-length functional protein translates without internal stop codons!
        let protein = res.locus.translate_consensus().unwrap();
        assert_eq!(protein.sequence, b"MRPA");
    }

    #[test]
    fn test_phenotypic_evidence_commutativity_and_correlation_discounting() {
        let osteo1 = OsteologyEvidence {
            specimen_label: "MOR 1125 Femur".to_string(),
            histology: Some(BoneHistologyType::Fibrolamellar),
            has_pneumatic_foramina: true,
            has_enamel_dentition: true,
            estimated_body_mass_kg: Some(6800.0),
            preservation_reliability: 0.90,
        };
        let iso = IsotopeEvidence {
            delta13c_permil: Some(-21.5),
            delta15n_permil: Some(9.4),
            delta18o_permil: Some(17.2),
            reliability: 0.85,
        };

        // Order A -> B vs B -> A must be identical (Commutative)
        let prof_ab = PathwaySelectionProfile::from_evidences(&[&osteo1, &iso]);
        let prof_ba = PathwaySelectionProfile::from_evidences(&[&iso, &osteo1]);

        for pathway in FunctionalPathway::all() {
            let p_ab = prof_ab.posterior_for(pathway);
            let p_ba = prof_ba.posterior_for(pathway);
            assert!(
                (p_ab - p_ba).abs() < 1e-12,
                "Pathway {pathway:?} violated commutativity: {p_ab} vs {p_ba}"
            );
        }

        // Verify biological directions for T. rex (megafaunal endothermic hypercarnivore with enamel teeth)
        assert!(prof_ab.posterior_for(FunctionalPathway::OncogenicSuppressionDnaRepair) > 0.80);
        assert!(prof_ab.posterior_for(FunctionalPathway::EndothermicMetabolism) > 0.85);
        assert!(prof_ab.posterior_for(FunctionalPathway::CarnivoryProteolysis) > 0.80);
        assert!(prof_ab.posterior_for(FunctionalPathway::MineralizedDentitionEnamel) > 0.85);
        assert!(prof_ab.posterior_for(FunctionalPathway::HerbivoryCelluloseDigestion) < 0.25);
        assert!(prof_ab.posterior_for(FunctionalPathway::ColdAdaptationLipidStorage) < 0.25);

        // Verify correlation discounting: 2 identical osteo measurements should NOT double the log-odds
        let prof_single = PathwaySelectionProfile::from_evidences(&[&osteo1]);
        let prof_double = PathwaySelectionProfile::from_evidences(&[&osteo1, &osteo1]);
        let p1 = prof_single.posterior_for(FunctionalPathway::MineralizedDentitionEnamel);
        let p2 = prof_double.posterior_for(FunctionalPathway::MineralizedDentitionEnamel);
        let lo1 = (p1 / (1.0 - p1)).ln();
        let lo2 = (p2 / (1.0 - p2)).ln();
        // With 1/sqrt(2) discounting, 2 identical measurements multiply log-odds by sqrt(2) ~= 1.4142, not 2.0
        assert!(((lo2 / lo1) - std::f64::consts::SQRT_2).abs() < 1e-9);
    }

    #[test]
    fn test_end_to_end_paleo_genome_draft_and_chassis_edits() {
        let mut engine = PaleoReconstructionEngine::new("Tyrannosaurus rex", "Gallus gallus");

        // 1. Add COL1A1 reconstructed locus (differs from Chicken at codon 1: Ala GCC vs Pro CCC)
        let col1a1_sites = vec![
            BasePosterior::from_base(b'A', 0.99),
            BasePosterior::from_base(b'T', 0.99),
            BasePosterior::from_base(b'G', 0.99),
            BasePosterior::from_base(b'G', 0.98), // Non-synonymous C -> G (Pro -> Ala)
            BasePosterior::from_base(b'C', 0.98),
            BasePosterior::from_base(b'C', 0.90),
            BasePosterior::from_base(b'G', 0.99),
            BasePosterior::from_base(b'G', 0.99),
            BasePosterior::from_base(b'A', 0.85), // Synonymous T -> A (Gly -> Gly)
        ];
        let col1a1_locus = ReconstructedLocus::new(
            "COL1A1",
            "Tyrannosaurus rex",
            col1a1_sites,
            LocusProvenance::new_uniform(9, SiteSource::Paleoproteomics),
        );
        let chicken_col1a1 = Dna::new(b"ATGCCCGGT").unwrap();
        engine
            .add_reconstructed_locus(col1a1_locus, chicken_col1a1)
            .unwrap();

        // 2. Add resurrected ENAM tooth enamel pseudogene
        let chicken_enam = Dna::new(b"ATGTGATAA").unwrap(); // Premature stop TGA at codon 1
        let alligator_enam = ExtantSequence::new(
            "Alligator_mississippiensis",
            Dna::new(b"ATGCGATAA").unwrap(),
        );
        let res_enam = PseudogeneResurrector::resurrect_inactivated_gene(
            "ENAM",
            "Tyrannosaurus rex",
            &chicken_enam,
            68.0,
            &[(alligator_enam, 228.0)],
            &SubstitutionModel::default_vertebrate_coding(),
        )
        .unwrap();
        engine
            .add_resurrected_gene(res_enam, chicken_enam)
            .unwrap();

        // 3. Add osteological evidence
        engine.add_phenotypic_evidence(OsteologyEvidence {
            specimen_label: "FMNH PR 2081 (Sue)".to_string(),
            histology: Some(BoneHistologyType::Fibrolamellar),
            has_pneumatic_foramina: true,
            has_enamel_dentition: true,
            estimated_body_mass_kg: Some(8400.0),
            preservation_reliability: 0.95,
        });

        let draft = engine.build_draft().unwrap();
        assert_eq!(draft.summary_metrics.total_chassis_edits, 3);
        assert_eq!(draft.summary_metrics.stop_reversal_edits, 1);
        assert_eq!(draft.summary_metrics.non_synonymous_edits, 1);

        // Verify ranking: PseudogeneStopReversal (ENAM) > NonSynonymous (COL1A1) > Synonymous (COL1A1)
        assert_eq!(
            draft.chassis_edit_plan[0].effect,
            EditEffectType::PseudogeneStopReversal
        );
        assert!(matches!(
            draft.chassis_edit_plan[1].effect,
            EditEffectType::NonSynonymous { .. }
        ));
        assert!(matches!(
            draft.chassis_edit_plan[2].effect,
            EditEffectType::Synonymous { .. }
        ));
    }
}
