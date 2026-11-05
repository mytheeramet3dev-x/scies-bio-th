//! Demonstration of the Tri-Layer Paleo-Genomic & Paleoproteomic Reconstruction Engine (`scies_bio_th::paleo`).
//!
//! Showcases two scientifically grounded reconstruction scenarios:
//! 1. **Tier 1 (Late Pleistocene aDNA)**: *Mammuthus primigenius* (Woolly Mammoth) `TRPM8` cold-sensor locus
//!    using Asian/African Elephant bracketing + Briggs deamination-aware ancient DNA reads + permafrost isotopes.
//! 2. **Tier 2/3 (Mesozoic Deep-Time Proteomics + Pseudogene Resurrection)**: *Tyrannosaurus rex* `COL1A1`
//!    collagen peptide back-translation + `ENAM` tooth-enamel pseudogene resurrection using the Avian-Crocodilian
//!    Extant Phylogenetic Bracket.

use scies_bio_th::paleo::*;
use scies_bio_th::phylo::PhyloNode;
use scies_bio_th::seq::Dna;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("================================================================================");
    println!(" 🦖 scies-bio-th v1.0-paleo: Tri-Layer Paleo-Genomic Reconstruction Demonstration");
    println!("================================================================================\n");

    run_mammoth_pleistocene_demo()?;
    println!();
    run_theropod_deep_time_demo()?;

    Ok(())
}

fn run_mammoth_pleistocene_demo() -> Result<(), Box<dyn std::error::Error>> {
    println!("--------------------------------------------------------------------------------");
    println!(" [Case 1] Tier 1 (Pleistocene ~40 ka): Woolly Mammoth (Mammuthus primigenius)");
    println!("          Chassis: Asian Elephant (Elephas maximus) | Target Gene: TRPM8");
    println!("--------------------------------------------------------------------------------");

    // 1. Phylogenetic tree: ((Elephas_maximus:6.0, Loxodonta_africana:7.5):18.0, Trichechus_manatus:55.0)
    let elephantids = PhyloNode::internal(
        vec![
            PhyloNode::leaf("Elephas_maximus", 6.0),
            PhyloNode::leaf("Loxodonta_africana", 7.5),
        ],
        18.0,
    );
    let sirenia = PhyloNode::leaf("Trichechus_manatus", 55.0);
    let tree = PhyloNode::internal(vec![elephantids, sirenia], 0.0);

    let elephas_trpm8 = Dna::new(b"ATGGCCCGCTTCCTGGGCACCGAG")?;
    let loxodonta_trpm8 = Dna::new(b"ATGGCCCGCTTCCTGGGCACCGAA")?;
    let manatee_trpm8 = Dna::new(b"ATGGCTCGCTTCCTGGGCACCGAG")?;

    let extant = vec![
        ExtantSequence::new("Elephas_maximus", elephas_trpm8.clone()),
        ExtantSequence::new("Loxodonta_africana", loxodonta_trpm8),
        ExtantSequence::new("Trichechus_manatus", manatee_trpm8),
    ];

    let model = SubstitutionModel::default_vertebrate_coding();
    let mut trpm8_locus = FelsensteinAsr::reconstruct_locus(
        &tree,
        &extant,
        &model,
        "Mammuthus primigenius",
        "TRPM8",
        4.5,
    )?;

    // 2. Integrate Pleistocene permafrost aDNA reads (showing Mammoth-specific cold-adaptation missense mutation C->A at pos 4)
    // Read 1 also has a post-mortem C->T deamination artifact at its 5' base (pos 0: 'C' would deaminate, here pos 3 'G'->'A' at 3')
    let reads = vec![
        AncientDnaRead::new_pleistocene(0, b"ATGACCCGCTTCCTGGGCACCGAG"),
        AncientDnaRead::new_pleistocene(0, b"ATGACCCGCTTCCTGGGCACCGAG"),
    ];
    AncientDnaAligner::integrate_reads(&mut trpm8_locus, &reads)?;

    let mut engine = PaleoReconstructionEngine::new("Mammuthus primigenius", "Elephas maximus");
    engine.add_reconstructed_locus(trpm8_locus, elephas_trpm8)?;

    engine.add_phenotypic_evidence(OsteologyEvidence {
        specimen_label: "Yukagir Mammoth Femur".to_string(),
        histology: Some(BoneHistologyType::Fibrolamellar),
        has_pneumatic_foramina: false,
        has_enamel_dentition: true,
        estimated_body_mass_kg: Some(5500.0),
        preservation_reliability: 0.96,
    });
    engine.add_phenotypic_evidence(IsotopeEvidence {
        delta13c_permil: Some(-21.2),
        delta15n_permil: Some(4.8),
        delta18o_permil: Some(-18.4), // Extreme Siberian glacial water signal
        reliability: 0.94,
    });

    let draft = engine.build_draft()?;
    println!(
        "  • Expected Base Accuracy : {:.2}% (Mean Entropy: {:.4} bits)",
        draft.summary_metrics.overall_base_accuracy_expected * 100.0,
        draft.summary_metrics.mean_shannon_entropy_bits
    );
    println!(
        "  • Cold Adaptation Pathway: {:.1}% posterior selection probability",
        draft
            .selection_profile
            .posterior_for(FunctionalPathway::ColdAdaptationLipidStorage)
            * 100.0
    );
    println!(
        "  • Peto's Paradox (TP53)  : {:.1}% posterior selection probability",
        draft
            .selection_profile
            .posterior_for(FunctionalPathway::OncogenicSuppressionDnaRepair)
            * 100.0
    );
    println!(
        "  • Chassis Edits Required : {} edit(s) against Elephas maximus",
        draft.chassis_edit_plan.len()
    );
    for edit in &draft.chassis_edit_plan {
        println!(
            "    -> [{}] pos {}: {} -> {} ({:?}, conf={:.1}%, priority={:.2})",
            edit.gene_symbol,
            edit.locus_position,
            edit.chassis_base as char,
            edit.target_base as char,
            edit.effect,
            edit.posterior_confidence * 100.0,
            edit.priority_score
        );
    }

    Ok(())
}

fn run_theropod_deep_time_demo() -> Result<(), Box<dyn std::error::Error>> {
    println!("--------------------------------------------------------------------------------");
    println!(" [Case 2] Tier 2/3 (Cretaceous 68 Ma): Tyrannosaurus rex");
    println!("          Chassis: Chicken (Gallus gallus) | Loci: COL1A1 (Proteome) & ENAM (Pseudogene)");
    println!("--------------------------------------------------------------------------------");

    // 1. Archosauria Phylogenetic Tree: ((Gallus_gallus:68, Falco_peregrinus:68):160, Alligator_mississippiensis:228)
    let aves = PhyloNode::internal(
        vec![
            PhyloNode::leaf("Gallus_gallus", 68.0),
            PhyloNode::leaf("Falco_peregrinus", 68.0),
        ],
        160.0,
    );
    let croc = PhyloNode::leaf("Alligator_mississippiensis", 228.0);
    let tree = PhyloNode::internal(vec![aves, croc], 0.0);

    let chicken_col1a1 = Dna::new(b"GGTCCCCCTGGTCCTGCTGGCCCTCCT")?;
    let falcon_col1a1 = Dna::new(b"GGTCCCCCTGGTCCTGCTGGCCCTCCA")?;
    let alligator_col1a1 = Dna::new(b"GGTCCTCCTGGTCCTGCCGGCCCCCCT")?;

    let extant_col1a1 = vec![
        ExtantSequence::new("Gallus_gallus", chicken_col1a1.clone()),
        ExtantSequence::new("Falco_peregrinus", falcon_col1a1),
        ExtantSequence::new("Alligator_mississippiensis", alligator_col1a1.clone()),
    ];

    let model = SubstitutionModel::default_vertebrate_coding();
    let mut col1a1_locus = FelsensteinAsr::reconstruct_locus(
        &tree,
        &extant_col1a1,
        &model,
        "Tyrannosaurus rex",
        "COL1A1",
        66.0,
    )?;

    // 2. Anchor COL1A1 with MOR 1125 fossil collagen peptide ("GPPGPAGPP", using 'O' for hydroxyproline)
    let fossil_peptide = PeptideFragment::new(0, b"GOOGPAGOO", 0.96);
    PaleoproteomeTranslator::anchor_with_peptides(
        &mut col1a1_locus,
        &[fossil_peptide],
        Some(&chicken_col1a1),
    )?;

    // 3. Resurrect inactivated avian tooth-enamel pseudogene (ENAM) using Crocodilian + Testudine bracket
    // Chicken ENAM has a premature TGA stop at codon 2 and TAG stop at codon 5
    let chicken_enam_pseudo = Dna::new(b"ATGCCCTGAGGCAGCTAGGGTCCCTAA")?;
    let alligator_enam_func = ExtantSequence::new(
        "Alligator_mississippiensis",
        Dna::new(b"ATGCCCCGAGGCAGCCAGGGTCCCTAA")?,
    );
    let turtle_enam_func = ExtantSequence::new(
        "Chelonia_mydas",
        Dna::new(b"ATGCCCCGAGGCAGCCAGGGTCCCTAA")?,
    );

    let resurrected_enam = PseudogeneResurrector::resurrect_inactivated_gene(
        "ENAM",
        "Tyrannosaurus rex",
        &chicken_enam_pseudo,
        68.0,
        &[(alligator_enam_func, 228.0), (turtle_enam_func, 250.0)],
        &model,
    )?;

    let mut engine = PaleoReconstructionEngine::new("Tyrannosaurus rex", "Gallus gallus");
    engine.add_reconstructed_locus(col1a1_locus, chicken_col1a1)?;
    engine.add_resurrected_gene(resurrected_enam, chicken_enam_pseudo)?;

    engine.add_phenotypic_evidence(OsteologyEvidence {
        specimen_label: "MOR 1125 / FMNH PR 2081".to_string(),
        histology: Some(BoneHistologyType::Fibrolamellar),
        has_pneumatic_foramina: true,
        has_enamel_dentition: true,
        estimated_body_mass_kg: Some(7800.0),
        preservation_reliability: 0.92,
    });
    engine.add_phenotypic_evidence(IsotopeEvidence {
        delta13c_permil: Some(-22.1),
        delta15n_permil: Some(9.2),
        delta18o_permil: Some(17.8),
        reliability: 0.88,
    });

    let draft = engine.build_draft()?;
    println!(
        "  • Expected Base Accuracy : {:.2}% (Proteomic: {:.0}%, Pseudogene Repair: {:.0}%)",
        draft.summary_metrics.overall_base_accuracy_expected * 100.0,
        draft.summary_metrics.fraction_anchored_by_proteomics * 100.0,
        draft.summary_metrics.fraction_pseudogene_repaired * 100.0,
    );
    println!(
        "  • Enamel Dentition Prior : {:.1}% posterior selection probability",
        draft
            .selection_profile
            .posterior_for(FunctionalPathway::MineralizedDentitionEnamel)
            * 100.0
    );
    println!(
        "  • Resurrected ENAM Events: {} premature stop codons reversed -> translated protein: {}",
        draft.resurrected_genes[0].repairs.len(),
        String::from_utf8_lossy(
            &draft.resurrected_genes[0]
                .locus
                .translate_consensus()?
                .sequence
        )
    );
    println!(
        "  • Ranked Chassis Edits   : {} total ({} stop reversals, {} missense)",
        draft.summary_metrics.total_chassis_edits,
        draft.summary_metrics.stop_reversal_edits,
        draft.summary_metrics.non_synonymous_edits
    );
    for edit in &draft.chassis_edit_plan {
        println!(
            "    -> [{}] pos {:2}: {} -> {} ({:?}, conf={:.1}%, priority={:.2})",
            edit.gene_symbol,
            edit.locus_position,
            edit.chassis_base as char,
            edit.target_base as char,
            edit.effect,
            edit.posterior_confidence * 100.0,
            edit.priority_score
        );
    }

    Ok(())
}
