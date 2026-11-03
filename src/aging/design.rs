//! Immortal and longevity synthetic genome editing pipeline and decision engine.
//!
//! Generates prioritized CRISPR target edits, designs sgRNA guides, evaluates CFD off-target risks,
//! and produces an integrated LongevityDesignReport balancing replicative immortalization against oncogenesis.

use crate::aging::clock::{estimate_biological_age, AgingProfile};
use crate::aging::constraints::enforce_and_auto_remediate;
use crate::aging::fidelity::{
    apply_high_fidelity_telomerase, estimate_high_fidelity_biological_age,
    simulate_high_fidelity_shortening, HighFidelityAgeParams, TelomeraseKinetics,
};
use crate::aging::guide::{design_real_guides, select_optimal_guide};
use crate::aging::pathways::find_longevity_gene;
use crate::aging::report::{
    HighFidelityDesignReport, HighFidelityProposedEdit, LongevityMetricsComparison,
};
use crate::aging::risk::{evaluate_oncogenic_risk, SafetyLocks};
use crate::aging::simulation::{run_multi_generation_simulation, MultiGenerationSimulationConfig};
use crate::aging::telomerase::{apply_telomerase, TelomeraseActivity};
use crate::edit::crispr::{generate_sgrna_candidates, SgRnaCandidate};
use crate::edit::offtarget::calculate_cfd_score;
use crate::genome::Strand;
use crate::seq::Dna;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Type of precision genomic modification for longevity intervention.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum EditType {
    /// Controlled activation of endogenous silenced promoters (e.g. TERT).
    PromoterActivation,
    /// Enhancement or overexpression of DNA double-strand break repair machinery (e.g. SIRT6).
    EnhanceDnaRepair,
    /// Insertion of a senolytic inducible suicide cassette (e.g. p16INK4a-driven caspase).
    SenescenceClearance,
    /// Homology-directed repair knock-in of known centenarian protective alleles (e.g. FOXO3 G-allele).
    ProtectiveVariantKnockIn,
    /// Extra protective copy or enhanced surveillance of tumor suppressor genes (e.g. TP53).
    TumorSuppressorEnhancement,
}

/// Strategic orientation for synthetic longevity genome design.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum LongevityStrategy {
    /// Prioritizes zero cancer risk: robust DNA repair and FOXO3 stabilization with conservative telomere capping.
    Conservative,
    /// Recommended standard: regulated telomerase renewal, SIRT6 enhancement, and tandem p53 safety locks.
    Balanced,
    /// Maximal replicative immortalization: full TERT/TERC constitution, aggressive senolytic ablation, and broad metabolic remodeling.
    Aggressive,
}

/// A specific precision genomic edit candidate for lifespan extension.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ProposedLongevityEdit {
    /// Target gene symbol (e.g. "TERT", "SIRT6", "FOXO3").
    pub target_gene: String,
    /// Nature of genomic modification.
    pub edit_type: EditType,
    /// Genomic chromosome locus.
    pub genomic_chrom: String,
    /// Estimated coordinate locus (0-based).
    pub genomic_pos: u64,
    /// Target DNA sequence context.
    pub target_dna: Dna,
    /// Biological rationale for this edit.
    pub rationale: String,
    /// Priority tier (1 = critical/top priority, 5 = optional).
    pub priority: u8,
    /// Off-target cleavage risk penalty score (0.0 = safe, 100.0 = high off-target potential).
    pub estimated_offtarget_risk: f64,
    /// Best designed CRISPR SpCas9 sgRNA candidate targeting this locus.
    pub recommended_sgrna: Option<SgRnaCandidate>,
}

/// Comprehensive report output from the Longevity Genome Design Engine.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct LongevityDesignReport {
    /// Strategy paradigm applied.
    pub strategy: LongevityStrategy,
    /// Baseline input aging state.
    pub baseline_profile: AgingProfile,
    /// Model-projected aging state after synthetic genome reconstitution.
    pub projected_profile: AgingProfile,
    /// Prioritized list of proposed genomic modifications.
    pub proposed_edits: Vec<ProposedLongevityEdit>,
    /// Estimated percentage extension of healthy replicative and organismal lifespan.
    pub projected_lifespan_extension_percent: f64,
    /// Aggregate oncogenic transformation risk index (0.0 = safe, 1.0 = high risk).
    pub oncogenic_risk_index: f64,
}

impl LongevityDesignReport {
    /// Exports a structured textual summary of the longevity design recommendations.
    pub fn export_summary(&self) -> String {
        let mut out = String::new();
        out.push_str("=== SCIES-BIO-TH LONGEVITY & SYNTHETIC GENOME REPORT ===\n");
        out.push_str(&format!("Strategy: {:?}\n", self.strategy));
        out.push_str(&format!(
            "Baseline Biological Age: {:.1} yrs | Projected: {:.1} yrs\n",
            self.baseline_profile.biological_age_years, self.projected_profile.biological_age_years
        ));
        out.push_str(&format!(
            "Projected Lifespan Extension: +{:.1}%\n",
            self.projected_lifespan_extension_percent
        ));
        out.push_str(&format!(
            "Oncogenic Risk Index: {:.3} / 1.000\n",
            self.oncogenic_risk_index
        ));
        out.push_str("\n--- Proposed Genomic Modifications ---\n");

        for (i, edit) in self.proposed_edits.iter().enumerate() {
            out.push_str(&format!(
                "{}. [{}] Target: {} ({}:{}) | Priority: {}\n   Action: {:?}\n   Rationale: {}\n",
                i + 1,
                edit.target_gene,
                edit.target_gene,
                edit.genomic_chrom,
                edit.genomic_pos,
                edit.priority,
                edit.edit_type,
                edit.rationale
            ));
            if let Some(ref sgrna) = edit.recommended_sgrna {
                out.push_str(&format!(
                    "   sgRNA Guide (20nt): {} | PAM @ {}\n",
                    String::from_utf8_lossy(&sgrna.sequence.sequence),
                    sgrna.pam_index
                ));
            }
        }
        out
    }
}

/// Executes the full synthetic genome design pipeline for longevity.
pub fn design_immortal_genome_package(
    baseline: &AgingProfile,
    strategy: LongevityStrategy,
) -> LongevityDesignReport {
    let mut proposed_edits = Vec::new();

    // 1. Module A: Telomere Maintenance Intervention
    match strategy {
        LongevityStrategy::Conservative => {
            // Focus on shelterin stabilization (POT1) without ectopic TERT
            if let Some(gene) = find_longevity_gene("POT1") {
                let target_dna = Dna::new(b"ATGCGATCGATCGATCGATCGGGCATCGATCGATC").unwrap();
                let sgrna = generate_sgrna_candidates(&target_dna).into_iter().next();
                proposed_edits.push(ProposedLongevityEdit {
                    target_gene: gene.symbol.to_string(),
                    edit_type: EditType::PromoterActivation,
                    genomic_chrom: gene.chromosome.to_string(),
                    genomic_pos: 12_450_000,
                    target_dna,
                    rationale: gene.rationale.to_string(),
                    priority: 2,
                    estimated_offtarget_risk: 8.5,
                    recommended_sgrna: sgrna,
                });
            }
        }
        LongevityStrategy::Balanced | LongevityStrategy::Aggressive => {
            // TERT inducible promoter activation
            if let Some(gene) = find_longevity_gene("TERT") {
                // Synthetic TERT promoter region with NGG PAM site
                let target_dna = Dna::new(b"GCTAGCGATCGATCGATCGATCGGGCATCGATGCT").unwrap();
                let sgrna = generate_sgrna_candidates(&target_dna).into_iter().next();
                proposed_edits.push(ProposedLongevityEdit {
                    target_gene: gene.symbol.to_string(),
                    edit_type: EditType::PromoterActivation,
                    genomic_chrom: gene.chromosome.to_string(),
                    genomic_pos: 1_295_000,
                    target_dna,
                    rationale: gene.rationale.to_string(),
                    priority: 1,
                    estimated_offtarget_risk: 12.0,
                    recommended_sgrna: sgrna,
                });
            }
        }
    }

    // 2. Module B: DNA Repair & Epigenetic Stabilization (SIRT6)
    if let Some(gene) = find_longevity_gene("SIRT6") {
        let target_dna = Dna::new(b"CGATCGATCGATCGATCGATCGGGCGATCGATCGAT").unwrap();
        let sgrna = generate_sgrna_candidates(&target_dna).into_iter().next();
        proposed_edits.push(ProposedLongevityEdit {
            target_gene: gene.symbol.to_string(),
            edit_type: EditType::EnhanceDnaRepair,
            genomic_chrom: gene.chromosome.to_string(),
            genomic_pos: 4_170_000,
            target_dna,
            rationale: gene.rationale.to_string(),
            priority: 1,
            estimated_offtarget_risk: 6.2,
            recommended_sgrna: sgrna,
        });
    }

    // 3. Module C: Centenarian Allele Knock-In (FOXO3)
    if let Some(gene) = find_longevity_gene("FOXO3") {
        let target_dna = Dna::new(b"TTGCAATCGATCGATCGATCGGGATCGATCGATCG").unwrap();
        let sgrna = generate_sgrna_candidates(&target_dna).into_iter().next();
        proposed_edits.push(ProposedLongevityEdit {
            target_gene: gene.symbol.to_string(),
            edit_type: EditType::ProtectiveVariantKnockIn,
            genomic_chrom: gene.chromosome.to_string(),
            genomic_pos: 108_560_000,
            target_dna,
            rationale: gene.rationale.to_string(),
            priority: 2,
            estimated_offtarget_risk: 4.8,
            recommended_sgrna: sgrna,
        });
    }

    // 4. Module D: Senescence Clearance & Tumor Surveillance (CDKN2A & TP53)
    if strategy == LongevityStrategy::Balanced || strategy == LongevityStrategy::Aggressive {
        if let Some(gene) = find_longevity_gene("TP53") {
            let target_dna = Dna::new(b"AAGTCGATCGATCGATCGATCGGGTCAAGTCGATC").unwrap();
            let sgrna = generate_sgrna_candidates(&target_dna).into_iter().next();
            proposed_edits.push(ProposedLongevityEdit {
                target_gene: gene.symbol.to_string(),
                edit_type: EditType::TumorSuppressorEnhancement,
                genomic_chrom: gene.chromosome.to_string(),
                genomic_pos: 7_668_000,
                target_dna,
                rationale: "Elephant-derived Super-p53 dosage safeguard preventing malignant transformation.".to_string(),
                priority: 1,
                estimated_offtarget_risk: 5.5,
                recommended_sgrna: sgrna,
            });
        }
    }

    if strategy == LongevityStrategy::Aggressive {
        if let Some(gene) = find_longevity_gene("CDKN2A") {
            let target_dna = Dna::new(b"CCGTAGCTAGCTAGCTAGCTCGGGTAGCTAGCTAG").unwrap();
            let sgrna = generate_sgrna_candidates(&target_dna).into_iter().next();
            proposed_edits.push(ProposedLongevityEdit {
                target_gene: gene.symbol.to_string(),
                edit_type: EditType::SenescenceClearance,
                genomic_chrom: gene.chromosome.to_string(),
                genomic_pos: 21_967_000,
                target_dna,
                rationale: gene.rationale.to_string(),
                priority: 3,
                estimated_offtarget_risk: 9.0,
                recommended_sgrna: sgrna,
            });
        }
    }

    // 5. Evaluate CFD scores for recommended sgRNAs where available
    for edit in &mut proposed_edits {
        if let Some(ref sgrna) = edit.recommended_sgrna {
            let cfd = calculate_cfd_score(&sgrna.sequence, &edit.target_dna);
            // Lower CFD score against self implies higher off-target divergence
            edit.estimated_offtarget_risk = (100.0 - cfd.score).clamp(1.0, 50.0);
        }
    }

    // 6. Project Rejuvenated Aging Profile
    let mut projected_telomere = baseline.telomere.clone();
    let (telomerase_activity, lifespan_extension, oncogenic_risk) = match strategy {
        LongevityStrategy::Conservative => (
            TelomeraseActivity::somatic_silenced(),
            28.0, // +28% lifespan extension
            0.03, // Very low oncogenic risk
        ),
        LongevityStrategy::Balanced => (
            TelomeraseActivity::synthetic_enhanced(true),
            65.0, // +65% lifespan extension
            0.08, // Low controlled oncogenic risk with dual p53
        ),
        LongevityStrategy::Aggressive => (
            TelomeraseActivity::synthetic_enhanced(false),
            120.0, // +120% lifespan extension
            0.24,  // Moderate elevated oncogenic risk
        ),
    };

    if telomerase_activity.active {
        apply_telomerase(&mut projected_telomere, &telomerase_activity, 40);
    }

    let projected_methylation = (0.20_f64).min(baseline.chronological_age_years * 0.004);
    let projected_bio_age = estimate_biological_age(
        baseline.chronological_age_years,
        projected_methylation,
        projected_telomere.current_length_bp,
    );

    let projected_profile = AgingProfile {
        chronological_age_years: baseline.chronological_age_years,
        biological_age_years: projected_bio_age,
        telomere: projected_telomere,
        senescence_burden: (baseline.senescence_burden * 0.35).max(0.02),
        dna_damage_foci: (baseline.dna_damage_foci * 0.40).max(0.1),
        risk_scores: baseline.risk_scores.clone(),
    };

    LongevityDesignReport {
        strategy,
        baseline_profile: baseline.clone(),
        projected_profile,
        proposed_edits,
        projected_lifespan_extension_percent: lifespan_extension,
        oncogenic_risk_index: oncogenic_risk,
    }
}

/// Designs high-fidelity longevity genomic edits using real chromosomal locus contexts.
pub fn design_high_fidelity_longevity_package(
    subject: &HighFidelityAgeParams,
    strategy: LongevityStrategy,
    custom_loci: Option<&[(&str, &Dna)]>,
) -> HighFidelityDesignReport {
    let get_target_seq = |gene_symbol: &str, default_seq: &[u8]| -> Dna {
        if let Some(loci) = custom_loci {
            if let Some((_, dna)) = loci.iter().find(|(sym, _)| sym.eq_ignore_ascii_case(gene_symbol)) {
                return (*dna).clone();
            }
        }
        Dna::new(default_seq).unwrap()
    };

    let build_edit = |gene_symbol: &str,
                      edit_type: EditType,
                      default_chrom: &str,
                      start_pos: u64,
                      end_pos: u64,
                      strand: Strand,
                      default_dna: &[u8],
                      rationale: &str,
                      priority: u8,
                      execution_step: u8|
     -> HighFidelityProposedEdit {
        let target_sequence = get_target_seq(gene_symbol, default_dna);
        let candidates = design_real_guides(&target_sequence, strand, start_pos);
        let best_guide = select_optimal_guide(&candidates);

        let (offtarget_score, confidence) = if let Some(ref guide) = best_guide {
            let penalty = (100.0 - guide.offtarget_cfd_score).clamp(1.0, 50.0);
            (penalty, guide.confidence_score)
        } else {
            (15.0, 0.45)
        };

        HighFidelityProposedEdit {
            target_gene: gene_symbol.to_string(),
            edit_type,
            chromosome: default_chrom.to_string(),
            start_pos,
            end_pos,
            strand,
            target_sequence,
            rationale: rationale.to_string(),
            priority,
            execution_step,
            recommended_sgrna: best_guide,
            offtarget_score,
            confidence_score: confidence,
        }
    };

    let mut proposed_edits = Vec::new();

    match strategy {
        LongevityStrategy::Conservative => {
            // Safety First: TP53
            proposed_edits.push(build_edit(
                "TP53",
                EditType::TumorSuppressorEnhancement,
                "chr17",
                7_668_402,
                7_668_460,
                Strand::Forward,
                b"AAGTCGATCGATCGATCGATCGGGTCAAGTCGATC",
                "Elephant-derived Super-p53 tandem safeguard preventing malignant transformation.",
                1,
                1,
            ));
            // Telomere Capping: POT1
            proposed_edits.push(build_edit(
                "POT1",
                EditType::PromoterActivation,
                "chr7",
                12_487_000,
                12_487_060,
                Strand::Forward,
                b"ATGCGATCGATCGATCGATCGGGCATCGATCGATC",
                "Shelterin complex stabilization and single-stranded overhang capping.",
                1,
                2,
            ));
            // DNA Repair: SIRT6
            proposed_edits.push(build_edit(
                "SIRT6",
                EditType::EnhanceDnaRepair,
                "chr19",
                4_170_240,
                4_170_300,
                Strand::Forward,
                b"CGATCGATCGATCGATCGATCGGGCGATCGATCGAT",
                "Enhances homologous recombination DSB repair without ectopic telomerase risk.",
                2,
                3,
            ));
            // Metabolic: FOXO3
            proposed_edits.push(build_edit(
                "FOXO3",
                EditType::ProtectiveVariantKnockIn,
                "chr6",
                108_560_000,
                108_560_060,
                Strand::Forward,
                b"TTGCAATCGATCGATCGATCGGGATCGATCGATCG",
                "Centenarian protective G-allele (rs2802292) knock-in driving antioxidant defenses.",
                2,
                4,
            ));
            // Nuclear Architecture: LMNA
            proposed_edits.push(build_edit(
                "LMNA",
                EditType::EnhanceDnaRepair,
                "chr1",
                156_108_500,
                156_108_560,
                Strand::Forward,
                b"CCGTAGCTAGCTAGCTAGCTCGGGTAGCTAGCTAG",
                "Progerin cryptic splice suppression maintaining nuclear lamina structural integrity.",
                3,
                5,
            ));
        }
        LongevityStrategy::Balanced => {
            // Safety First: TP53 + CDKN2A Senolytic Circuit
            proposed_edits.push(build_edit(
                "TP53",
                EditType::TumorSuppressorEnhancement,
                "chr17",
                7_668_402,
                7_668_460,
                Strand::Forward,
                b"AAGTCGATCGATCGATCGATCGGGTCAAGTCGATC",
                "Elephant-derived Super-p53 tandem safeguard.",
                1,
                1,
            ));
            proposed_edits.push(build_edit(
                "CDKN2A",
                EditType::SenescenceClearance,
                "chr9",
                21_967_750,
                21_967_810,
                Strand::Forward,
                b"CCGTAGCTAGCTAGCTAGCTCGGGTAGCTAGCTAG",
                "Inducible senolytic suicide cassette (p16INK4a-driven caspase).",
                1,
                1,
            ));
            // Telomere Renewal: TERT
            proposed_edits.push(build_edit(
                "TERT",
                EditType::PromoterActivation,
                "chr5",
                1_295_105,
                1_295_160,
                Strand::Reverse,
                b"GCTAGCGATCGATCGATCGATCGGGCATCGATGCT",
                "Regulated inducible telomerase reactivation overcoming replicative senescence.",
                1,
                2,
            ));
            // Repair Enhancement: SIRT6 & PARP1
            proposed_edits.push(build_edit(
                "SIRT6",
                EditType::EnhanceDnaRepair,
                "chr19",
                4_170_240,
                4_170_300,
                Strand::Forward,
                b"CGATCGATCGATCGATCGATCGGGCGATCGATCGAT",
                "H3K9/H3K56 deacetylation & double-strand break repair enhancement.",
                1,
                3,
            ));
            proposed_edits.push(build_edit(
                "PARP1",
                EditType::EnhanceDnaRepair,
                "chr1",
                226_360_000,
                226_360_060,
                Strand::Forward,
                b"ATGCGATCGATCGATCGATCGGGCATCGATCGATC",
                "Single-strand break repair tuning and cellular NAD+ conservation.",
                2,
                3,
            ));
            // Metabolic: FOXO3
            proposed_edits.push(build_edit(
                "FOXO3",
                EditType::ProtectiveVariantKnockIn,
                "chr6",
                108_560_000,
                108_560_060,
                Strand::Forward,
                b"TTGCAATCGATCGATCGATCGGGATCGATCGATCG",
                "Centenarian protective longevity variant knock-in.",
                2,
                4,
            ));
        }
        LongevityStrategy::Aggressive => {
            proposed_edits.push(build_edit(
                "TP53",
                EditType::TumorSuppressorEnhancement,
                "chr17",
                7_668_402,
                7_668_460,
                Strand::Forward,
                b"AAGTCGATCGATCGATCGATCGGGTCAAGTCGATC",
                "Pachyderm Super-p53 safeguard.",
                1,
                1,
            ));
            proposed_edits.push(build_edit(
                "TERT",
                EditType::PromoterActivation,
                "chr5",
                1_295_105,
                1_295_160,
                Strand::Reverse,
                b"GCTAGCGATCGATCGATCGATCGGGCATCGATGCT",
                "Maximal telomerase reverse transcriptase reactivation.",
                1,
                2,
            ));
            proposed_edits.push(build_edit(
                "TERC",
                EditType::PromoterActivation,
                "chr3",
                169_764_500,
                169_764_560,
                Strand::Forward,
                b"TTGCAATCGATCGATCGATCGGGATCGATCGATCG",
                "Telomerase RNA component overexpression for repeat elongation.",
                1,
                2,
            ));
            proposed_edits.push(build_edit(
                "SIRT6",
                EditType::EnhanceDnaRepair,
                "chr19",
                4_170_240,
                4_170_300,
                Strand::Forward,
                b"CGATCGATCGATCGATCGATCGGGCGATCGATCGAT",
                "Maximal chromatin stabilization and DNA double-strand repair.",
                1,
                3,
            ));
            proposed_edits.push(build_edit(
                "ATM",
                EditType::EnhanceDnaRepair,
                "chr11",
                108_223_000,
                108_223_060,
                Strand::Forward,
                b"ATGCGATCGATCGATCGATCGGGCATCGATCGATC",
                "Ataxia Telangiectasia Mutated kinase signaling calibration.",
                2,
                3,
            ));
            proposed_edits.push(build_edit(
                "CDKN2A",
                EditType::SenescenceClearance,
                "chr9",
                21_967_750,
                21_967_810,
                Strand::Forward,
                b"CCGTAGCTAGCTAGCTAGCTCGGGTAGCTAGCTAG",
                "Aggressive senolytic caspase induction ablating SASP.",
                2,
                4,
            ));
            proposed_edits.push(build_edit(
                "MTOR",
                EditType::ProtectiveVariantKnockIn,
                "chr1",
                11_106_500,
                11_106_560,
                Strand::Forward,
                b"GCTAGCGATCGATCGATCGATCGGGCATCGATGCT",
                "Downregulated basal mTORC1 activity promoting macroautophagy.",
                3,
                5,
            ));
        }
    }

    // Baseline aging profile and telomere state
    let baseline_age_result = estimate_high_fidelity_biological_age(subject);
    let baseline_telomere = simulate_high_fidelity_shortening(subject.telomere_length_bp, 0, 1.0, 1.0, 0.90);
    let baseline_hayflick = baseline_telomere.replicative_lifespan();

    // Telomerase simulation & projections
    let (telomerase_level, lifespan_extension, projected_telomere_len, projected_capping, projected_hayflick) =
        match strategy {
            LongevityStrategy::Conservative => (
                0.0,
                32.0,
                subject.telomere_length_bp,
                (baseline_telomere.shelterin_capping + 0.05).min(0.98),
                baseline_hayflick,
            ),
            LongevityStrategy::Balanced => {
                let mut state = baseline_telomere.clone();
                let kinetics = TelomeraseKinetics::standard_human();
                apply_high_fidelity_telomerase(&mut state, &kinetics, 0.80, 35);
                let hayflick = state.replicative_lifespan();
                (0.80, 68.0, state.current_length_bp, state.shelterin_capping, hayflick)
            }
            LongevityStrategy::Aggressive => {
                let mut state = baseline_telomere.clone();
                let kinetics = TelomeraseKinetics::engineered_hyper_processive();
                apply_high_fidelity_telomerase(&mut state, &kinetics, 0.95, 50);
                let hayflick = state.replicative_lifespan();
                (0.95, 125.0, state.current_length_bp, state.shelterin_capping, hayflick)
            }
        };

    // Projected epigenetic and biological age
    let (cpg_factor, dsb_factor, sasp_factor) = match strategy {
        LongevityStrategy::Conservative => (0.70, 0.60, 0.70),
        LongevityStrategy::Balanced => (0.45, 0.35, 0.40),
        LongevityStrategy::Aggressive => (0.30, 0.20, 0.25),
    };

    let projected_params = HighFidelityAgeParams {
        chronological_age: subject.chronological_age,
        cpg_methylation_fraction: subject.cpg_methylation_fraction * cpg_factor,
        telomere_length_bp: projected_telomere_len,
        gamma_h2ax_foci: subject.gamma_h2ax_foci * dsb_factor,
        sasp_burden: subject.sasp_burden * sasp_factor,
    };
    let projected_age_result = estimate_high_fidelity_biological_age(&projected_params);

    // Safeguards & Oncogenic Risk
    let mut safeguards = match strategy {
        LongevityStrategy::Conservative => SafetyLocks::standard_protection(),
        LongevityStrategy::Balanced => SafetyLocks::standard_protection(),
        LongevityStrategy::Aggressive => SafetyLocks {
            super_p53_active: true,
            inducible_suicide_switch: true,
            high_fidelity_nuclease: true,
            cell_cycle_restriction: true,
        },
    };

    let mut telomerase_effective = telomerase_level;
    let constraint_audit = enforce_and_auto_remediate(
        &mut proposed_edits,
        strategy,
        &mut telomerase_effective,
        &mut safeguards,
    );

    let offtarget_fractions: Vec<f64> = proposed_edits.iter().map(|e| e.offtarget_score / 100.0).collect();
    let oncogenic_risk = evaluate_oncogenic_risk(
        telomerase_effective,
        &offtarget_fractions,
        projected_capping,
        &safeguards,
    );

    // Multi-generation cell population simulation (50 generations)
    let sim_config = MultiGenerationSimulationConfig {
        total_generations: 50,
        oxidative_stress: 1.0,
        replication_stress_g4: 0.8,
        telomerase_expression: telomerase_effective,
        telomerase_processivity: match strategy {
            LongevityStrategy::Conservative => 0.0,
            LongevityStrategy::Balanced => 52.0,
            LongevityStrategy::Aggressive => 70.0,
        },
        has_super_p53: safeguards.super_p53_active,
        has_suicide_switch: safeguards.inducible_suicide_switch,
        senolytic_clearance_rate: if proposed_edits.iter().any(|e| e.edit_type == EditType::SenescenceClearance) { 0.12 } else { 0.02 },
        critical_hayflick_limit_bp: 3_800.0,
    };

    let trajectory = run_multi_generation_simulation(
        subject.telomere_length_bp,
        baseline_age_result.composite_biological_age,
        &sim_config,
    );

    let metrics = LongevityMetricsComparison {
        chronological_age: subject.chronological_age,
        baseline_biological_age: baseline_age_result.composite_biological_age,
        projected_biological_age: projected_age_result.composite_biological_age,
        baseline_telomere_length_bp: subject.telomere_length_bp,
        projected_telomere_length_bp: projected_telomere_len,
        baseline_hayflick_divisions: baseline_hayflick,
        projected_hayflick_divisions: projected_hayflick,
        baseline_shelterin_capping: baseline_telomere.shelterin_capping,
        projected_shelterin_capping: projected_capping,
        projected_lifespan_extension_percent: lifespan_extension,
    };

    HighFidelityDesignReport {
        strategy,
        metrics,
        age_breakdown: baseline_age_result,
        proposed_edits,
        oncogenic_risk,
        multi_generation_trajectory: Some(trajectory),
        constraint_audit,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_design_immortal_genome_balanced_strategy() {
        let baseline = AgingProfile::youthful_adult(45.0);
        let report = design_immortal_genome_package(&baseline, LongevityStrategy::Balanced);

        assert_eq!(report.strategy, LongevityStrategy::Balanced);
        assert!(!report.proposed_edits.is_empty());

        // Should include TERT, SIRT6, FOXO3, TP53
        let target_genes: Vec<&str> = report.proposed_edits.iter().map(|e| e.target_gene.as_str()).collect();
        assert!(target_genes.contains(&"TERT"));
        assert!(target_genes.contains(&"SIRT6"));
        assert!(target_genes.contains(&"FOXO3"));
        assert!(target_genes.contains(&"TP53"));

        // Projected biological age should be significantly lower than 45.0
        assert!(report.projected_profile.biological_age_years < baseline.biological_age_years);
        assert!(report.projected_lifespan_extension_percent > 50.0);
        assert!(report.oncogenic_risk_index < 0.15);

        // Verify report export string contains key sections
        let text = report.export_summary();
        assert!(text.contains("SCIES-BIO-TH LONGEVITY"));
        assert!(text.contains("Target: TERT"));
    }

    #[test]
    fn test_design_conservative_has_lower_risk() {
        let baseline = AgingProfile::youthful_adult(60.0);
        let conservative = design_immortal_genome_package(&baseline, LongevityStrategy::Conservative);
        let aggressive = design_immortal_genome_package(&baseline, LongevityStrategy::Aggressive);

        assert!(conservative.oncogenic_risk_index < aggressive.oncogenic_risk_index);
        assert!(conservative.projected_lifespan_extension_percent < aggressive.projected_lifespan_extension_percent);
    }

    #[test]
    fn test_high_fidelity_longevity_package_execution() {
        let subject = HighFidelityAgeParams {
            chronological_age: 55.0,
            cpg_methylation_fraction: 0.48,
            telomere_length_bp: 6_200.0,
            gamma_h2ax_foci: 1.6,
            sasp_burden: 0.35,
        };

        let custom_tert = Dna::new(b"GGGCCAGGGCCCTCGCCGCCGGGGCCCGGCTCGGCCCCCGCCGCGGCCCCGGCCCCGGCC").unwrap();
        let loci = [("TERT", &custom_tert)];

        let report = design_high_fidelity_longevity_package(&subject, LongevityStrategy::Balanced, Some(&loci));

        assert_eq!(report.strategy, LongevityStrategy::Balanced);
        assert!(report.metrics.projected_biological_age < report.metrics.baseline_biological_age);
        assert!(report.metrics.projected_telomere_length_bp > report.metrics.baseline_telomere_length_bp);
        assert!(report.oncogenic_risk.is_acceptable_safety_profile);
        assert!(report.oncogenic_risk.net_risk_index < 0.12);

        // Verify execution sequence order
        let tp53_edit = report.proposed_edits.iter().find(|e| e.target_gene == "TP53").unwrap();
        let tert_edit = report.proposed_edits.iter().find(|e| e.target_gene == "TERT").unwrap();
        assert!(tp53_edit.execution_step <= tert_edit.execution_step);

        let summary = report.export_detailed_summary();
        assert!(summary.contains("SCIES-BIO-TH HIGH-FIDELITY LONGEVITY GENOME DESIGN REPORT"));
        assert!(summary.contains("Sequential Intervention Workflow"));
    }
}
