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
    SequenceProvenance,
};
use crate::aging::risk::{evaluate_oncogenic_risk, SafetyLocks};
use crate::aging::simulation::{run_multi_generation_simulation, MultiGenerationSimulationConfig};
use crate::aging::telomerase::{apply_telomerase, TelomeraseActivity};
use crate::aging::telomere::predict_replicative_lifespan;
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
    /// Estimated stem cell pool fraction.
    pub stem_pool_fraction: f64,
    /// Somatic lineage identity retention fraction.
    pub identity_retention: f64,
    /// Recommended autonomous guard action.
    pub recommended_guard_action: crate::safety::GuardAction,
}

impl LongevityDesignReport {
    /// Extracts the standardized intervention safety coupling payload.
    pub fn safety_coupling(&self) -> crate::safety::InterventionSafetyCoupling {
        crate::safety::InterventionSafetyCoupling {
            oncogenic_risk: self.oncogenic_risk_index,
            stem_pool_fraction: self.stem_pool_fraction,
            identity_retention: self.identity_retention,
            recommended_guard_action: self.recommended_guard_action,
        }
    }

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

/// Computes projected healthy lifespan extension (%) from physiological age rejuvenation,
/// replicative Hayflick division gain, DNA repair/senolytic pathway coverage, and net oncogenic risk.
///
/// Replaces static strategy constants so the projection responds dynamically to subject baseline
/// parameters, telomerase kinetics, and off-target hazard penalties.
fn compute_projected_lifespan_extension_percent(
    baseline_bio_age: f64,
    projected_bio_age: f64,
    baseline_hayflick: u32,
    projected_hayflick: u32,
    dna_repair_and_protective_edits: usize,
    senolytic_active: bool,
    net_oncogenic_risk: f64,
) -> f64 {
    // Normalized against reference 80-year human healthspan so absolute years of age reversal scale proportionally
    let age_rejuvenation_pct = ((baseline_bio_age - projected_bio_age).max(0.0) / 80.0) * 100.0;
    // Normalized against canonical Hayflick & Moorhead (1961) human diploid fibroblast limit (60 divisions)
    let hayflick_gain_pct =
        ((projected_hayflick.saturating_sub(baseline_hayflick)) as f64 / 60.0) * 65.0;
    let pathway_resilience_pct = (dna_repair_and_protective_edits as f64 * 6.5)
        + if senolytic_active { 8.0 } else { 0.0 };
    let oncogenic_penalty_pct = net_oncogenic_risk * 25.0;

    (0.95 * age_rejuvenation_pct + 0.65 * hayflick_gain_pct + pathway_resilience_pct
        - oncogenic_penalty_pct)
        .clamp(0.0, 200.0)
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
                    genomic_pos: 124_929_500,
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
                    genomic_pos: 1_295_104,
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
            genomic_pos: 4_182_500,
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
            genomic_pos: 108_559_652,
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
                genomic_pos: 7_687_400,
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
                genomic_pos: 21_974_800,
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

    // 6. Project Rejuvenated Aging Profile & Compute Risk + Lifespan Extension Dynamically
    let mut projected_telomere = baseline.telomere.clone();
    let (telomerase_activity, shelterin_capping, safeguards) = match strategy {
        LongevityStrategy::Conservative => (
            TelomeraseActivity::somatic_silenced(),
            0.94,
            SafetyLocks {
                super_p53_active: false,
                inducible_suicide_switch: false,
                high_fidelity_nuclease: true,
                cell_cycle_restriction: true,
            },
        ),
        LongevityStrategy::Balanced => (
            TelomeraseActivity::synthetic_enhanced(true),
            0.88,
            SafetyLocks::standard_protection(),
        ),
        LongevityStrategy::Aggressive => (
            TelomeraseActivity::synthetic_enhanced(false),
            0.74,
            SafetyLocks {
                super_p53_active: proposed_edits.iter().any(|e| e.target_gene == "TP53"),
                inducible_suicide_switch: false,
                high_fidelity_nuclease: false,
                cell_cycle_restriction: false,
            },
        ),
    };

    if telomerase_activity.active {
        let cycles = if telomerase_activity.safety_lock_present { 40 } else { 55 };
        apply_telomerase(&mut projected_telomere, &telomerase_activity, cycles);
    }

    let projected_methylation = (0.20_f64).min(baseline.chronological_age_years * 0.004);
    let projected_bio_age = estimate_biological_age(
        baseline.chronological_age_years,
        projected_methylation,
        projected_telomere.current_length_bp,
    );

    let offtarget_fractions: Vec<f64> = proposed_edits
        .iter()
        .map(|e| e.estimated_offtarget_risk / 100.0)
        .collect();
    let telomerase_expr = if telomerase_activity.active {
        telomerase_activity.expression_level
    } else {
        0.0
    };
    let oncogenic_risk = evaluate_oncogenic_risk(
        telomerase_expr,
        &offtarget_fractions,
        shelterin_capping,
        &safeguards,
    )
    .net_risk_index;

    let repair_edits = proposed_edits
        .iter()
        .filter(|e| {
            matches!(
                e.edit_type,
                EditType::EnhanceDnaRepair | EditType::ProtectiveVariantKnockIn
            )
        })
        .count();
    let senolytic_active = proposed_edits
        .iter()
        .any(|e| e.edit_type == EditType::SenescenceClearance);
    let lifespan_extension = compute_projected_lifespan_extension_percent(
        baseline.biological_age_years,
        projected_bio_age,
        predict_replicative_lifespan(&baseline.telomere),
        predict_replicative_lifespan(&projected_telomere),
        repair_edits,
        senolytic_active,
        oncogenic_risk,
    );

    let projected_profile = AgingProfile {
        chronological_age_years: baseline.chronological_age_years,
        biological_age_years: projected_bio_age,
        telomere: projected_telomere,
        senescence_burden: (baseline.senescence_burden * 0.35).max(0.02),
        dna_damage_foci: (baseline.dna_damage_foci * 0.40).max(0.1),
        risk_scores: baseline.risk_scores.clone(),
    };

    let stem_pool_fraction = match strategy {
        LongevityStrategy::Conservative => 0.05,
        LongevityStrategy::Balanced => 0.08,
        LongevityStrategy::Aggressive => 0.12,
    };
    let identity_retention = match strategy {
        LongevityStrategy::Conservative => 0.98,
        LongevityStrategy::Balanced => 0.95,
        LongevityStrategy::Aggressive => 0.90,
    };
    let recommended_guard_action = if oncogenic_risk >= crate::aging::MAX_SAFE_ONCOGENIC_RISK_THRESHOLD {
        crate::safety::GuardAction::AbortSuicideSwitch
    } else if oncogenic_risk > 0.08 {
        crate::safety::GuardAction::TitrateDown
    } else {
        crate::safety::GuardAction::Continue
    };

    LongevityDesignReport {
        strategy,
        baseline_profile: baseline.clone(),
        projected_profile,
        proposed_edits,
        projected_lifespan_extension_percent: lifespan_extension,
        oncogenic_risk_index: oncogenic_risk,
        stem_pool_fraction,
        identity_retention,
        recommended_guard_action,
    }
}

/// Designs high-fidelity longevity genomic edits using real chromosomal locus contexts.
pub fn design_high_fidelity_longevity_package(
    subject: &HighFidelityAgeParams,
    strategy: LongevityStrategy,
    custom_loci: Option<&[(&str, &Dna)]>,
) -> HighFidelityDesignReport {
    let get_target_seq = |gene_symbol: &str, default_seq: &[u8]| -> (Dna, SequenceProvenance) {
        if let Some(loci) = custom_loci {
            if let Some((_, dna)) = loci.iter().find(|(sym, _)| sym.eq_ignore_ascii_case(gene_symbol)) {
                return ((*dna).clone(), SequenceProvenance::RealReference);
            }
        }
        (
            Dna::new(default_seq).unwrap(),
            SequenceProvenance::SyntheticPlaceholder,
        )
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
        let (target_sequence, sequence_provenance) = get_target_seq(gene_symbol, default_dna);
        let candidates = design_real_guides(&target_sequence, strand, start_pos);
        let best_guide = select_optimal_guide(&candidates);

        let (offtarget_score, confidence) = if let Some(ref guide) = best_guide {
            let penalty = (100.0 - guide.offtarget_cfd_score).clamp(1.0, 50.0);
            (penalty, guide.confidence_score)
        } else {
            (15.0, 0.45)
        };

        HighFidelityProposedEdit {
            sequence_provenance,
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
                7_687_001,
                7_688_000,
                Strand::Reverse,
                b"AAGTCGATCGATCGATCGATCGGGTCAAGTCGATC",
                "Elephant-derived Super-p53 tandem safeguard preventing malignant transformation.",
                1,
                1,
            ));
            // Telomere Capping: POT1 (GRCh38 chr7:124929401-124930200)
            proposed_edits.push(build_edit(
                "POT1",
                EditType::PromoterActivation,
                "chr7",
                124_929_401,
                124_930_200,
                Strand::Reverse,
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
                4_182_001,
                4_182_800,
                Strand::Reverse,
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
                108_559_501,
                108_560_300,
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
                156_138_201,
                156_139_000,
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
                7_687_001,
                7_688_000,
                Strand::Reverse,
                b"AAGTCGATCGATCGATCGATCGGGTCAAGTCGATC",
                "Elephant-derived Super-p53 tandem safeguard.",
                1,
                1,
            ));
            proposed_edits.push(build_edit(
                "CDKN2A",
                EditType::SenescenceClearance,
                "chr9",
                21_974_401,
                21_975_400,
                Strand::Reverse,
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
                1_294_801,
                1_295_600,
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
                4_182_001,
                4_182_800,
                Strand::Reverse,
                b"CGATCGATCGATCGATCGATCGGGCGATCGATCGAT",
                "H3K9/H3K56 deacetylation & double-strand break repair enhancement.",
                1,
                3,
            ));
            proposed_edits.push(build_edit(
                "PARP1",
                EditType::EnhanceDnaRepair,
                "chr1",
                226_407_601,
                226_408_400,
                Strand::Reverse,
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
                108_559_501,
                108_560_300,
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
                7_687_001,
                7_688_000,
                Strand::Reverse,
                b"AAGTCGATCGATCGATCGATCGGGTCAAGTCGATC",
                "Pachyderm Super-p53 safeguard.",
                1,
                1,
            ));
            proposed_edits.push(build_edit(
                "TERT",
                EditType::PromoterActivation,
                "chr5",
                1_294_801,
                1_295_600,
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
                169_764_401,
                169_765_200,
                Strand::Reverse,
                b"TTGCAATCGATCGATCGATCGGGATCGATCGATCG",
                "Telomerase RNA component overexpression for repeat elongation.",
                1,
                2,
            ));
            proposed_edits.push(build_edit(
                "SIRT6",
                EditType::EnhanceDnaRepair,
                "chr19",
                4_182_001,
                4_182_800,
                Strand::Reverse,
                b"CGATCGATCGATCGATCGATCGGGCGATCGATCGAT",
                "Maximal chromatin stabilization and DNA double-strand repair.",
                1,
                3,
            ));
            proposed_edits.push(build_edit(
                "ATM",
                EditType::EnhanceDnaRepair,
                "chr11",
                108_222_601,
                108_223_400,
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
                21_974_401,
                21_975_400,
                Strand::Reverse,
                b"CCGTAGCTAGCTAGCTAGCTCGGGTAGCTAGCTAG",
                "Aggressive senolytic caspase induction ablating SASP.",
                2,
                4,
            ));
            proposed_edits.push(build_edit(
                "MTOR",
                EditType::ProtectiveVariantKnockIn,
                "chr1",
                11_262_001,
                11_262_800,
                Strand::Reverse,
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
    let (telomerase_level, projected_telomere_len, projected_capping, projected_hayflick) =
        match strategy {
            LongevityStrategy::Conservative => (
                0.0,
                subject.telomere_length_bp,
                (baseline_telomere.shelterin_capping + 0.05).min(0.98),
                baseline_hayflick,
            ),
            LongevityStrategy::Balanced => {
                let mut state = baseline_telomere.clone();
                let kinetics = TelomeraseKinetics::standard_human();
                apply_high_fidelity_telomerase(&mut state, &kinetics, 0.80, 35);
                let hayflick = state.replicative_lifespan();
                (0.80, state.current_length_bp, state.shelterin_capping, hayflick)
            }
            LongevityStrategy::Aggressive => {
                let mut state = baseline_telomere.clone();
                let kinetics = TelomeraseKinetics::engineered_hyper_processive();
                apply_high_fidelity_telomerase(&mut state, &kinetics, 0.95, 50);
                let hayflick = state.replicative_lifespan();
                (0.95, state.current_length_bp, state.shelterin_capping, hayflick)
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

    let repair_edits = proposed_edits
        .iter()
        .filter(|e| {
            matches!(
                e.edit_type,
                EditType::EnhanceDnaRepair | EditType::ProtectiveVariantKnockIn
            )
        })
        .count();
    let senolytic_active = proposed_edits
        .iter()
        .any(|e| e.edit_type == EditType::SenescenceClearance);
    let lifespan_extension = compute_projected_lifespan_extension_percent(
        baseline_age_result.composite_biological_age,
        projected_age_result.composite_biological_age,
        baseline_hayflick,
        projected_hayflick,
        repair_edits,
        senolytic_active,
        oncogenic_risk.net_risk_index,
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

    let stem_pool_fraction = match strategy {
        LongevityStrategy::Conservative => 0.05,
        LongevityStrategy::Balanced => 0.08,
        LongevityStrategy::Aggressive => 0.12,
    };
    let identity_retention = match strategy {
        LongevityStrategy::Conservative => 0.98,
        LongevityStrategy::Balanced => 0.95,
        LongevityStrategy::Aggressive => 0.90,
    };
    let recommended_guard_action = if oncogenic_risk.net_risk_index >= crate::aging::MAX_SAFE_ONCOGENIC_RISK_THRESHOLD {
        crate::safety::GuardAction::AbortSuicideSwitch
    } else if oncogenic_risk.net_risk_index > 0.08 {
        crate::safety::GuardAction::TitrateDown
    } else {
        crate::safety::GuardAction::Continue
    };

    HighFidelityDesignReport {
        strategy,
        metrics,
        age_breakdown: baseline_age_result,
        proposed_edits,
        oncogenic_risk,
        stem_pool_fraction,
        identity_retention,
        recommended_guard_action,
        multi_generation_trajectory: Some(trajectory),
        constraint_audit,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn load_testdata_locus(file: &str) -> Dna {
        let path = format!("{}/testdata/{}", env!("CARGO_MANIFEST_DIR"), file);
        let text = std::fs::read_to_string(&path).expect("testdata FASTA present");
        let seq: Vec<u8> = text
            .lines()
            .filter(|l| !l.starts_with('>'))
            .flat_map(|l| l.trim().bytes())
            .map(|b| b.to_ascii_uppercase())
            .filter(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
            .collect();
        Dna::new(&seq).expect("valid ACGT locus")
    }

    fn subject() -> HighFidelityAgeParams {
        HighFidelityAgeParams {
            chronological_age: 55.0,
            cpg_methylation_fraction: 0.48,
            telomere_length_bp: 6_200.0,
            gamma_h2ax_foci: 1.6,
            sasp_burden: 0.35,
        }
    }

    #[test]
    fn test_provenance_is_flagged_when_no_real_loci_supplied() {
        let report = design_high_fidelity_longevity_package(&subject(), LongevityStrategy::Balanced, None);
        assert!(report
            .proposed_edits
            .iter()
            .all(|e| e.sequence_provenance == SequenceProvenance::SyntheticPlaceholder));
        let summary = report.export_detailed_summary();
        assert!(summary.contains("SYNTHETIC PLACEHOLDER"));
        assert!(!summary.contains("REAL reference locus"));
    }

    #[test]
    fn test_real_fasta_loci_are_marked_real_and_only_for_supplied_genes() {
        let sirt6 = load_testdata_locus("t2t_chr19_sirt6.fasta");
        let tp53 = load_testdata_locus("t2t_chr17_tp53.fasta");
        assert!(sirt6.sequence.len() > 500, "SIRT6 locus should be real-length, not a toy string");
        let loci = [("SIRT6", &sirt6), ("TP53", &tp53)];

        let report = design_high_fidelity_longevity_package(
            &subject(),
            LongevityStrategy::Balanced,
            Some(&loci),
        );
        for edit in &report.proposed_edits {
            let expected = if edit.target_gene == "SIRT6" || edit.target_gene == "TP53" {
                SequenceProvenance::RealReference
            } else {
                SequenceProvenance::SyntheticPlaceholder
            };
            assert_eq!(edit.sequence_provenance, expected, "gene {}", edit.target_gene);
        }
        let sirt6_edit = report.proposed_edits.iter().find(|e| e.target_gene == "SIRT6").unwrap();
        assert_eq!(sirt6_edit.target_sequence, sirt6);
    }

    #[test]
    fn test_offtarget_score_is_computed_not_constant() {
        // Different real loci must yield different guides/scores; a hardcoded value would be identical.
        let sirt6 = load_testdata_locus("t2t_chr19_sirt6.fasta");
        let tp53 = load_testdata_locus("t2t_chr17_tp53.fasta");
        let loci = [("SIRT6", &sirt6), ("TP53", &tp53)];
        let report = design_high_fidelity_longevity_package(
            &subject(),
            LongevityStrategy::Balanced,
            Some(&loci),
        );
        let s = report.proposed_edits.iter().find(|e| e.target_gene == "SIRT6").unwrap();
        let t = report.proposed_edits.iter().find(|e| e.target_gene == "TP53").unwrap();
        assert!(s.recommended_sgrna.is_some() && t.recommended_sgrna.is_some());
        assert_ne!(
            s.recommended_sgrna.as_ref().unwrap().sequence,
            t.recommended_sgrna.as_ref().unwrap().sequence
        );
    }

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

    #[test]
    fn test_lifespan_extension_and_oncogenic_risk_respond_to_subject_state() {
        let young_long_telomere = HighFidelityAgeParams {
            chronological_age: 30.0,
            cpg_methylation_fraction: 0.24,
            telomere_length_bp: 11_200.0,
            gamma_h2ax_foci: 0.4,
            sasp_burden: 0.08,
        };
        let aged_depleted_telomere = HighFidelityAgeParams {
            chronological_age: 72.0,
            cpg_methylation_fraction: 0.68,
            telomere_length_bp: 4_600.0,
            gamma_h2ax_foci: 2.9,
            sasp_burden: 0.58,
        };

        let rep_young = design_high_fidelity_longevity_package(
            &young_long_telomere,
            LongevityStrategy::Balanced,
            None,
        );
        let rep_aged = design_high_fidelity_longevity_package(
            &aged_depleted_telomere,
            LongevityStrategy::Balanced,
            None,
        );

        // Lifespan extension % and net oncogenic risk must NOT be hardcoded constants across subjects
        assert!(
            (rep_young.metrics.projected_lifespan_extension_percent
                - rep_aged.metrics.projected_lifespan_extension_percent)
                .abs()
                > 10.0,
            "Lifespan extension should respond to subject baseline state (young={:.2}%, aged={:.2}%)",
            rep_young.metrics.projected_lifespan_extension_percent,
            rep_aged.metrics.projected_lifespan_extension_percent
        );
        assert!(
            (rep_young.oncogenic_risk.net_risk_index - rep_aged.oncogenic_risk.net_risk_index)
                .abs()
                > 1e-4,
            "Oncogenic risk should vary with shelterin capping and telomere integrity"
        );
    }
}
