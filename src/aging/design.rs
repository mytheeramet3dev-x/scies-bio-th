//! Immortal and longevity synthetic genome editing pipeline and decision engine.
//!
//! Generates prioritized CRISPR target edits, designs sgRNA guides, evaluates CFD off-target risks,
//! and produces an integrated LongevityDesignReport balancing replicative immortalization against oncogenesis.

use crate::aging::clock::{estimate_biological_age, AgingProfile};
use crate::aging::pathways::find_longevity_gene;
use crate::aging::telomerase::{apply_telomerase, TelomeraseActivity};
use crate::edit::crispr::{generate_sgrna_candidates, SgRnaCandidate};
use crate::edit::offtarget::calculate_cfd_score;
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
}
