//! Reconstituted Ovarian Organoid (`rOvary`) and 3-Stage In Vitro Oogenesis Engine
//! (Hayashi et al. 2016; Saitou et al. 2018/2023).
//!
//! Models the somatic-germline niche interactions between Primordial Germ Cell-Like Cells (PGCLCs)
//! and Fetal Ovarian Somatic Cell-Like Cells (FOSLCs / Granulosa & Theca lineages) across:
//! 1. **IVD (In Vitro Differentiation)**: Primordial $\rightarrow$ Secondary follicle assembly,
//!    paracrine `GDF9/BMP15` $\leftrightarrow$ `KITL (SCF)` signaling, cytoplasmic growth, and mtDNA amplification.
//! 2. **IVG (In Vitro Growth)**: FSH-dependent granulosa proliferation, antral cavity formation,
//!    and de novo Maternal Imprinting Control Region (ICR) methylation by DNMT3A/3L.
//! 3. **IVM (In Vitro Maturation)**: LH + EGF/Amphiregulin surge driving Germinal Vesicle Breakdown (GVBD),
//!    Meiosis I spindle assembly, first polar body extrusion, and Metaphase II (MII) arrest.

use crate::error::{BioError, Result};
use crate::ivg::epigenetics::{EpigeneticProfile, ImprintPattern};
use crate::ivg::gamete::{
    derive_homologous_partner, recombine_homologous_chromosomes,
    simulate_karyotype_meiotic_crossovers, species_bivalent_pairs, GameteQuality, GameteType,
    PgclcCulture, SyntheticGamete, ZonaPellucida,
};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Developmental stage of an ovarian follicle inside a reconstituted ovarian organoid (`rOvary`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum FollicleStage {
    /// Dormant or newly assembled primordial follicle (oocyte ~20 µm, flat pre-granulosa).
    Primordial,
    /// Primary follicle (cuboidal granulosa monolayer, Zona Pellucida assembly begins).
    Primary,
    /// Secondary pre-antral follicle (multilayered granulosa + theca, rapid cytoplasmic growth).
    Secondary,
    /// Antral Graafian follicle (fluid-filled antrum, cumulus-oocyte complex formed, ~110 µm).
    Antral,
    /// Mature Metaphase II (MII) oocyte with extruded first polar body, ready for IVF/ICSI.
    PreOvulatoryMII,
}

/// Culture and hormonal parameters for a 3-stage `rOvary` organoid oogenesis pipeline.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct OvarianOrganoidConfig {
    /// Ratio of supporting FOSLC (granulosa/theca) somatic cells per PGCLC (optimal ~5:1 to 10:1).
    pub foslc_to_pgclc_ratio: f64,
    /// Oocyte-secreted GDF9 + BMP15 paracrine signaling activity in $[0.0, 1.0]$.
    pub gdf9_bmp15_signaling: f64,
    /// Granulosa-secreted KIT Ligand (SCF) trophic support in $[0.0, 1.0]$.
    pub kitl_scf_support: f64,
    /// Follicle-Stimulating Hormone (FSH) concentration during IVG phase (mIU/mL, optimal ~10–25).
    pub fsh_miu_ml: f64,
    /// Luteinizing Hormone (LH / hCG) + EGF surge during IVM phase (mIU/mL, optimal ~50–100).
    pub lh_surge_miu_ml: f64,
    /// Duration of In Vitro Differentiation (IVD) phase in days (optimal ~21 days).
    pub ivd_days: u32,
    /// Duration of In Vitro Growth (IVG) phase in days (optimal ~11 days).
    pub ivg_days: u32,
    /// Duration of In Vitro Maturation (IVM) phase in hours (optimal ~16–24 hours).
    pub ivm_hours: u32,
}

impl OvarianOrganoidConfig {
    /// Standard high-competence mammalian/human reconstituted ovarian organoid protocol.
    pub fn standard_rovary() -> Self {
        Self {
            foslc_to_pgclc_ratio: 7.5,
            gdf9_bmp15_signaling: 0.92,
            kitl_scf_support: 0.94,
            fsh_miu_ml: 15.0,
            lh_surge_miu_ml: 75.0,
            ivd_days: 21,
            ivg_days: 11,
            ivm_hours: 18,
        }
    }
}

/// Comprehensive biophysical and developmental report from `simulate_rovary_folliculogenesis`.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct OogenesisOrganoidReport {
    /// Final follicular maturation stage achieved.
    pub final_stage: FollicleStage,
    /// Final oocyte diameter in micrometers ($\mu\text{m}$, mature mammalian MII is $105\text{--}120\ \mu\text{m}$).
    pub oocyte_diameter_um: f64,
    /// Number of homologous meiotic chiasmata crossovers formed during Prophase I.
    pub meiotic_chiasmata_count: usize,
    /// Whether Germinal Vesicle Breakdown (GVBD) succeeded during IVM.
    pub gvbd_completed: bool,
    /// Completeness of DNMT3A/3L maternal imprinting (`SNRPN`, `IGF2R`, `KCNQ1OT1`) in $[0.0, 1.0]$.
    pub maternal_icr_methylation: f64,
    /// Resulting mature synthetic oocyte construct.
    pub oocyte: SyntheticGamete,
}

/// Simulates the complete 3-stage (`IVD -> IVG -> IVM`) reconstituted ovarian organoid (`rOvary`)
/// folliculogenesis pipeline from a [`PgclcCulture`] to a Metaphase II [`SyntheticGamete`].
pub fn simulate_rovary_folliculogenesis(
    pgclc: &PgclcCulture,
    config: &OvarianOrganoidConfig,
    meiotic_seed: u64,
) -> Result<OogenesisOrganoidReport> {
    if pgclc.induction_efficiency < 0.30 {
        return Err(BioError::ReprogrammingFailed(
            "PGCLC induction efficiency (< 0.30) is insufficient to seed rOvary organoids"
                .to_string(),
        ));
    }
    if config.foslc_to_pgclc_ratio < 1.0 {
        return Err(BioError::InvalidInput(
            "rOvary organoid requires at least 1:1 supporting FOSLC granulosa cells per PGCLC"
                .to_string(),
        ));
    }

    // 1. Stage 1 — IVD (Primordial -> Secondary follicle & cytoplasmic expansion)
    let niche_synergy = (0.50 * config.gdf9_bmp15_signaling.clamp(0.0, 1.0)
        + 0.50 * config.kitl_scf_support.clamp(0.0, 1.0))
        * (config.foslc_to_pgclc_ratio / 6.0).clamp(0.4, 1.15);
    let ivd_progress = ((config.ivd_days as f64) / 21.0).clamp(0.1, 1.2) * niche_synergy;

    // 2. Homologous Meiotic Recombination (Prophase I pachytene chiasmata across full karyotype)
    let homolog_partner = derive_homologous_partner(&pgclc.donor_dna, meiotic_seed);
    let (haploid_dna, _contig_crossovers) =
        recombine_homologous_chromosomes(&pgclc.donor_dna, &homolog_partner, meiotic_seed)?;
    let meiotic_competence =
        (0.55 * pgclc.induction_efficiency + 0.45 * niche_synergy.min(1.0)).clamp(0.1, 1.0);
    let (chiasmata, achiasmate_bivalents) = simulate_karyotype_meiotic_crossovers(
        pgclc.tax_id,
        meiotic_competence,
        meiotic_seed,
    );
    let expected_bivalents = species_bivalent_pairs(pgclc.tax_id);

    // 3. Stage 2 — IVG (FSH-driven antral cavity & DNMT3A/3L maternal ICR establishment)
    let fsh_response = (config.fsh_miu_ml / 15.0).clamp(0.1, 1.2);
    let ivg_progress = ((config.ivg_days as f64) / 11.0).clamp(0.1, 1.2) * fsh_response;

    let oocyte_diameter_um = (20.0 + 55.0 * ivd_progress.min(1.0) + 40.0 * ivg_progress.min(1.0))
        .clamp(20.0, 125.0);

    // Mitochondrial biogenesis scales with cytoplasmic volume from ~2,000 in PGCLC to ~120,000 in MII oocyte
    let vol_ratio = ((oocyte_diameter_um / 115.0).powi(3)).clamp(0.02, 1.2);
    let mitochondrial_count = (115_000.0 * vol_ratio).round() as u64;

    // Size-dependent de novo maternal imprinting (requires oocyte diameter >= 70 um)
    let maternal_icr_methylation = if oocyte_diameter_um >= 70.0 {
        (pgclc.induction_efficiency * 0.55 + 0.45 * ivg_progress.min(1.0)).clamp(0.0, 0.99)
    } else {
        (0.30 * ivd_progress).clamp(0.0, 0.50)
    };

    // 4. Stage 3 — IVM (LH surge -> GVBD -> MI spindle -> First Polar Body extrusion -> MII arrest)
    let lh_competent = config.lh_surge_miu_ml >= 25.0 && config.ivm_hours >= 12;
    let gvbd_completed = oocyte_diameter_um >= 80.0 && lh_competent;
    let polar_body_extruded = gvbd_completed && chiasmata >= expected_bivalents;

    let final_stage = if polar_body_extruded {
        FollicleStage::PreOvulatoryMII
    } else if oocyte_diameter_um >= 85.0 {
        FollicleStage::Antral
    } else if oocyte_diameter_um >= 55.0 {
        FollicleStage::Secondary
    } else if oocyte_diameter_um >= 35.0 {
        FollicleStage::Primary
    } else {
        FollicleStage::Primordial
    };

    // Obligatory chiasmata across all bivalents and proper niche signaling suppress Meiosis I non-disjunction
    // Reference: Hassold & Hunt (2001) Nat Rev Genet 2:280; Broman et al. (1998) AJHG 63:861.
    let target_chiasmata = expected_bivalents * 2;
    let crossover_deficit_penalty = if chiasmata < target_chiasmata {
        ((target_chiasmata - chiasmata) as f64 / target_chiasmata as f64) * 0.35
    } else {
        0.0
    };
    let achiasmate_penalty = achiasmate_bivalents as f64 * 0.18;
    let aneuploidy_risk = (0.04
        + achiasmate_penalty
        + crossover_deficit_penalty
        + (1.0 - niche_synergy.min(1.0)) * 0.15)
        .clamp(0.01, 0.95);

    let overall_score = (maternal_icr_methylation * 45.0
        + (1.0 - aneuploidy_risk) * 35.0
        + (if polar_body_extruded { 20.0 } else { 0.0 }))
        .clamp(0.0, 100.0);

    let mut ep_profile = EpigeneticProfile::mature_oocyte();
    ep_profile.imprint_pattern = if maternal_icr_methylation >= 0.75 {
        ImprintPattern::Maternal
    } else {
        ImprintPattern::Aberrant
    };

    let zona = if pgclc.tax_id == 9606 {
        Some(ZonaPellucida::human_default())
    } else {
        Some(ZonaPellucida::cave_lion_default())
    };

    let oocyte = SyntheticGamete {
        gamete_type: GameteType::Oocyte,
        haploid_genome: haploid_dna,
        epigenetic_profile: ep_profile,
        zona_pellucida: zona,
        quality: GameteQuality {
            aneuploidy_risk,
            imprint_fidelity: maternal_icr_methylation,
            mitochondrial_count,
            polar_body_extruded,
            overall_score,
        },
        tax_id: pgclc.tax_id,
    };

    Ok(OogenesisOrganoidReport {
        final_stage,
        oocyte_diameter_um,
        meiotic_chiasmata_count: chiasmata,
        gvbd_completed,
        maternal_icr_methylation,
        oocyte,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ivg::gamete::InductionProtocol;
    use crate::seq::Dna;

    #[test]
    fn test_rovary_folliculogenesis_to_mii_oocyte() {
        let donor_dna = Dna::new(b"ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATC").unwrap();
        let pgclc = PgclcCulture::induce(&donor_dna, 9606, InductionProtocol::Hayashi2016, 6).unwrap();
        let cfg = OvarianOrganoidConfig::standard_rovary();

        let report = simulate_rovary_folliculogenesis(&pgclc, &cfg, 2026).unwrap();
        assert_eq!(report.final_stage, FollicleStage::PreOvulatoryMII);
        assert!(report.oocyte_diameter_um >= 100.0);
        assert!(report.gvbd_completed);
        assert!(
            report.meiotic_chiasmata_count >= 45,
            "Human MII oocyte should have >=45 chiasmata across 23 bivalents, got {}",
            report.meiotic_chiasmata_count
        );
        assert_ne!(
            report.oocyte.haploid_genome, donor_dna,
            "Recombinant haploid genome should differ from parental donor_dna"
        );
        assert!(report.maternal_icr_methylation >= 0.85);
        assert!(report.oocyte.quality.polar_body_extruded);
        assert!(report.oocyte.quality.mitochondrial_count >= 80_000);
        assert!(report.oocyte.quality.overall_score >= 85.0);
    }

    #[test]
    fn test_suboptimal_meiotic_niche_penalizes_crossovers_and_aneuploidy() {
        let donor_dna = Dna::new(b"ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATC").unwrap();
        let mut poor_pgclc =
            PgclcCulture::induce(&donor_dna, 9606, InductionProtocol::ChemicallyDefinedFeederFree, 2)
                .unwrap();
        poor_pgclc.induction_efficiency = 0.35; // Barely above 0.30 survival floor
        let mut poor_cfg = OvarianOrganoidConfig::standard_rovary();
        poor_cfg.gdf9_bmp15_signaling = 0.35;
        poor_cfg.kitl_scf_support = 0.35;

        let report = simulate_rovary_folliculogenesis(&poor_pgclc, &poor_cfg, 2026).unwrap();
        assert!(report.meiotic_chiasmata_count < 23);
        assert!(report.oocyte.quality.aneuploidy_risk > 0.35);
        assert!(report.oocyte.quality.overall_score < 60.0);
    }
}
