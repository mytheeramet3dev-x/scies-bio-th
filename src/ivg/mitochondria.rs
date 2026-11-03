//! Mitochondrial DNA (`mtDNA`) heteroplasmy, 4,977-bp deletion burden, reactive oxygen species (`ROS`),
//! and purifying selection (`PINK1`/`Parkin` mitophagy) in somatic rejuvenation and In Vitro Gametogenesis.
//!
//! References:
//! - Cortopassi & Arnheim (1990) *Nucleic Acids Res* 18:6927 (Age-dependent 4,977-bp mtDNA common deletion).
//! - Taylor & Turnbull (2005) *Nat Rev Genet* 6:389 (Clinical OXPHOS biochemical threshold ~60%).
//! - Wai et al. (2008) *Nat Genet* 40:1484 (Postnatal germline genetic bottleneck of mammalian mtDNA).
//! - Narendra et al. (2008) *J Cell Biol* 183:795 (Parkin is recruited selectively to impaired mitochondria).
//! - Lieber et al. (2019) *Nature* 570:380 (Selective removal of deleterious mtDNA via mitophagy).

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Bioenergetic, mutational, and oxidative state of a cellular mitochondrial population.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct MitochondrialPopulationState {
    /// Total combined mutant mtDNA heteroplasmy fraction ($0.0$ to $1.0$).
    pub heteroplasmy_fraction: f64,
    /// Total mtDNA copy number per cell (~2,200 in somatic cells; ~115,000 in mature MII oocytes).
    pub mtdna_copy_number: u64,
    /// Mitochondrial Reactive Oxygen Species (`ROS`: superoxide $\text{O}_2^{\bullet-}$ / $\text{H}_2\text{O}_2$) level ($0.0$ to $1.0$).
    pub ros_level: f64,
    /// Alias of `mtdna_copy_number` preserved for germline bottleneck compatibility.
    pub total_copy_number: u64,
    /// Fraction of mtDNA molecules carrying pathogenic OXPHOS point mutations ($0.0$ to $1.0$).
    pub point_mutation_heteroplasmy: f64,
    /// Fraction of mtDNA carrying the age-associated 4,977-bp common deletion ($0.0$ to $1.0$).
    pub common_deletion_4977bp_fraction: f64,
    /// Inner mitochondrial membrane potential $\Delta\Psi_m$ in millivolts (healthy $\approx 175\text{ mV}$).
    pub membrane_potential_delta_psi_mv: f64,
}

/// Type alias for backward compatibility.
pub type MitochondrialPopulation = MitochondrialPopulationState;

impl MitochondrialPopulationState {
    /// Estimates a baseline somatic fibroblast mitochondrial state as a function of chronological donor age.
    pub fn from_somatic_donor_age(age_years: f64) -> Self {
        let age = age_years.clamp(0.0, 120.0);
        // Point mutations and 4,977-bp deletions accumulate non-linearly after age 35
        let point_het = (0.008 + 0.000035 * age.powf(1.85)).clamp(0.005, 0.45);
        let deletion_het = (0.002 + 0.000022 * age.powf(1.90)).clamp(0.001, 0.35);
        let total_burden = (point_het + deletion_het).clamp(0.006, 0.80);
        let delta_psi = (178.0 - total_burden * 95.0).clamp(85.0, 180.0);
        let copies = (2_400.0 - age * 8.5).clamp(900.0, 2_500.0) as u64;
        // Electron transport chain leakage increases with mutant heteroplasmy and age
        let ros = (0.04 + total_burden * 1.45 + age * 0.0015).clamp(0.02, 0.95);

        Self {
            heteroplasmy_fraction: total_burden,
            mtdna_copy_number: copies,
            ros_level: ros,
            total_copy_number: copies,
            point_mutation_heteroplasmy: point_het,
            common_deletion_4977bp_fraction: deletion_het,
            membrane_potential_delta_psi_mv: delta_psi,
        }
    }

    /// Creates a youthful somatic baseline state (equivalent to neonate / young adult).
    pub fn youthful_somatic() -> Self {
        Self::from_somatic_donor_age(20.0)
    }

    /// Creates a high-copy mature MII oocyte mitochondrial state.
    pub fn mature_oocyte_default() -> Self {
        Self {
            heteroplasmy_fraction: 0.008,
            mtdna_copy_number: 115_000,
            ros_level: 0.03,
            total_copy_number: 115_000,
            point_mutation_heteroplasmy: 0.006,
            common_deletion_4977bp_fraction: 0.002,
            membrane_potential_delta_psi_mv: 177.0,
        }
    }

    /// Total combined mutant mtDNA fraction (`heteroplasmy_fraction`).
    pub fn total_mutant_heteroplasmy(&self) -> f64 {
        self.heteroplasmy_fraction.clamp(0.0, 1.0)
    }

    /// Relative ATP oxidative phosphorylation capacity ($0.0$ to $1.0$).
    ///
    /// Reflects the classic biochemical threshold effect (Taylor & Turnbull 2005): OXPHOS capacity
    /// remains buffered until mutant heteroplasmy approaches $\sim 50\text{--}60\%$, above which ATP collapses.
    pub fn oxphos_atp_capacity(&self) -> f64 {
        let h = self.total_mutant_heteroplasmy();
        // Sigmoidal biochemical threshold around h = 0.60
        let threshold_factor = 1.0 / (1.0 + ((h - 0.58) / 0.10).exp());
        let psi_factor = ((self.membrane_potential_delta_psi_mv - 80.0) / 95.0).clamp(0.1, 1.0);
        (0.75 * threshold_factor + 0.25 * psi_factor).clamp(0.05, 1.0)
    }
}

/// Simulates somatic `PINK1`/`Parkin`-mediated mitophagy coupled with `PGC-1alpha`/`TFAM` biogenesis.
///
/// Depolarized mitochondria carrying 4,977-bp deletions or pathogenic point mutations stabilize
/// `PINK1` on the outer mitochondrial membrane, recruiting `Parkin` E3 ubiquitin ligase for selective
/// autophagic clearance, while `TFAM` biogenesis replenishes wild-type `mtDNA` copy number.
pub fn simulate_pink1_parkin_mitophagy(
    state: &MitochondrialPopulationState,
    pink1_parkin_activity: f64,
    tfam_biogenesis_boost: f64,
) -> MitochondrialPopulationState {
    let mitophagy = pink1_parkin_activity.clamp(0.0, 1.0);
    let biogenesis = tfam_biogenesis_boost.clamp(0.0, 2.0);

    // 4,977-bp deletions cause steeper depolarization and are cleared more aggressively than point mutations
    let deletion_retention = (1.0 - 0.82 * mitophagy).clamp(0.08, 1.0);
    let point_retention = (1.0 - 0.72 * mitophagy).clamp(0.12, 1.0);

    let post_point = (state.point_mutation_heteroplasmy * point_retention).clamp(0.002, 1.0);
    let post_deletion =
        (state.common_deletion_4977bp_fraction * deletion_retention).clamp(0.0005, 1.0);
    let post_het = (post_point + post_deletion).clamp(0.0025, 1.0);

    // Biogenesis restores wild-type mtDNA copy number toward youthful somatic homeostasis (~2,400 copies)
    let cleared_fraction = (state.heteroplasmy_fraction - post_het).max(0.0);
    let surviving_copies = (state.mtdna_copy_number as f64 * (1.0 - cleared_fraction * 0.5)).max(500.0);
    let new_copies = (surviving_copies * (1.0 + 0.45 * biogenesis)).clamp(800.0, 3_200.0) as u64;

    let post_delta_psi = (178.0 - post_het * 95.0).clamp(85.0, 180.0);
    let post_ros = (state.ros_level * (1.0 - 0.68 * mitophagy) + post_het * 0.15).clamp(0.02, 0.95);

    MitochondrialPopulationState {
        heteroplasmy_fraction: post_het,
        mtdna_copy_number: new_copies,
        ros_level: post_ros,
        total_copy_number: new_copies,
        point_mutation_heteroplasmy: post_point,
        common_deletion_4977bp_fraction: post_deletion,
        membrane_potential_delta_psi_mv: post_delta_psi,
    }
}

/// Outcome of germline mitochondrial bottleneck constriction and Balbiani body mitophagy selection.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct MitochondrialBottleneckReport {
    /// Effective mtDNA bottleneck size $N_b$ in PGCLCs (typically $150\text{--}250$ copies).
    pub bottleneck_size_copies: u64,
    /// Wright-Fisher binomial drift standard deviation $\sqrt{p(1-p)/N_b}$ across sibling oocytes.
    pub heteroplasmy_drift_std: f64,
    /// Resulting mitochondrial population in the mature MII oocyte after mitophagy and amplification.
    pub post_selection_population: MitochondrialPopulationState,
    /// Whether the mature oocyte satisfies clinical safety criteria (total heteroplasmy $< 0.10$
    /// and OXPHOS capacity $\ge 0.85$).
    pub passes_clinical_oxphos_threshold: bool,
}

/// Simulates the germline mitochondrial bottleneck ($N_b \approx 200$ copies) coupled with
/// `PINK1`/`Parkin`-mediated mitophagy purifying selection during PGCLC $\rightarrow$ MII oocyte maturation.
pub fn simulate_germline_mitochondrial_bottleneck(
    donor_mito: &MitochondrialPopulationState,
    pink1_parkin_mitophagy_activity: f64,
    mature_target_copies: u64,
) -> MitochondrialBottleneckReport {
    let mitophagy = pink1_parkin_mitophagy_activity.clamp(0.0, 1.0);
    let bottleneck_size = 200u64;

    // Large 4,977-bp deletions cause steeper depolarization and are cleared more strongly than point mutations
    // Reference: Lieber et al. (2019) Nature 570:380.
    let deletion_retention = (1.0 - 0.85 * mitophagy).clamp(0.08, 1.0);
    let point_retention = (1.0 - 0.75 * mitophagy).clamp(0.12, 1.0);

    let post_point = (donor_mito.point_mutation_heteroplasmy * point_retention).clamp(0.001, 1.0);
    let post_deletion =
        (donor_mito.common_deletion_4977bp_fraction * deletion_retention).clamp(0.0002, 1.0);
    let post_total = (post_point + post_deletion).min(1.0);

    let drift_var = (post_total * (1.0 - post_total)) / (bottleneck_size as f64);
    let heteroplasmy_drift_std = drift_var.sqrt();

    let post_delta_psi = (178.0 - post_total * 95.0).clamp(85.0, 180.0);
    let post_ros = (0.03 + post_total * 1.2).clamp(0.02, 0.95);
    let post_pop = MitochondrialPopulationState {
        heteroplasmy_fraction: post_total,
        mtdna_copy_number: mature_target_copies,
        ros_level: post_ros,
        total_copy_number: mature_target_copies,
        point_mutation_heteroplasmy: post_point,
        common_deletion_4977bp_fraction: post_deletion,
        membrane_potential_delta_psi_mv: post_delta_psi,
    };

    let passes_clinical_oxphos_threshold =
        post_pop.total_mutant_heteroplasmy() < 0.10 && post_pop.oxphos_atp_capacity() >= 0.85;

    MitochondrialBottleneckReport {
        bottleneck_size_copies: bottleneck_size,
        heteroplasmy_drift_std,
        post_selection_population: post_pop,
        passes_clinical_oxphos_threshold,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_somatic_age_accumulates_mtdna_heteroplasmy_and_deletions() {
        let young = MitochondrialPopulationState::from_somatic_donor_age(22.0);
        let aged = MitochondrialPopulationState::from_somatic_donor_age(75.0);

        assert!(young.total_mutant_heteroplasmy() < 0.03);
        assert!(aged.total_mutant_heteroplasmy() > 0.15);
        assert!(aged.ros_level > young.ros_level);
        assert!(aged.mtdna_copy_number < young.mtdna_copy_number);
        assert!(aged.common_deletion_4977bp_fraction > young.common_deletion_4977bp_fraction * 5.0);
        assert!(young.oxphos_atp_capacity() > aged.oxphos_atp_capacity());
    }

    #[test]
    fn test_pink1_parkin_mitophagy_purifies_aged_donor_mitochondria() {
        let aged_donor = MitochondrialPopulationState::from_somatic_donor_age(72.0);
        assert!(aged_donor.total_mutant_heteroplasmy() > 0.14);

        // Somatic PINK1/Parkin mitophagy + TFAM biogenesis
        let rejuvenated_somatic = simulate_pink1_parkin_mitophagy(&aged_donor, 0.88, 1.0);
        assert!(rejuvenated_somatic.heteroplasmy_fraction < aged_donor.heteroplasmy_fraction * 0.40);
        assert!(rejuvenated_somatic.ros_level < aged_donor.ros_level * 0.60);
        assert!(rejuvenated_somatic.mtdna_copy_number > aged_donor.mtdna_copy_number);

        // Without mitophagy purifying selection (0.0), aged donor heteroplasmy fails clinical threshold
        let unpurified = simulate_germline_mitochondrial_bottleneck(&aged_donor, 0.0, 115_000);
        assert!(!unpurified.passes_clinical_oxphos_threshold);

        // With active Balbiani body PINK1/Parkin mitophagy (0.90), mutant mtDNA is purged below 10%
        let purified = simulate_germline_mitochondrial_bottleneck(&aged_donor, 0.90, 115_000);
        assert!(purified.passes_clinical_oxphos_threshold);
        assert!(
            purified.post_selection_population.total_mutant_heteroplasmy()
                < aged_donor.total_mutant_heteroplasmy() * 0.35
        );
    }
}
