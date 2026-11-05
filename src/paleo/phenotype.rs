//! Layer 3: Commutative, Group-Discounted Bayesian Phenotypic & Environmental Selection.
//!
//! Translates non-molecular fossil evidence (osteology, bone histology, body-mass scaling,
//! and stable isotopes $\delta^{13}\text{C}, \delta^{15}\text{N}, \delta^{18}\text{O}$) into
//! pathway-level natural selection posteriors ($dN/dS$ priority weights).
//!
//! ### Mathematical Properties
//! - **Commutative**: Evidence order never affects the posterior distribution.
//! - **Soft Bayesian Pooling**: No hard-cut fragility; conflicting evidence smoothly shifts log-odds.
//! - **Correlation Damping**: Evidences within the same [`CorrelationGroup`] $g$ are scaled by
//!   $\alpha_g = 1 / \sqrt{|g|}$ to prevent overconfidence from non-independent fossil measurements.

use std::collections::HashMap;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Correlation group for discounting non-independent paleontological measurements.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum CorrelationGroup {
    /// Bone histology, growth rings (LAGs), pneumaticity, and dentition from skeletal remains.
    OsteologyHistology,
    /// Stable isotope ratios ($\delta^{13}\text{C}, \delta^{15}\text{N}, \delta^{18}\text{O}$) from enamel/collagen.
    GeochemicalIsotope,
    /// Allometric body mass and biomechanical muscle attachment scaling.
    MorphometricScaling,
    /// Depositional environment and paleoclimate proxies.
    PaleoEnvironment,
}

/// Biologically grounded functional gene pathways modulated by paleontological evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum FunctionalPathway {
    /// Peto's Paradox megafauna tumor suppression & DNA repair (`TP53`, `LIF`, `ATM`, `SIRT6`).
    OncogenicSuppressionDnaRepair,
    /// Endothermic/mesothermic oxidative metabolism & thermogenesis (`UCP1`, `PPARGC1A`, `COX1`).
    EndothermicMetabolism,
    /// Cryogenic thermal sensing, subcutaneous lipid storage & hair keratin (`TRPM8`, `LEP`, `KRT71`).
    ColdAdaptationLipidStorage,
    /// Hypercarnivorous proteolysis & umami taste receptor retention (`PRSS1`, `CTRC`, `TAS1R1`).
    CarnivoryProteolysis,
    /// Herbivorous plant/xenobiotic detoxification & starch digestion (`AMY1`, `CYP2E1`, `TAS2R`).
    HerbivoryCelluloseDigestion,
    /// Biomineralized tooth enamel & dentin matrix formation (`ENAM`, `AMELX`, `MMP20`).
    MineralizedDentitionEnamel,
    /// Rapid radial osteogenesis & postcranial skeletal pneumaticity (`IGF1`, `RUNX2`, `BMP4`).
    SkeletalPneumaticityGrowth,
}

impl FunctionalPathway {
    /// Returns all canonical pathways in deterministic order.
    pub fn all() -> [Self; 7] {
        [
            Self::OncogenicSuppressionDnaRepair,
            Self::EndothermicMetabolism,
            Self::ColdAdaptationLipidStorage,
            Self::CarnivoryProteolysis,
            Self::HerbivoryCelluloseDigestion,
            Self::MineralizedDentitionEnamel,
            Self::SkeletalPneumaticityGrowth,
        ]
    }

    /// Maps a gene symbol (case-insensitive) to its primary `FunctionalPathway`.
    pub fn from_gene_symbol(symbol: &str) -> Option<Self> {
        match symbol.to_ascii_uppercase().as_str() {
            "TP53" | "LIF" | "ATM" | "SIRT6" | "BRCA1" => Some(Self::OncogenicSuppressionDnaRepair),
            "UCP1" | "PPARGC1A" | "COX1" | "COX2" => Some(Self::EndothermicMetabolism),
            "TRPM8" | "LEP" | "FADS2" | "KRT71" => Some(Self::ColdAdaptationLipidStorage),
            "PRSS1" | "CTRC" | "TAS1R1" | "CELA1" => Some(Self::CarnivoryProteolysis),
            "AMY1" | "CYP2E1" | "TAS2R" => Some(Self::HerbivoryCelluloseDigestion),
            "ENAM" | "AMELX" | "MMP20" | "DSPP" => Some(Self::MineralizedDentitionEnamel),
            "IGF1" | "RUNX2" | "BMP4" | "COL1A1" | "COL1A2" => {
                Some(Self::SkeletalPneumaticityGrowth)
            }
            _ => None,
        }
    }
}

/// Trait implemented by all non-molecular paleontological evidence sources.
pub trait PhenotypicEvidence {
    /// Descriptive identifier for audit logging.
    fn name(&self) -> &'static str;

    /// Correlation group used for $1/\sqrt{|g|}$ non-independence discounting.
    fn correlation_group(&self) -> CorrelationGroup;

    /// Specimen preservation and measurement reliability in `[0.0, 1.0]`.
    fn sample_reliability(&self) -> f64;

    /// Returns signed log-likelihood shifts $\ell_e(M)$ for relevant functional pathways.
    fn pathway_log_likelihoods(&self) -> Vec<(FunctionalPathway, f64)>;
}

/// Bone histology and structural density classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum BoneHistologyType {
    /// Rapid, vascularized woven-fibered bone characteristic of tachymetabolic endotherms (mammals, theropods, birds).
    Fibrolamellar,
    /// Slow, cyclical lamellar-zonal bone characteristic of ectothermic reptiles.
    LamellarZonal,
}

/// Skeletal osteology, histology, and body-mass scaling evidence.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct OsteologyEvidence {
    /// Specimen label (e.g. `"MOR 1125 Femur"`).
    pub specimen_label: String,
    /// Histological microstructure.
    pub histology: Option<BoneHistologyType>,
    /// Presence of postcranial skeletal pneumaticity (air-sac foramina in vertebrae).
    pub has_pneumatic_foramina: bool,
    /// Presence of mineralized enamel-coated teeth in alveoli.
    pub has_enamel_dentition: bool,
    /// Estimated adult body mass in kg (from stylopodial circumference scaling).
    pub estimated_body_mass_kg: Option<f64>,
    /// Specimen preservation quality `[0.0, 1.0]`.
    pub preservation_reliability: f64,
}

impl PhenotypicEvidence for OsteologyEvidence {
    fn name(&self) -> &'static str {
        "OsteologyHistologyEvidence"
    }

    fn correlation_group(&self) -> CorrelationGroup {
        CorrelationGroup::OsteologyHistology
    }

    fn sample_reliability(&self) -> f64 {
        self.preservation_reliability.clamp(0.05, 1.0)
    }

    fn pathway_log_likelihoods(&self) -> Vec<(FunctionalPathway, f64)> {
        let mut shifts = Vec::new();

        // 1. Histology -> Endothermic metabolism & osteogenesis rate
        if let Some(hist) = self.histology {
            match hist {
                BoneHistologyType::Fibrolamellar => {
                    shifts.push((FunctionalPathway::EndothermicMetabolism, 1.8));
                    shifts.push((FunctionalPathway::SkeletalPneumaticityGrowth, 1.4));
                }
                BoneHistologyType::LamellarZonal => {
                    shifts.push((FunctionalPathway::EndothermicMetabolism, -1.6));
                }
            }
        }

        // 2. Skeletal pneumaticity -> avian-style respiratory air sacs & high aerobic scope
        if self.has_pneumatic_foramina {
            shifts.push((FunctionalPathway::SkeletalPneumaticityGrowth, 1.6));
            shifts.push((FunctionalPathway::EndothermicMetabolism, 1.2));
        }

        // 3. Dentition -> functional enamel matrix genes (ENAM, AMELX, MMP20)
        if self.has_enamel_dentition {
            shifts.push((FunctionalPathway::MineralizedDentitionEnamel, 2.5));
        } else {
            shifts.push((FunctionalPathway::MineralizedDentitionEnamel, -2.0));
        }

        // 4. Peto's Paradox allometric body mass scaling (Sulak et al. 2016)
        if let Some(mass_kg) = self.estimated_body_mass_kg {
            if mass_kg >= 1_000.0 {
                let log10_mass = mass_kg.log10(); // 3.0 for 1 ton, 3.85 for 7 tons
                let peto_shift = ((log10_mass - 2.0) * 1.1).clamp(0.5, 3.0);
                shifts.push((FunctionalPathway::OncogenicSuppressionDnaRepair, peto_shift));
                shifts.push((FunctionalPathway::SkeletalPneumaticityGrowth, 1.2));
            }
        }

        shifts
    }
}

/// Geochemical stable isotope ratios ($\delta^{13}\text{C}, \delta^{15}\text{N}, \delta^{18}\text{O}$) from bioapatite or collagen.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct IsotopeEvidence {
    /// $\delta^{13}\text{C}$ (‰ VPDB): dietary carbon source.
    pub delta13c_permil: Option<f64>,
    /// $\delta^{15}\text{N}$ (‰ AIR): trophic level enrichment (+3 to +5‰ per trophic step).
    pub delta15n_permil: Option<f64>,
    /// $\delta^{18}\text{O}$ (‰ VSMOW): meteoric water & ambient paleotemperature proxy.
    pub delta18o_permil: Option<f64>,
    /// Geochemical sample integrity (e.g. C:N atomic ratio inside `2.9..3.6` indicates low diagenesis).
    pub reliability: f64,
}

impl PhenotypicEvidence for IsotopeEvidence {
    fn name(&self) -> &'static str {
        "GeochemicalIsotopeEvidence"
    }

    fn correlation_group(&self) -> CorrelationGroup {
        CorrelationGroup::GeochemicalIsotope
    }

    fn sample_reliability(&self) -> f64 {
        self.reliability.clamp(0.05, 1.0)
    }

    fn pathway_log_likelihoods(&self) -> Vec<(FunctionalPathway, f64)> {
        let mut shifts = Vec::new();

        // 1. Nitrogen-15 trophic enrichment -> Carnivory vs Herbivory
        if let Some(d15n) = self.delta15n_permil {
            if d15n >= 8.0 {
                shifts.push((FunctionalPathway::CarnivoryProteolysis, 2.1));
                shifts.push((FunctionalPathway::HerbivoryCelluloseDigestion, -1.8));
            } else if d15n <= 5.5 {
                shifts.push((FunctionalPathway::HerbivoryCelluloseDigestion, 1.9));
                shifts.push((FunctionalPathway::CarnivoryProteolysis, -1.7));
            }
        }

        // 2. Carbon-13 plant diet signal when paired with non-predatory N-15
        if let Some(d13c) = self.delta13c_permil {
            let is_high_trophic = self.delta15n_permil.map(|n| n >= 8.0).unwrap_or(false);
            if (-29.0..=-18.0).contains(&d13c) && !is_high_trophic {
                shifts.push((FunctionalPathway::HerbivoryCelluloseDigestion, 1.2));
            }
        }

        // 3. Oxygen-18 paleoclimate temperature signal
        if let Some(d18o) = self.delta18o_permil {
            if d18o <= -10.0 {
                // Sub-arctic / glacial steppe environment
                shifts.push((FunctionalPathway::ColdAdaptationLipidStorage, 2.4));
                shifts.push((FunctionalPathway::EndothermicMetabolism, 1.5));
            } else if d18o >= 16.0 {
                // Warm Mesozoic greenhouse regime
                shifts.push((FunctionalPathway::ColdAdaptationLipidStorage, -1.8));
            }
        }

        shifts
    }
}

/// Posterior selection probabilities and $dN/dS$ ($\omega$) modulation across functional pathways.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PathwaySelectionProfile {
    /// Posterior probability $P(M \mid \mathcal{E}) \in (0, 1)$ for each pathway.
    pub pathway_posteriors: HashMap<FunctionalPathway, f64>,
    /// Number of evidence items pooled per correlation group.
    pub group_counts: HashMap<CorrelationGroup, usize>,
}

impl PathwaySelectionProfile {
    /// Computes the commutative, group-discounted Bayesian posterior across all pathways.
    ///
    /// Adding evidences in any permutation yields the exact same `PathwaySelectionProfile`.
    pub fn from_evidences(evidences: &[&dyn PhenotypicEvidence]) -> Self {
        let mut group_counts: HashMap<CorrelationGroup, usize> = HashMap::new();
        let mut group_accum: HashMap<(CorrelationGroup, FunctionalPathway), f64> = HashMap::new();

        for ev in evidences {
            let g = ev.correlation_group();
            *group_counts.entry(g).or_insert(0) += 1;
            let r = ev.sample_reliability();
            for (pathway, ll_shift) in ev.pathway_log_likelihoods() {
                if ll_shift.is_finite() {
                    *group_accum.entry((g, pathway)).or_insert(0.0) += r * ll_shift;
                }
            }
        }

        let mut pathway_posteriors = HashMap::new();
        for pathway in FunctionalPathway::all() {
            let mut total_log_odds = 0.0; // Neutral prior P_0 = 0.50 -> logit(P_0) = 0.0
            for (&g, &count) in &group_counts {
                if count > 0 {
                    let discount = 1.0 / (count as f64).sqrt();
                    let raw_sum = group_accum.get(&(g, pathway)).copied().unwrap_or(0.0);
                    total_log_odds += discount * raw_sum;
                }
            }
            // Logistic sigmoid: P = 1 / (1 + exp(-log_odds))
            let post = 1.0 / (1.0 + (-total_log_odds).exp());
            pathway_posteriors.insert(pathway, post);
        }

        Self {
            pathway_posteriors,
            group_counts,
        }
    }

    /// Returns the posterior selection probability for `pathway` (defaults to neutral `0.50`).
    pub fn posterior_for(&self, pathway: FunctionalPathway) -> f64 {
        self.pathway_posteriors
            .get(&pathway)
            .copied()
            .unwrap_or(0.50)
    }

    /// Returns the selection priority multiplier for a given gene symbol based on its pathway posterior.
    pub fn gene_priority_weight(&self, gene_symbol: &str) -> f64 {
        self.selection_intensity_omega(gene_symbol)
    }

    /// Computes the selection intensity multiplier $\omega$ ($dN/dS$ weighting factor) for a gene.
    ///
    /// Values range from $0.5$ (relaxed/pseudogenizing) through $1.0$ (neutral baseline) up to $1.5$ (accelerated/adaptive).
    ///
    /// > **Architectural Invariant**: Phenotypic evidence strictly modulates selection intensity ($\omega$)
    /// > and locus presence/absence priors. It never directly overrides site-level nucleotide probabilities $P(A,C,G,T)$.
    pub fn selection_intensity_omega(&self, gene_symbol: &str) -> f64 {
        match FunctionalPathway::from_gene_symbol(gene_symbol) {
            Some(pathway) => 0.5 + self.posterior_for(pathway),
            None => 1.0,
        }
    }

    /// Evaluates the macro-evolutionary gene presence/absence prior $P(\text{Gene is functional} \mid \mathcal{E})$.
    ///
    /// For example:
    /// - Toothed dentition fossil evidence produces high prior ($> 0.85$) for enamel genes (`ENAM`, `AMELX`).
    /// - Arctic permafrost evidence produces high prior ($> 0.80$) for thermogenic adaptation (`UCP1`, `TRPM8`).
    pub fn gene_presence_absence_prior(&self, gene_symbol: &str) -> f64 {
        match FunctionalPathway::from_gene_symbol(gene_symbol) {
            Some(pathway) => {
                let post = self.posterior_for(pathway);
                post.clamp(0.05, 0.95)
            }
            None => 0.50,
        }
    }
}
