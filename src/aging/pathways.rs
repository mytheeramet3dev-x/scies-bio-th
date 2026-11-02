//! Canonical longevity pathways, aging hallmarks, and curated target gene registries.
//!
//! Provides structured metadata for critical anti-aging and genomic maintenance genes
//! across DNA repair, telomere reconstitution, nutrient sensing, and senescent cell clearance.

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Biological hallmark pathway governing mammalian longevity and healthspan.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum LongevityPathway {
    /// Telomere maintenance, end-capping, and reverse transcription (TERT, TERC, POT1).
    TelomereMaintenance,
    /// Double-strand break repair, chromatin deacetylation, and genomic stability (SIRT6, SIRT1, PARP1).
    DnaRepairGenomicStability,
    /// Insulin/IGF-1 signaling, AMPK activation, and nutrient sensing (FOXO3, AMPK, MTOR).
    NutrientSensingMetabolism,
    /// Senescence-associated secretory phenotype (SASP) and cell-cycle arrest (CDKN2A/p16, TP53).
    SenescenceControl,
    /// Nuclear envelope integrity, chromatin anchoring, and laminopathies (LMNA).
    NuclearArchitecture,
    /// Mitochondrial biogenesis, mitophagy, and reactive oxygen species buffering (PGC1A, TFAM).
    MitochondrialHomeostasis,
}

/// Metadata and therapeutic genetic modification profile for a longevity target gene.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct LongevityGeneInfo {
    /// Official HUGO Gene Nomenclature Committee (HGNC) symbol.
    pub symbol: &'static str,
    /// Full descriptive biological gene name.
    pub name: &'static str,
    /// Associated aging pathway.
    pub pathway: LongevityPathway,
    /// Primary chromosome locus (GRCh38/T2T coordinate reference).
    pub chromosome: &'static str,
    /// NCBI Entrez Gene ID.
    pub entrez_id: u32,
    /// Recommended synthetic genetic intervention.
    pub recommended_action: &'static str,
    /// Biological rationale for longevity enhancement.
    pub rationale: &'static str,
}

/// Curated registry of high-impact mammalian longevity genes.
pub const CANONICAL_LONGEVITY_GENES: &[LongevityGeneInfo] = &[
    LongevityGeneInfo {
        symbol: "TERT",
        name: "Telomerase Reverse Transcriptase",
        pathway: LongevityPathway::TelomereMaintenance,
        chromosome: "chr5",
        entrez_id: 7015,
        recommended_action: "Controlled Inducible Reactivation / Promoter Optimization",
        rationale: "Synthesizes TTAGGG telomeric repeats to overcome somatic replicative senescence.",
    },
    LongevityGeneInfo {
        symbol: "TERC",
        name: "Telomerase RNA Component",
        pathway: LongevityPathway::TelomereMaintenance,
        chromosome: "chr3",
        entrez_id: 7012,
        recommended_action: "Sustained Baseline Expression",
        rationale: "Serves as the vital RNA template for telomere repeat elongation by TERT.",
    },
    LongevityGeneInfo {
        symbol: "POT1",
        name: "Protection of Telomeres 1",
        pathway: LongevityPathway::TelomereMaintenance,
        chromosome: "chr7",
        entrez_id: 25913,
        recommended_action: "Structural Stabilization",
        rationale: "Binds single-stranded telomeric 3' overhangs, preventing ATR kinase activation.",
    },
    LongevityGeneInfo {
        symbol: "SIRT6",
        name: "Sirtuin 6 (H3K9/H3K56 Deacetylase & Mono-ADP-ribosyltransferase)",
        pathway: LongevityPathway::DnaRepairGenomicStability,
        chromosome: "chr19",
        entrez_id: 51548,
        recommended_action: "Moderate Transcriptional Upregulation / Gain-of-Function",
        rationale: "Enhances non-homologous end joining and homologous recombination DSB repair; extends lifespan in mammals.",
    },
    LongevityGeneInfo {
        symbol: "SIRT1",
        name: "Sirtuin 1 (NAD-dependent Protein Deacetylase)",
        pathway: LongevityPathway::DnaRepairGenomicStability,
        chromosome: "chr10",
        entrez_id: 23411,
        recommended_action: "Enhance NAD+ Substrate Sensitivity",
        rationale: "Promotes mitochondrial biogenesis and deacetylates PGC-1alpha and FOXO factors.",
    },
    LongevityGeneInfo {
        symbol: "FOXO3",
        name: "Forkhead Box O3",
        pathway: LongevityPathway::NutrientSensingMetabolism,
        chromosome: "chr6",
        entrez_id: 2309,
        recommended_action: "Centenarian G-Allele (rs2802292) Knock-In",
        rationale: "Master transcription factor driving antioxidant defenses, autophagy, and proteostasis.",
    },
    LongevityGeneInfo {
        symbol: "CDKN2A",
        name: "Cyclin Dependent Kinase Inhibitor 2A (p16INK4a / p14ARF)",
        pathway: LongevityPathway::SenescenceControl,
        chromosome: "chr9",
        entrez_id: 1029,
        recommended_action: "Senolytic Inducible Suicide Circuit (INK-ATTAC style)",
        rationale: "Enables targeted ablation of high-p16 senescent cells secreting toxic inflammatory SASP factors.",
    },
    LongevityGeneInfo {
        symbol: "TP53",
        name: "Tumor Protein P53",
        pathway: LongevityPathway::SenescenceControl,
        chromosome: "chr17",
        entrez_id: 7157,
        recommended_action: "Pachyderm Extra Copy / Super-p53 Dosage Elevation",
        rationale: "Provides robust cancer surveillance and apoptosis clearance for cells acquiring oncogenic lesions.",
    },
    LongevityGeneInfo {
        symbol: "ATM",
        name: "ATM Serine/Threonine Kinase",
        pathway: LongevityPathway::DnaRepairGenomicStability,
        chromosome: "chr11",
        entrez_id: 472,
        recommended_action: "Kinase Phosphorylation Optimization",
        rationale: "Master coordinator of DNA double-strand break response and chromatin relaxation.",
    },
    LongevityGeneInfo {
        symbol: "ATR",
        name: "ATR Serine/Threonine Kinase",
        pathway: LongevityPathway::DnaRepairGenomicStability,
        chromosome: "chr3",
        entrez_id: 545,
        recommended_action: "Replication Stress Signaling Regulation",
        rationale: "Responds to single-stranded DNA and stalled replication forks, preventing telomere fragility.",
    },
    LongevityGeneInfo {
        symbol: "PARP1",
        name: "Poly(ADP-Ribose) Polymerase 1",
        pathway: LongevityPathway::DnaRepairGenomicStability,
        chromosome: "chr1",
        entrez_id: 142,
        recommended_action: "NAD+ Sparing & Single-Strand Break Repair Tuning",
        rationale: "Detects single-strand DNA breaks and initiates base excision repair while conserving cellular NAD+.",
    },
    LongevityGeneInfo {
        symbol: "LMNA",
        name: "Lamin A/C (Nuclear Lamina Structural Component)",
        pathway: LongevityPathway::NuclearArchitecture,
        chromosome: "chr1",
        entrez_id: 4000,
        recommended_action: "Cryptic Splice Site Correction (Progerin Resistance)",
        rationale: "Prevents aberrant Hutchinson-Gilford progerin accumulation and maintains nuclear membrane rigidity.",
    },
    LongevityGeneInfo {
        symbol: "MTOR",
        name: "Mechanistic Target Of Rapamycin Kinase",
        pathway: LongevityPathway::NutrientSensingMetabolism,
        chromosome: "chr1",
        entrez_id: 2475,
        recommended_action: "Downregulated Basal mTORC1 Activity (Rapamycin Mimetic)",
        rationale: "Inhibits cap-dependent translation while activating macroautophagy and lysosomal biogenesis.",
    },
    LongevityGeneInfo {
        symbol: "AMPK",
        name: "5'-AMP-Activated Protein Kinase (PRKAA1/Alpha 1)",
        pathway: LongevityPathway::NutrientSensingMetabolism,
        chromosome: "chr5",
        entrez_id: 5562,
        recommended_action: "Constitutive Low-Level Energy Sensor Activation",
        rationale: "Phosphorylates ULK1 and FOXO3, promoting mitochondrial quality control and lipid oxidation.",
    },
];

/// Looks up curated longevity gene information by official gene symbol.
pub fn find_longevity_gene(symbol: &str) -> Option<&'static LongevityGeneInfo> {
    CANONICAL_LONGEVITY_GENES.iter().find(|g| g.symbol.eq_ignore_ascii_case(symbol))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_canonical_genes() {
        let tert = find_longevity_gene("TERT").unwrap();
        assert_eq!(tert.entrez_id, 7015);
        assert_eq!(tert.pathway, LongevityPathway::TelomereMaintenance);

        let sirt6 = find_longevity_gene("sirt6").unwrap(); // Case-insensitive
        assert_eq!(sirt6.entrez_id, 51548);
        assert_eq!(sirt6.pathway, LongevityPathway::DnaRepairGenomicStability);

        let lmna = find_longevity_gene("LMNA").unwrap();
        assert_eq!(lmna.entrez_id, 4000);
        assert_eq!(lmna.pathway, LongevityPathway::NuclearArchitecture);

        let mtor = find_longevity_gene("mtor").unwrap();
        assert_eq!(mtor.entrez_id, 2475);
        assert_eq!(mtor.pathway, LongevityPathway::NutrientSensingMetabolism);

        assert!(find_longevity_gene("UNKNOWN_GENE").is_none());
    }
}
