# scies-bio-th 🧬

> A high-performance, dependency-free bioinformatics & computational biology toolkit for Rust.
> ส่วนหนึ่งของชุดเครื่องมือวิทยาศาสตร์การคำนวณ **scies-series**

---

## 🏛️ Features & Architecture

`scies-bio-th` provides high-performance, production-grade biological algorithms written in **100% Pure Rust with Zero External Dependencies**:

```text
scies_bio_th::
├── seq            — DNA, RNA, and Protein validated sequence primitives & K-mers
├── align          — Local (Smith-Waterman), Global (Needleman-Wunsch), Banded NW/SW, Levenshtein & Hamming
├── codon          — Universal genetic code translation, ORF & genomic region finder
├── edit           — CRISPR-Cas9/Cas12a PAM discovery, Doench 2016 CFD matrix, Prime Editing pegRNA, and ABE8e Base Editing
├── vaccine        — Personalized cancer neoantigen identification, HLA/MHC Class I binding & confidence tiers
├── motif          — Position Weight Matrix (PWM / PSSM) scanning and consensus pattern matching
├── hmm            — Hidden Markov Models (Viterbi decoding, sequence emission probabilities)
├── genome         — Genomic intervals and Interval Tree coordinate lookup
├── variant        — VCF parsing and genomic variation models (SNP, INDEL)
├── annotation     — GFF3 genomic feature parsing
├── io             — Zero-copy FASTA, FASTQ, GFF3, VCF, and PDB/mmCIF streaming readers
├── phylo          — Phylogenetic distance matrices, UPGMA tree clustering, and Newick tree export
├── popgen         — Allele & genotype frequencies, Heterozygosity, and Chi-square Hardy-Weinberg testing
├── genomics       — K-mer spectrum counting and De Bruijn Graph unitig sequence assembly
├── transcriptomics— Expression normalization: Transcripts Per Million (TPM), RPKM, and CPM
├── proteomics     — In silico proteolytic cleavage (Trypsin, Chymotrypsin, LysC), monoisotopic mass & pI
├── metagenomics   — Alpha diversity metrics: Shannon-Wiener index, Simpson, Evenness, and Chao1
│
├── anatomy        — Multi-organ physiology circuits (Cardiovascular, Cardiopulmonary, Enterohepatic, RAAS, Neuroendocrine, Immune, Energy, ECM)
├── structure      — [feature = "structure"] PDB/mmCIF parser, DSSP secondary structure, Kabsch RMSD/TM-score, Contact maps
├── network        — [feature = "network"] Biological graphs, PPI networks & hubs, KEGG pathway DAG, Brandes centrality
├── stats          — [feature = "stats"] Mann-Whitney U, Wilcoxon, Kruskal-Wallis, Fisher's Exact, Chi-square, Benjamini-Hochberg FDR
├── docking        — [feature = "docking"] 3D grid cavity detection, Lennard-Jones 6-12 Van der Waals & Coulombic electrostatics
├── epigenome      — [feature = "epigenome"] Gardiner-Garden CpG islands, Beta/M-values, ATAC-seq accessibility, ChromHMM states
├── regulatory     — [feature = "regulatory"] Gene Regulatory Networks (GRN), Enhancers, Promoters, JASPAR TFBS motifs
├── multiomics     — [feature = "multiomics"] NMF cell-type deconvolution, multi-omics biomarker disease signatures
├── singlecell     — [feature = "singlecell"] Diffusion pseudotime cell ordering, RNA velocity induction & repression
├── pharma         — [feature = "pharma"] PK/PD pharmacokinetic simulation, Hill equation dynamics, therapeutic windows
├── safety         — [feature = "safety"] p53 tumor suppressor DNA damage checkpoints & oncogenic lockout
├── aging          — [feature = "aging"] Longevity engineering, Telomeres, Telomerase, Horvath Clock, OSK reprogramming, Immune rejuvenation
├── imaging        — [feature = "imaging"] Cell segmentation and bioimaging metadata primitives
├── paleo          — [feature = "paleo"] Paleo-genomic & paleoproteomic reconstruction, ASR Felsenstein pruning, Briggs aDNA, pseudogene resurrection, chassis editing
└── ivg            — [feature = "experimental"] In Vitro Gametogenesis, SCNT cloning, Ovarian organoid folliculogenesis
```

---

## 📦 Package Information

- **Crate Name**: `scies-bio-th`
- **Rust Import Path**: `scies_bio_th`
- **Default Build**: Zero external dependencies (Pure Rust)
- **Repository**: [https://github.com/mytheeramet3dev-x/scies-bio-th](https://github.com/mytheeramet3dev-x/scies-bio-th)

---

## 🚀 Quick Start Examples

### 1. CRISPR-Cas9 Guide RNA Design & Off-Target Scoring
```rust
use scies_bio_th::edit::crispr::generate_sgrna_candidates;
use scies_bio_th::edit::offtarget::calculate_cfd_score;
use scies_bio_th::seq::Dna;

let target_dna = Dna::new(b"ATGCGATCGATCGATCGATCGGGCATCGATCGATC")?;
let candidates = generate_sgrna_candidates(&target_dna);

for sgrna in candidates {
    let cfd = calculate_cfd_score(&sgrna.sequence, &sgrna.sequence);
    println!("Candidate PAM @ {} | CFD Score: {:.2}", sgrna.pam_index, cfd.score);
}
```

### 2. In Vitro Gametogenesis (IVG) & Synthetic Oocyte Modeling
```rust
use scies_bio_th::ivg::gamete::{PgclcCulture, GameteType, InductionProtocol};
use scies_bio_th::seq::Dna;

let stem_cell_dna = Dna::new(b"ATGCGATCGATCGATC")?;
let pgclc = PgclcCulture::induce(&stem_cell_dna, 9606, InductionProtocol::Hayashi2016, 5)?;
let synthetic_oocyte = pgclc.mature_gamete(GameteType::Oocyte, 42)?;

println!("Polar Body Extruded: {}", synthetic_oocyte.quality.polar_body_extruded);
println!("Quality Score: {:.1} / 100.0", synthetic_oocyte.quality.overall_score);
```

### 3. Personalized mRNA Cancer Vaccine Design
```rust
use scies_bio_th::seq::Protein;
use scies_bio_th::vaccine::mhc::{predict_mhc_binding, HlaAllele};
use scies_bio_th::vaccine::mrna_design::MrnaVaccineConstruct;
use scies_bio_th::vaccine::neoantigen::Neoantigen;

// Somatic KRAS G12D mutation
let kras_wt = Protein::new(b"MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSY")?;
let neo = Neoantigen::from_mutation("KRAS", &kras_wt, 12, 'D', 6)?;

let vaccine = MrnaVaccineConstruct::build("PT-2026-001", &[neo])?;
println!("Codon-optimized GC%: {:.2}%", vaccine.orf_gc_content);
println!("mRNA Length: {} nt", vaccine.total_length_nt);
```

### 4. Phylogenetics & Newick Tree Export
```rust
use scies_bio_th::phylo::{DistanceMatrix, upgma};

let names = vec!["Human".into(), "Chimp".into(), "Gorilla".into()];
let mat = vec![
    0.0, 1.6, 2.4,
    1.6, 0.0, 2.4,
    2.4, 2.4, 0.0,
];
let dm = DistanceMatrix::new(names, mat)?;
let tree = upgma(&dm)?;
println!("Newick: {}", tree.to_newick());
```

### 5. Paleo-Genomic & Paleoproteomic Reconstruction (`paleo`)
```rust
use scies_bio_th::paleo::*;
use scies_bio_th::phylo::PhyloNode;
use scies_bio_th::seq::Dna;

// Archosaurian bracket for Tyrannosaurus rex COL1A1
let birds = PhyloNode::internal(vec![PhyloNode::leaf("Chicken", 68.0)], 160.0);
let tree = PhyloNode::internal(vec![birds, PhyloNode::leaf("Alligator", 228.0)], 0.0);
let seqs = vec![
    ExtantSequence::new("Chicken", Dna::new(b"ATGGGTCCCGCT")?),
    ExtantSequence::new("Alligator", Dna::new(b"ATGGGTCCCGCA")?),
];

let model = SubstitutionModel::default_vertebrate_coding(); // HKY85 + discrete Gamma
let mut locus = FelsensteinAsr::reconstruct_locus(&tree, &seqs, &model, "T_rex", "COL1A1", 66.0)?;

// Anchor with fossil collagen peptide fragment
let peptide = PeptideFragment::new(0, b"MGP", 0.98);
PaleoproteomeTranslator::anchor_with_peptides(&mut locus, &[peptide], None)?;

println!("Consensus Sequence: {}", locus.to_consensus_dna().as_str());
println!("Anchored Accuracy: {:.1}% | Pure ASR Accuracy: {:.1}%",
    locus.biomolecule_mean_confidence.unwrap_or(0.0) * 100.0,
    locus.asr_mean_confidence.unwrap_or(0.0) * 100.0);
```

### 5. In Silico Trypsin Digestion & Peptide Mass
```rust
use scies_bio_th::proteomics::{digest, Enzyme};
use scies_bio_th::seq::Protein;

let protein = Protein::new(b"MKTLLLLAVGKSAR")?;
let peptides = digest(&protein, Enzyme::Trypsin, 0);

for pep in peptides {
    println!("Peptide: {} | Mass: {:.2} Da", 
        String::from_utf8_lossy(&pep.sequence.sequence), pep.mass_da);
}
```

### 6. Real Guide Design, Multi-Organ Coupling & Phase A Hardening (v0.8.0)
```rust
use scies_bio_th::aging::{
    design_high_fidelity_longevity_package, load_real_genomic_loci_from_file,
    validate_longevity_package, HighFidelityAgeParams, LongevityStrategy,
};
use scies_bio_th::anatomy::{
    audit_clinical_reference_anchors, propagate_cross_organ_coupling,
    AnatomicalSpatialAtlas, ClinicalBiomarkerValues, WholeBodySystemicCircuits,
};

// 1. Load real genomic target loci with cryptographic SHA-256 hash
let loci = load_real_genomic_loci_from_file("testdata/t2t_chr5_telomere_tert.fasta")?;
let tert_locus = loci.iter().find(|l| l.id.contains("tert")).unwrap();

// 2. High-Fidelity multi-factor longevity intervention design
let subject = HighFidelityAgeParams {
    chronological_age: 56.0,
    cpg_methylation_fraction: 0.54,
    telomere_length_bp: 5_900.0,
    gamma_h2ax_foci: 2.2,
    sasp_burden: 0.45,
};
let report = design_high_fidelity_longevity_package(
    &subject,
    LongevityStrategy::Balanced,
    Some(&[("TERT", &tert_locus.sequence)]),
);

// 3. Pre-flight clinical constraints & safety coupling verification
let validation = validate_longevity_package(&report);
assert!(validation.is_valid);
let safety = report.safety_coupling();
println!("Stem pool: {:.1}% | Guard: {:?}", safety.stem_pool_fraction * 100.0, safety.recommended_guard_action);

// 4. Synchronized cross-organ coupling propagation
let mut atlas = AnatomicalSpatialAtlas::from_somatic_age(56.0);
let mut circuits = WholeBodySystemicCircuits::new();
let coupling = propagate_cross_organ_coupling(&mut atlas, &mut circuits);
println!("Systemic SASP Burden: {:.2}", coupling.systemic_sasp_burden);
```

### 5. Paleo-Genomic Reconstruction & De-Extinction Chassis Editing
```rust
use scies_bio_th::paleo::*;
use scies_bio_th::phylo::PhyloNode;
use scies_bio_th::seq::Dna;

// 1. Ancestral Sequence Reconstruction (ASR) via Felsenstein pruning
let birds = PhyloNode::internal(
    vec![PhyloNode::leaf("Gallus_gallus", 68.0), PhyloNode::leaf("Falco_peregrinus", 68.0)],
    160.0,
);
let tree = PhyloNode::internal(vec![birds, PhyloNode::leaf("Alligator_mississippiensis", 228.0)], 0.0);
let seqs = vec![
    ExtantSequence::new("Gallus_gallus", Dna::new(b"GGTCCCCCT")?),
    ExtantSequence::new("Falco_peregrinus", Dna::new(b"GGACCCCCT")?),
    ExtantSequence::new("Alligator_mississippiensis", Dna::new(b"GGGCCTCCA")?),
];

let mut locus = FelsensteinAsr::reconstruct_locus(
    &tree, &seqs, &SubstitutionModel::default_vertebrate_coding(), "Tyrannosaurus rex", "COL1A1", 66.0,
)?;

// 2. Anchor with fossil collagen peptide fragment (Gly-Hyp-Ala)
let peptide = PeptideFragment::new(0, b"GOA", 0.98);
PaleoproteomeTranslator::anchor_with_peptides(&mut locus, &[peptide], None)?;

// 3. Synthesize chassis edit plan against modern chicken
let mut engine = PaleoReconstructionEngine::new("Tyrannosaurus rex", "Gallus gallus");
engine.add_reconstructed_locus(locus, Dna::new(b"GGTCCCCCT")?)?;
let draft = engine.build_draft()?;
println!("Target Chassis Edits: {}", draft.chassis_edit_plan.len());
```

---

## 🧪 Testing & Verification

Run the full test suite (194+ tests passing across all features):
```bash
cargo test --all-features
```

Run paleo-genomic reconstruction test suite:
```bash
cargo test --features paleo
```

Run paleo de-extinction demonstration:
```bash
cargo run --release --example paleo_reconstruction_demo --features paleo
```

Run longevity & aging test suite:
```bash
cargo test --features aging
```

Run Criterion microbenchmarks:
```bash
cargo bench
```

---

## 📄 License

Licensed under the [MIT License](LICENSE).
