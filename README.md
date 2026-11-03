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
├── edit           — CRISPR-Cas9/Cas12a PAM discovery, Doench 2016 CFD matrix, and PCR primer design
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
├── aging          — Longevity engineering, Telomere dynamics (TTAGGG), Telomerase (TERT), and Multi-generation simulation
├── metagenomics   — Alpha diversity metrics: Shannon-Wiener index, Simpson, Evenness, and Chao1
│
├── structure      — [feature = "structure"] PDB/mmCIF parser, DSSP secondary structure, Kabsch RMSD/TM-score, Contact maps
├── network        — [feature = "network"] Biological graphs, PPI networks & hubs, KEGG pathway DAG, Brandes centrality
├── stats          — [feature = "stats"] Mann-Whitney U, Fisher's Exact (2x2), Chi-square, Benjamini-Hochberg FDR, Vector metrics
├── docking        — [feature = "docking"] 3D grid cavity detection, Lennard-Jones 6-12 Van der Waals & Coulombic electrostatics
├── epigenome      — [feature = "epigenome"] Gardiner-Garden CpG islands, Beta/M-values, ATAC-seq accessibility, ChromHMM states
└── ivg            — [feature = "experimental"] In Vitro Gametogenesis (IVG), Epigenetic reprogramming, SCNT cloning
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

### 6. Real Guide Design, Simulation & Hard Constraints (v0.6.0)
```rust
use scies_bio_th::aging::{
    design_high_fidelity_longevity_package, load_loci_from_file,
    validate_longevity_package, HighFidelityAgeParams, LongevityStrategy,
};

// 1. Load real genomic target loci from T2T-CHM13
let tert_loci = load_loci_from_file("testdata/t2t_chr5_telomere_tert.fasta")?;
let tert_dna = &tert_loci[1].1;

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
    Some(&[("TERT", tert_dna)]),
);

// 3. Pre-flight clinical constraints validation
let validation = validate_longevity_package(&report);
assert!(validation.is_valid);
println!("{}", report.export_detailed_summary());
```

---

## 🧪 Testing & Verification

Run the full test suite (75+ tests passing):
```bash
cargo test
```

Run Criterion microbenchmarks:
```bash
cargo bench
```

---

## 📄 License

Licensed under the [MIT License](LICENSE).
