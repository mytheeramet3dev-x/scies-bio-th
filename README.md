# scies-bio-th 🧬

> A high-performance, dependency-free bioinformatics & computational biology toolkit for Rust.
> ส่วนหนึ่งของชุดเครื่องมือวิทยาศาสตร์การคำนวณ **scies-series**

---

## 🏛️ Features & Architecture

`scies-bio-th` provides high-performance, production-grade biological algorithms written in **100% Pure Rust with Zero External Dependencies**:

```text
scies_bio_th::
├── seq            — DNA, RNA, and Protein validated sequence primitives
├── align          — Local (Smith-Waterman), Global (Needleman-Wunsch), Levenshtein & Hamming distances
├── codon          — Universal genetic code translation and open reading frame (ORF) finder
├── edit           — CRISPR-Cas9/Cas12a PAM discovery, sgRNA design, CFD off-target scoring, and PCR primer design
├── ivg            — In Vitro Gametogenesis (IVG), Epigenetic reprogramming, SCNT cloning, Zona Pellucida & Artificial Womb
├── vaccine        — Personalized cancer neoantigen identification, HLA/MHC Class I binding prediction, and mRNA vaccine design
├── motif          — Position Weight Matrix (PWM / PSSM) scanning and consensus pattern matching
├── hmm            — Hidden Markov Models (Viterbi decoding, sequence emission probabilities)
├── genome         — Genomic intervals and Interval Tree coordinate lookup
├── variant        — VCF parsing and genomic variation models (SNP, INDEL)
├── annotation     — GFF3 genomic feature parsing
└── io             — Zero-copy FASTA, FASTQ, GFF3, and VCF streaming readers
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

---

## 🧪 Testing & Verification

Run the full test suite (48 tests passing):
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
