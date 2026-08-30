# scies-bio-th — v0.1 Design Document

> A high-performance bio toolkit for Rust.
> ส่วนหนึ่งของซีรีส์ `scies` (ต่อจาก `scies-math`)

---

## 1. เป้าหมายของโปรเจกต์

เป้าหมายระยะยาวคือทำให้ `scies-bio-th` ครอบคลุมงานด้าน biology ให้กว้างที่สุดเท่าที่เป็นไปได้ในรูปแบบ library โดยไม่พยายามยัดทุกอย่างไว้ใน release แรก

ขอบเขตใหญ่ของโปรเจกต์นี้แบ่งเป็น:
- **bioinformatics / computational biology**: sequence, alignment, annotation, omics, phylo, popgen
- **molecular biology tooling**: CRISPR helpers, primer design, cloning/plasmid utilities, assay helpers
- **structural / systems biology**: protein structure, pathways, networks, simulation
- **clinical / translational bio**: variant interpretation, cohort-level workflows, biomarker helpers
- **ecology / evolution**: taxonomy, metagenomics, ecological and evolutionary analysis

v0.1 ยังเป็น **core computation library** ไม่ใช่ฐานข้อมูล และไม่ใช่ wet-lab protocol engine
แต่ต้องวางฐานข้อมูล/ชนิดข้อมูล/อัลกอริทึมพื้นฐานให้ดีพอ เพื่อให้ v0.2+ ขยายไปยังหมวด biology อื่น ๆ ได้โดยไม่ต้องรื้อ API

### Design Goals
- **Zero dependencies** ใน default build (pure `std`)
- **Generic** over element type ที่จำเป็น (เช่น scoring matrix เป็น `i32`/`f64`)
- **Const generics** สำหรับ fixed-size k-mers / codons
- **Lazy / streaming iterators** สำหรับไฟล์ genomic ขนาดใหญ่ (ไม่โหลดทั้งไฟล์เข้า memory)
- **Benchmarked** ด้วย Criterion ทุก module ที่มี hot path
- Doc-comment (`///`) เป็นทั้งเอกสารและ doctest ในตัว แทนการแยกเขียนเอกสารต่างหาก
- ถ้าต้องใช้ math utility ให้ดึงจาก `scies-math` แทนการเขียนซ้ำใน crate นี้

---

## 2. Module Tree (v0.1 scope)

```
scies_bio_th::
├── seq          — Dna, Rna, Protein types, k-mer, reverse complement, translation, GC content
├── io           — zero-dep parsers: FASTA, FASTQ (+quality scores)
├── align        — Needleman-Wunsch, Smith-Waterman, edit distance
├── codon        — codon tables (standard genetic code), ORF finding
└── edit
    ├── crispr        — PAM site finding (SpCas9 NGG, Cas12a TTTV), sgRNA candidate generation
    ├── offtarget      — off-target scoring (mismatch-based, CFD-style)
    └── primer         — PCR primer design: Tm calculation, GC clamp, dimer/hairpin check
```

> โมดูลอื่นที่ออกแบบไว้ในภาพรวมทั้งซีรีส์ (`hmm`, `motif`, `phylo`, `popgen`, `epi`,
> `struct_bio`, `expr`, `graph_bio`, `vector`, `crispr_repair`, `variant`, `annotation`, `metagenomics`,
> `cloning`, `pathway`, `network`, `taxonomy`, `imaging`, `cell`, `microbio`, `ecology`, `clinical` ฯลฯ)
> **ไม่อยู่ใน scope ของ v0.1**
> จะทยอยเพิ่มใน v0.2+ ตามลำดับที่ระบุในหัวข้อ 5

### หมายเหตุขอบเขต `edit`
โมดูลนี้เป็นเครื่องมือ **คำนวณ/ค้นหา sequence เท่านั้น** (หา PAM site, คำนวณ Tm,
สร้าง candidate list) คล้ายเครื่องมือมาตรฐานอย่าง CRISPOR / Benchling / SnapGene
**ไม่แตะ wet-lab protocol หรือ synthesis** เพื่อให้ scope เหมาะสมกับสิ่งที่เป็น
"library คำนวณ" จริง ๆ

---

## 3. Feature Flags

| Flag | Enables |
|---|---|
| `serde` | Serialize/Deserialize สำหรับ core types (`Dna`, `Rna`, `Protein`) |
| `parallel` | multithreaded alignment (ใช้ `std::thread` เท่านั้น ไม่ดึง `rayon`) |
| `interop` | bridge เชื่อมกับ `scies-math` (optional dep) สำหรับ use cases ที่ต้องใช้ statistical / numerical helpers |

Default build: **ไม่เปิด feature ใด ๆ** (`default = []`)

---

## 4. Cargo.toml Skeleton

```toml
[package]
name = "scies-bio-th"
version = "0.1.0"
edition = "2021"
license = "MIT"
description = "A high-performance, dependency-free bioinformatics toolkit for Rust."
repository = "https://github.com/mytheeramet3dev-x/scies-bio-th"

[dependencies]
scies-math = { version = "0.2", optional = true }
serde = { version = "1", features = ["derive"], optional = true }

[features]
default = []
serde = ["dep:serde"]
parallel = []
interop = ["dep:scies-math"]

[dev-dependencies]
criterion = "0.5"

[[bench]]
name = "align_bench"
harness = false
```

---

## 5. ลำดับการ Implement

โค้ดก่อนเอกสารเต็ม แต่ร่าง public API แบบ signature-only (`todo!()`) ก่อนลง implementation จริงในทุก module

1. **`seq`** — ฐานของทุก module อื่น (`Dna`, `Rna`, `Protein`, k-mer, reverse complement, translation, GC content)
2. **`io::fasta`** (+ `io::fastq`) — ต้องมี data เข้าเพื่อทดสอบ `seq` จริงจัง
3. **`edit::crispr`** — ใช้ `seq` โดยตรง ทดสอบง่าย ไม่ต้องพึ่ง alignment
4. **`align`** (Needleman-Wunsch → Smith-Waterman → edit distance)
5. **`codon`** — codon table + ORF finding (ใช้ `seq`)
6. **`edit::offtarget`** — ต้องมี `align` มาก่อน เพื่อคำนวณ mismatch กับ genome
7. **`edit::primer`** — Tm calculation, GC clamp, hairpin/dimer check (ใช้ `seq` + `align` บางส่วน)

แต่ละ module: เขียน signature → `todo!()` → unit test → implement จริง → doc-comment + doctest คู่ขนาน → benchmark (เฉพาะ hot path)

### Long-term Biology Coverage Roadmap

#### v0.1 Core Biology Foundation
- `seq`, `io`, `align`, `codon`, `edit`
- immutable/borrowed sequence views
- metadata types สำหรับ species, chromosome, strand, locus
- error model กลางสำหรับ parser, alignment, translation, design tools

#### v0.2 Sequence Biology Expansion
- `hmm`, `motif`
- `annotation` primitives
- `variant` primitives สำหรับ SNP/indel representation
- `genome` utilities สำหรับ reference segments, intervals, masks
- `interop` กับ `scies-math` สำหรับ scoring, statistics, optimization, and numerical routines

#### v0.3 Omics and Comparative Biology
- `genomics`: assembly helpers, read mapping helpers, repeat handling
- `transcriptomics`: splicing, isoform handling, RNA-seq summaries
- `proteomics`: peptide utilities, mass-aware helpers, protein feature extraction
- `metagenomics`: taxonomic assignment helpers, abundance tables
- `phylo`: trees, distances, ancestral reconstruction primitives
- `popgen`: allele frequency, Hardy-Weinberg, drift/selection helpers

#### v0.4 Structural and Functional Biology
- `struct_bio`: secondary structure, domain/feature maps, simple folding interfaces
- `pathway`: pathway graphs, reaction networks, flux-style abstractions
- `network`: gene regulatory / interaction network primitives
- `expr`: expression matrices, normalization helpers, differential-expression building blocks

#### v0.5 Clinical and Translational Biology
- `clinical`: sample/cohort metadata, phenotype associations, biomarker utilities
- `variant`: annotation layers, effect models, consequence grouping
- `diagnostics` helpers for reporting and quality checks
- `evidence` model for linking findings to sources and confidence

#### v0.6 Broader Biology Coverage
- `cell`: cell-state and lineage abstractions
- `microbio`: microbial and strain-level features
- `imaging`: bioimage metadata and analysis adapters
- `ecology`: community and ecosystem abstractions
- `taxonomy`: classification and naming helpers
- `education`/`demo` examples only if they stay out of core API

---

## 6. Out of Scope (v0.1)

- ฐานข้อมูล organism/genome (`scies-bio-db`) — เป็นโปรเจกต์แยกในอนาคต จะ ingest ผ่าน type ของ lib นี้
- Wet-lab protocol, synthesis order, HDR template builder (`crispr_repair`)
- Structural biology (`struct_bio`), phylogenetics (`phylo`), population genetics (`popgen`),
  epidemiology (`epi`), gene expression (`expr`), variant annotation, omics modules, clinical modules — ตามด้วยใน v0.2+
- GPU/SIMD acceleration (พิจารณาเมื่อ core API นิ่งแล้ว)

---

## 7. Roadmap ถัดจาก v0.1

| Version | Scope |
|---|---|
| v0.1 | `seq`, `io`, `align`, `codon`, `edit::crispr/offtarget/primer` |
| v0.2 | `hmm`, `motif`, `annotation`, `variant` primitives, `genome` intervals, `interop` กับ `scies-math` |
| v0.3 | `genomics`, `transcriptomics`, `proteomics`, `metagenomics`, `phylo`, `popgen` |
| v0.4 | `struct_bio`, `pathway`, `network`, `expr` |
| v0.5 | `clinical`, `diagnostics`, `evidence`, `taxonomy`, `cell`, `microbio`, `imaging`, `ecology` |
| v0.6+ | `scies-bio-db` (โปรเจกต์แยก) — ดึงจาก NCBI/Ensembl/UniProt, ใช้ type จาก v0.1-v0.5 เป็น core |
