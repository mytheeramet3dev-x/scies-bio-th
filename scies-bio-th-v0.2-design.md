# scies-bio-th — v0.2 Design Document (Sequence Biology Expansion)

> A high-performance bio toolkit for Rust.
> ส่วนต่อขยาย (Sequence Biology Expansion) จากฐานที่วางไว้ใน v0.1

---

## 1. เป้าหมายของโปรเจกต์ใน v0.2

หลังจาก v0.1 ได้วางรากฐานของ Data types, IO, Alignment, Codon, และพื้นฐานการ Edit ไปแล้ว ใน **v0.2** จะมุ่งเน้นที่ **Sequence Biology Expansion** ซึ่งหมายถึงการวิเคราะห์ลวดลาย (Patterns), โมเดลทางสถิติ, และการระบุตำแหน่งบนจีโนม (Annotation/Variant)

เป้าหมายหลักของ v0.2 ประกอบด้วย:
1. **Motif & HMM:** การวิเคราะห์รูปแบบซ้ำๆ ของลำดับเบสหรือโปรตีน (เช่น Transcription Factor Binding Sites) ผ่าน Position Weight Matrix (PWM) และ Hidden Markov Models (HMM) พื้นฐาน
2. **Genomic Intervals & Annotation:** การจัดการระบบพิกัดบนจีโนม และโครงสร้างของยีน (Exon, Intron, UTR) รวมถึงการค้นหา Overlap ระหว่างช่วง
3. **Variant Representation:** การเป็นตัวแทนของการกลายพันธุ์พื้นฐาน (SNP, Indel) เพื่อเตรียมพร้อมรับข้อมูลทางการแพทย์
4. **Integration with `scies-math`:** เริ่มใช้ประโยชน์จาก `scies-math` ในการคำนวณทางสถิติและความน่าจะเป็นสำหรับ HMM และ Motif scoring

---

## 2. Module Tree (v0.2 scope)

โครงสร้างเดิมจะถูกขยายความสามารถ และเพิ่มโมดูลใหม่ดังนี้:

```text
scies_bio_th::
├── (เดิม) seq, io, align, codon, edit
├── io           — (เพิ่ม) GFF3/GTF parser สำหรับ annotation, VCF parser อย่างง่าย
├── motif        — Position Frequency Matrix (PFM), Position Weight Matrix (PWM), Sequence Scanning
├── hmm          — Profile HMM พื้นฐาน, Viterbi algorithm (หา path ที่น่าจะเป็นไปได้ที่สุด)
├── genome       — Interval, IntervalTree (สำหรับหา Overlap), Chromosome masking
└── variant      — Snp, Indel, Variant effect abstractions
```

---

## 3. Feature Flags

ใน v0.2 จะเริ่มดึง `scies-math` เข้ามามีบทบาทสำคัญมากขึ้น:

| Flag | Enables |
|---|---|
| `interop` | (มีตั้งแต่ v0.1) ใน v0.2 จะผูก `hmm` และ `motif` เข้ากับ Probability distributions และ Matrix operations จาก `scies-math` |

---

## 4. โครงสร้าง Cargo.toml (ปรับปรุงจาก v0.1)

```toml
[dependencies]
# ยังคงเป็น optional เพื่อรักษาความเป็น Zero-dependency ใน default build
scies-math = { version = "0.2", optional = true }
serde = { version = "1", features = ["derive"], optional = true }
```
*(กรณีที่ไม่เปิด `interop` โค้ดใน HMM/Motif จะใช้อัลกอริทึมทางคณิตศาสตร์แบบ Basic ที่เขียนขึ้นเองภายใน)*

---

## 5. ลำดับการ Implement

การทำงานยังคงยึดหลัก **Signature-first workflow** เช่นเดียวกับ v0.1:

1. **`genome`** — สร้าง Data structure สำหรับ `Interval` (start, end, strand) และอัลกอริทึมเช็ค Overlap (เช่น Interval Tree พื้นฐาน) เพื่อเป็นฐานให้โมดูลอื่น
2. **`annotation` & `io::gff`** — สร้าง Types สำหรับ Feature (Gene, Transcript, Exon) และตัวอ่านไฟล์ GFF3
3. **`variant` & `io::vcf`** — สร้าง Types แทน SNP และ Indel และตัวอ่านไฟล์ VCF อย่างง่าย (ยังไม่เจาะลึกข้อมูล genotype เต็มรูปแบบ)
4. **`motif`** — สร้าง Position Weight Matrix และระบบ Scan หาคะแนนบน `Dna` sequence
5. **`hmm`** — สร้าง Profile HMM เบื้องต้น และ Implement Viterbi algorithm สำหรับทำ Sequence Alignment แบบอิงสถิติ

แต่ละขั้นตอนจะต้องมี Doc-comment, Doctest และ `cargo bench` กำกับในจุดที่เป็น Hot path เสมอ

---

## 6. Out of Scope (v0.2)

- การประกอบจีโนม (Genome Assembly)
- การแมป Read ปริมาณมหาศาล (Short-read mapping / BWA-style)
- โครงสร้างโปรตีน 3 มิติ (PDB), Network Biology, Phylogenetics
- ฐานข้อมูลเต็มรูปแบบ (ปล่อยให้ `scies-bio-db` จัดการ)

---

## 7. Roadmap ถัดจาก v0.2

| Version | Scope |
|---|---|
| **v0.2** | `hmm`, `motif`, `genome` intervals, `annotation`, `variant` primitives |
| **v0.3** | Omics and Comparative Biology (`genomics`, `transcriptomics`, `proteomics`, `metagenomics`, `phylo`, `popgen`) |
| **v0.4** | Structural and Functional Biology (`struct_bio`, `pathway`, `network`, `expr`) |
| **v0.5+**| Clinical, Diagnostics, Cell, Ecology, Imaging |
