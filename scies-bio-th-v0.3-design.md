# scies-bio-th — v0.3 Design Document (Omics & Comparative Biology)

> A high-performance bio toolkit for Rust.
> ส่วนขยายระดับระบบ (Systems & Comparative) สำหรับการวิเคราะห์โอมิกส์และวิวัฒนาการ

---

## 1. เป้าหมายของโปรเจกต์ใน v0.3

ใน v0.3 เราจะขยับจากการจัดการ Sequence/Feature ระดับเส้นตรง (1D) ไปสู่การวิเคราะห์ข้อมูลเชิงประชากร วิวัฒนาการ และข้อมูลโอมิกส์เบื้องต้น โดยยังคงรักษาเอกลักษณ์ **Zero-dependency** สำหรับ Core features

เป้าหมายหลักประกอบด้วย:
1. **Phylogenetics (`phylo`)**: การสร้างต้นไม้วิวัฒนาการ (Phylogenetic tree) พื้นฐานด้วย UPGMA และรองรับรูปแบบ Newick
2. **Population Genetics (`popgen`)**: คอร์สำหรับคำนวณความถี่อัลลีล และทดสอบ Hardy-Weinberg Equilibrium
3. **Genomics (`genomics`)**: โครงสร้างพื้นฐานสำหรับ De Bruijn Graph (k-mer counting, edges) สำหรับปูทางสู่ Assembly
4. **Transcriptomics (`transcriptomics`)**: ฟังก์ชันทางคณิตศาสตร์สำหรับหาค่า TPM (Transcripts Per Million) จาก Read counts
5. **Proteomics (`proteomics`)**: การจำลองการย่อยโปรตีนด้วยเอนไซม์ (Trypsin cleavage) และการคำนวณมวล (Mass calculation)
6. **Metagenomics (`metagenomics`)**: การคำนวณดัชนีความหลากหลายทางชีวภาพ (Alpha diversity: Shannon index)

---

## 2. Module Tree (v0.3 scope)

```text
scies_bio_th::
├── (เดิม) seq, io, align, codon, edit, genome, annotation, variant, motif, hmm
├── phylo           — Tree struct, Distance matrix, UPGMA, Newick format
├── popgen          — Allele frequency, Hardy-Weinberg calculation
├── genomics        — De Bruijn graph primitives (K-mer edges)
├── transcriptomics — Expression quantification (TPM)
├── proteomics      — Cleavage (Trypsin), Monoisotopic mass
└── metagenomics    — Diversity metrics (Shannon index)
```

---

## 3. ลำดับการ Implement

1. **`transcriptomics` & `metagenomics`**: โมดูลทางคณิตศาสตร์ที่เรียบง่ายที่สุด (TPM, Shannon Index)
2. **`proteomics`**: โมดูลย่อยโปรตีน และคำนวณมวล (อิงจาก `seq::Protein`)
3. **`popgen`**: โมดูลคำนวณ HWE จากความถี่จีโนไทป์
4. **`phylo`**: สร้างโครงสร้าง Node/Tree, Distance Matrix และอัลกอริทึม UPGMA พร้อม Export เป็น Newick
5. **`genomics`**: สร้าง K-mer counting และ Edge generation สำหรับ De Bruijn Graph

---

## 4. Feature Flags

ยังคงใช้โครงสร้างเดิม โดยสามารถเปิด feature `interop` เพื่อเชื่อมโยงกับ `scies-math` สำหรับฟังก์ชันทางสถิติที่ซับซ้อนขึ้นในอนาคตได้ (แต่สำหรับ v0.3 อัลกอริทึมพื้นฐานจะยังคงอยู่ในตัวเพื่อรักษา Zero-dependency)
