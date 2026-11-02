# scies-bio-th v0.5.0 – Real Genome Integration & High-Fidelity Longevity Design

## 1. Overview and Core Philosophy

Version **0.5.0** transitions `scies-bio-th` from abstract numerical simulations to **real-world telomere and genome analysis** directly compatible with assemblies such as **T2T-CHM13v2.0** and **GRCh38.p14**.

The primary goals of this release are:
1. **Real Chromosome Terminal Analysis**: Direct detection, quantification, and subtelomeric boundary transition mapping of real human telomeres (both 5' p-arm `(CCCTAA)_n` and 3' q-arm `(TTAGGG)_n` canonical and degenerate variant repeats).
2. **Selective / Indexed FASTA Streaming**: High-performance streaming of targeted chromosomal loci and terminal regions without allocating full gigabyte-scale genomic assemblies in RAM, maintaining **zero external dependencies** in the default build.
3. **High-Fidelity Simulation Models (`fidelity.rs`)**:
   - Telomere attrition driven by cell division, 8-oxo-dG oxidative guanine lesions, G-quadruplex replication stress, and Shelterin complex capping degradation.
   - Comprehensive multi-component Epigenetic/Biological Age estimation combining Horvath-like DNA methylation, telomere length, DNA double-strand break (DSB) foci ($\gamma$-H2AX / 53BP1), and SASP burden.
   - Realistic Telomerase kinetics modeling Repeat Addition Processivity (RAP), catalytic translocation, and shelterin-mediated steric accessibility.
4. **Pathway-Specific Oncogenic Transformation Risk (`risk.rs`)**:
   - Explicit quantification of neoplastic risk across 4 distinct failure modes:
     1. Telomerase-mediated crisis bypass without senescence.
     2. Accidental proto-oncogene transactivation.
     3. Checkpoint loss / tumor suppressor insufficiency.
     4. Dicentric chromosome formation via shelterin uncapping.
   - Multilayered genetic mitigations (Super-p53 gene duplication, inducible suicide switches, bi-allelic safeguards).
5. **Expanded Longevity Pathways (`pathways.rs`)**:
   - Added canonical longevity and maintenance genes: `ATM`, `ATR`, `PARP1`, `LMNA`, `MTOR`, `AMPK` (`PRKAA1`), alongside existing `TERT`, `TERC`, `POT1`, `SIRT6`, `SIRT1`, `FOXO3`, `CDKN2A`, `TP53`.
6. **Longevity Design Engine 2.0 & Structured Reporting (`report.rs`, `design.rs`)**:
   - Distinct, scientifically calibrated design strategies: `Conservative`, `Balanced`, and `Aggressive`.
   - Actionable execution priority ordering (Safety locks first $\rightarrow$ telomere stabilization $\rightarrow$ DNA repair enhancement $\rightarrow$ metabolic & senolytic modulation).
   - High-fidelity `LongevityDesignReport` with per-edit confidence scores, real genomic coordinates, target strand orientation, PAM loci, and candidate sgRNAs with CFD off-target risk assessment.

---

## 2. Module Architecture

```
src/
├── aging/
│   ├── mod.rs              // Public re-exports and module orchestration
│   ├── real_genome.rs      // Terminal chromosome loaders, repeat detection, subtelomere mapping
│   ├── fidelity.rs         // High-fidelity biophysical telomere & age simulation
│   ├── risk.rs             // Pathway-level oncogenic risk modeling & safety lock evaluation
│   ├── report.rs           // High-fidelity longevity report, confidence metrics, execution order
│   ├── pathways.rs         // Expanded longevity registry (14 canonical genes)
│   ├── design.rs           // Design Engine 2.0 producing prioritized CRISPR sgRNAs
│   ├── telomere.rs         // Core telomere states and backward-compatible simulation
│   ├── telomerase.rs       // Core telomerase activity and catalytic modeling
│   └── clock.rs            // Baseline biological age calculations
├── io/
│   ├── mod.rs              // Selective FASTA region scanning and buffered parsing
│   ├── fasta.rs / gff / vcf
```

---

## 3. High-Fidelity Mathematical Specifications

### 3.1 Telomere Shortening with Oxidative and Replication Stress
$$\Delta L = B \cdot \left(1.0 + \kappa_{\text{ox}} \cdot \text{ROS} + \kappa_{\text{g4}} \cdot \text{G4Stress}\right) \cdot \frac{1.0}{\max(0.2, \text{ShelterinIntegrity})}$$
Where:
- $B \approx 55\text{ bp/division}$ (basal end-replication problem loss).
- $\text{ROS} \in [0.0, 3.0]$: Reactive oxygen species producing 8-oxo-7,8-dihydroguanine.
- $\text{G4Stress} \in [0.0, 2.0]$: G-quadruplex stabilization impeding lagging strand synthesis.
- $\text{ShelterinIntegrity} \in [0.0, 1.0]$: TRF1/TRF2/POT1 complex density.

### 3.2 High-Fidelity Biological Age
$$\text{Age}_{\text{Bio}} = w_m \cdot M(\text{CpG}) + w_t \cdot T(L) + w_d \cdot D(\gamma\text{-H2AX}) + w_s \cdot S(\text{SASP}) + w_c \cdot \text{Age}_{\text{Chrono}}$$
Where:
- $w_m = 0.40$: DNA methylation epigenetic age component.
- $w_t = 0.25$: Telomeric attrition component.
- $w_d = 0.15$: Double-strand break DNA damage burden.
- $w_s = 0.10$: Senescence-Associated Secretory Phenotype burden.
- $w_c = 0.10$: Baseline chronological drift.

### 3.3 Oncogenic Risk Index and Mitigation
$$\text{Risk}_{\text{net}} = \left( \sum_{p \in \text{Pathways}} R_p \cdot w_p \right) \times \prod_{m \in \text{Mitigations}} (1.0 - E_m)$$
Where:
- $R_p$: Raw risk per pathway (Crisis bypass, Proto-oncogene transactivation, Checkpoint loss, Chromosome fusion).
- $E_m$: Efficacy factor of active safety locks (e.g. Super-p53: $0.65$, Inducible suicide switch: $0.70$).

---

## 4. Real Genome Integration Details

- **Terminal repeat recognition**:
  - Forward telomere: `TTAGGG` and variants (`TGAGGG`, `TCAGGG`, `TTGGGG`).
  - Reverse telomere: `CCCTAA` and variants (`CCCTCA`, `CCCTGA`, `CCCCAA`).
- **Subtelomeric Transition**:
  - Scanning algorithm calculates telomeric repeat density over sliding windows of 100 bp.
  - The boundary is detected when repeat density drops below 40%, demarcating the canonical telomere tract from the variable subtelomeric region.
- **Selective FASTA Extraction**:
  - Reader scans FASTA headers matching target sequence IDs or chromosomal identifiers and streams only relevant window slices without reading entire multimillion-base chromosomes into memory.

---

## 5. Verification and Quality Targets

- Default build remains pure `std` (Zero external crates).
- Full compatibility with existing v0.1-v0.4 APIs.
- 100% test pass rate across all modules.
- Zero clippy warnings with `--all-targets --all-features`.
- Working example `examples/design_from_real_genome.rs` demonstrating end-to-end execution.
