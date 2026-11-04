# scies-bio-th — v0.8 Design Document
**Phase A Foundation Hardening: Aging, Anatomy, Real-Genome Guides & Clinical Pharmacodynamics**

---

## 1. Context and Architectural Stance

`scies-bio-th` v0.8 delivers **Phase A Foundation Hardening**, establishing a mathematically sound, biophysically grounded substrate for synthetic biology, whole-body digital twins (`humanz`), and longevity intervention design.

Adhering strictly to our core mandate of **Zero External Dependencies in the Default Build (100% Pure Rust)**, v0.8 hardens four key biological modules:
1. **`aging`**: Strictly binds epigenetic rejuvenation efficacy to somatic lineage identity retention, couples all design reports with standardized autonomous safety guards, and prevents unrealistic zero-age clock resets.
2. **`anatomy`**: Provides complete biomarker sets, aging decay rates, and intervention response coefficients for 6 core organ systems used by `humanz`; introduces a synchronized cross-organ coupling engine (SASP circulation, vascular oxygenation, immune surveillance clearance); and anchors 8 clinical longevity biomarkers against published human reference ranges.
3. **`real_genome` & `guide`**: Implements end-to-end guide RNA design directly from real reference genome assemblies (e.g. T2T-CHM13, GRCh38), tracking provenance metadata (assembly, chrom, coordinates, cryptographic SHA-256 file hash), scanning SpCas9 and Cas12a PAMs, and enforcing hard biophysical constraints (poly-T, poly-G, GC 30–75%, telomeric repeat exclusion).
4. **`pharma`**: Expands pharmacokinetic/pharmacodynamic modeling to two-compartment kinetics with central/peripheral distribution volumes, hybrid half-lives, Hill-equation $E_{\max}$/$\text{EC}_{50}$/$\gamma$ dynamics, arbitrary intermittent dosing schedules, aging hallmark lever modulation, and automated therapeutic window decision rules.

---

## 2. Module Tree

```text
scies_bio_th::
├── (v0.1–v0.7) seq, io, align, codon, edit, genome, annotation, variant,
│               motif, hmm, phylo, popgen, genomics, transcriptomics,
│               proteomics, metagenomics, structure, network, stats,
│               docking, epigenome, imaging, ivg
│
├── seq
│   └── sha256         — Dependency-free streaming SHA-256 (FIPS 180-4) hasher
│
├── safety
│   ├── audit          — InterventionSafetyCoupling payload & GuardAction enum
│   └── checkpoint     — p53 DNA damage threshold and oncogenic lockouts
│
├── anatomy
│   ├── anchors        — 8 Clinical reference ranges (CKD-EPI eGFR, AHA BP, ACSM VO2max, ADA glucose, Flow-FISH telomeres, SASP, thymus, skin)
│   ├── coupling       — Synchronized cross-organ coupling (systemic SASP, O2 delivery, immune clearance)
│   ├── organ          — 20 Gray's Anatomy TA2 organ classes with OrganBiomarkerSet, OrganAgingRateParams, & OrganInterventionResponseCoefficients
│   └── cycles         — 7 Closed-loop whole-body systemic physiological circuits
│
├── aging
│   ├── reprogramming  — Identity-competence formula: P = D * f(lineage, chromatin)
│   ├── guide          — SpCas9 (NGG) & Cas12a (TTTV) real guide design, on-target Doench RS2 score, local CFD, hard constraint filters
│   ├── real_genome    — Real FASTA loader with RealLocusMetadata, provenance, & SHA-256 checksums
│   ├── clock          — Horvath 353-CpG epigenetic clock & multi-generation Gompertz aging
│   └── design         — High-fidelity longevity genome design report with InterventionSafetyCoupling
│
└── pharma
    └── pkpd           — Two-compartment PK (V1, V2, kel, k12, k21), Hill PD (Emax, EC50, gamma), intermittent dosing, aging levers & dosing actions
```

---

## 3. Mathematical Specifications

### 3.1 Lineage Identity Retention & Rejuvenation Competence Binding
Epigenetic partial reprogramming (OSK) cannot produce biological rejuvenation if the target cell loses lineage fidelity or undergoes chromatin destabilization:
$$\text{competence} = \text{lineage\_integrity} \times \left(1.0 - \frac{1.0}{1.0 + \exp(12.0 \times (\text{chromatin\_stability} - 0.50))}\right)$$
$$\text{productive\_dosage} = \text{nominal\_dosage} \times \text{competence}$$
If $\text{lineage\_integrity} \to 0$ or $\text{chromatin\_stability} < 0.35$, $\text{competence} \to 0$, rendering nominal OSK doses biologically unproductive and mitigating teratoma risk.

### 3.2 Synchronized Cross-Organ Coupling Dynamics
Inter-systemic physiological feedback couples the 6 core organ compartments:
1. **Systemic SASP Circulation**:
   $$\text{SASP}_{\text{sys}} = \sum_{i=1}^6 w_i \cdot \text{SASP}_i \quad (w_{\text{liver}}=0.30, w_{\text{art}}=0.20, w_{\text{skin}}=0.15, w_{\text{marrow}}=0.15, w_{\text{thym}}=0.10, w_{\text{brain}}=0.10)$$
2. **Vascular Oxygen Delivery Modulation**:
   $$p\text{O}_{2,\text{organ}} = p\text{O}_{2,\text{intrinsic}} \times (0.40 + 0.60 \times \eta_{\text{O2}})$$
   $$\eta_{\text{O2}} = 0.60 \times \text{Reserve}_{\text{vasc}} + 0.40 \times \frac{\text{CO}}{5.4}$$
3. **Immune Surveillance Senescent Clearance**:
   $$\text{Competence}_{\text{immune}} = 0.55 \times \frac{\text{Output}_{\text{thymic}}}{12.0} + 0.45 \times \text{Viability}_{\text{marrow}}$$
   If $\text{Competence}_{\text{immune}} > 0.80$, senescent AGE burden in distal organs is actively cleared:
   $$\Delta \text{AGE} = -0.04 \times (\text{Competence}_{\text{immune}} - 0.80)$$

### 3.3 Two-Compartment Pharmacokinetics & Hill-Equation Pharmacodynamics
Linear two-compartment distribution with central ($V_1$) and peripheral ($V_2$) compartments:
$$\frac{dC_1}{dt} = -(k_{10} + k_{12}) C_1 + k_{21} C_2 \frac{V_2}{V_1}$$
$$\frac{dC_2}{dt} = k_{12} C_1 \frac{V_1}{V_2} - k_{21} C_2$$
$$\alpha, \beta = \frac{(k_{10} + k_{12} + k_{21}) \pm \sqrt{(k_{10} + k_{12} + k_{21})^2 - 4 k_{10} k_{21}}}{2}$$
$$t_{1/2,\alpha} = \frac{\ln 2}{\alpha}, \quad t_{1/2,\beta} = \frac{\ln 2}{\beta}$$
Pharmacodynamic Hill equation response:
$$E(C) = E_{\max} \frac{C^\gamma}{\text{EC}_{50}^\gamma + C^\gamma}, \quad \text{Toxicity}(C) = T_{\max} \frac{C^{\gamma_{\text{tox}}}}{\text{TD}_{50}^{\gamma_{\text{tox}}} + C^{\gamma_{\text{tox}}}}$$

### 3.4 Clinical Reference Intervals (8 Longevity Biomarkers)
1. **eGFR**: $90.0\text{--}130.0\text{ mL/min/1.73m}^2$ (CKD-EPI Stage 1)
2. **Systolic BP**: $105.0\text{--}125.0\text{ mmHg}$ (AHA/ACC Normotensive)
3. **$\text{VO}_2\max$**: $40.0\text{--}55.0\text{ mL/kg/min}$ (ACSM Young Adult)
4. **Fasting Glucose**: $70.0\text{--}99.0\text{ mg/dL}$ (ADA Euglycemic)
5. **2h OGTT Glucose**: $70.0\text{--}139.0\text{ mg/dL}$ (ADA Normal Glucose Tolerance)
6. **Telomere Length**: Age-stratified bands (Young Adult: $8,500\text{--}11,000\text{ bp}$; Middle: $6,500\text{--}8,500\text{ bp}$; Aged: $4,500\text{--}6,500\text{ bp}$)
7. **SASP / CRP Index**: $0.01\text{--}0.15\text{ index}$ (Youthful Low Inflammaging)
8. **Thymic Naive T Output**: $8.0\text{--}15.0 \times 10^7\text{ cells/day}$ (Pre-involution)

---

## 4. Verification Suite

| Metric | Target | Result |
| :--- | :--- | :--- |
| `cargo test --features aging` | 100% pass | **174 / 174 passed (0 failed)** |
| `cargo test --all-features` | 100% pass | **188 / 188 passed (0 failed)** |
| `cargo check --all-targets --features aging` | 0 warnings | **0 errors, 0 warnings** |
| `cargo test -p humanz --lib --release` | 100% pass | **23 / 23 passed (0 failed)** |
| `examples/design_from_real_genome.rs` | Builds & runs | **Verified** |
