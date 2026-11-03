# scies-bio-th v0.6.0 – Real Guide Design, Multi-Generation Simulation & Hard Biological Constraints

## 1. Overview and Core Objectives

Version **0.6.0** of `scies-bio-th` introduces rigorous biological constraints and population dynamics modeling to the synthetic longevity genome engineering platform:

1. **Deterministic Real sgRNA Design (`guide.rs`)**:
   - Extraction of real 20-nt protospacer sequences from forward (`NGG`) and reverse (`CCN`) SpCas9 PAM sites on actual locus sequence contexts.
   - Stringent oligo synthesis and transcription filter heuristics:
     - Strict Pol III U6 transcription termination avoidance: rejects or heavily penalizes $\ge 4$ consecutive Thymines/Uracils (`TTTT` / `UUUU`).
     - G-quadruplex oligo avoidance: filters $\ge 4$ consecutive Guanines (`GGGG`).
     - GC content optimization within $40\% - 70\%$ range.
     - Avoidance of direct terminal telomeric repeat tracts (`TTAGGG` / `CCCTAA`) where DSBs cause terminal arm loss or dicentric end fusions.
   - Integrated CFD off-target penalty assessment and multi-parameter engineering confidence score.

2. **Multi-Generation Cell Population Simulation (`simulation.rs`)**:
   - Long-term stochastic/deterministic population modeling over 20–100 cell division generations post-editing.
   - Dynamic tracking of:
     - Average chromosomal telomere length (bp).
     - Epigenetic and multi-factorial biological age (years).
     - Cell population partitioning across 3 distinct compartments:
       1. **Proliferating Viable Fraction**: dividing cells maintaining homeostasis.
       2. **Replicative Senescent Fraction**: arrest triggered at Hayflick limit ($\le 3,800\text{ bp}$) or DNA damage foci.
       3. **Malignantly Transformed Fraction**: neoplastic clones bypassing M1/M2 crisis.
     - Cumulative oncogenic transformation probability over organismal life.

3. **Hard Biological Constraints & Safety Enforcement (`constraints.rs`)**:
   - **Constraint 1 (TERT Proscription Without Lock)**: Strict prohibition of TERT activation unless coupled with Pachyderm Super-p53 tandem duplication and inducible suicide switches (iCasp9).
   - **Constraint 2 (Checkpoint Primacy)**: Safety locks and tumor suppressor enhancements (TP53 / CDKN2A) must strictly precede or execute in Step 1 alongside any proliferative reactivation.
   - **Constraint 3 (Strategy Conformance)**: Conservative strategy strictly forbids TERT expression; Balanced restricts expression to physiological stem-cell levels with single/dual locks; Aggressive demands multi-layered active locks.
   - **Constraint 4 (Terminal Telomere DSB Guard)**: Rejects any guide cutting within the pure terminal telomere tract.
   - Automatic remediation and audit violation reporting (`ConstraintViolation`).

4. **Expanded Real Locus Genomic Datasets (`testdata/`)**:
   - Complete T2T-CHM13/GRCh38 canonical sequence files for `TP53`, `CDKN2A`, `FOXO3`, `ATM`, `TERT`, and `SIRT6`.

---

## 2. Mathematical and Simulation Framework

### 2.1 Multi-Generation Telomere and Cell Partitioning
For generation $g \in [1, N]$:
$$\Delta L_{\text{loss}}(g) = B_0 \cdot \left(1.0 + \kappa_{\text{ox}} \cdot \text{ROS} + \kappa_{\text{g4}} \cdot \text{G4}\right) \cdot \frac{1.0}{\max(0.2, \text{Shelterin}(g))}$$
$$\Delta L_{\text{gain}}(g) = \text{Expression} \cdot \text{Processivity} \cdot \text{Accessibility}(g)$$
$$L(g) = \min\left(18,500, \max\left(600, L(g-1) + \Delta L_{\text{gain}}(g) - \Delta L_{\text{loss}}(g)\right)\right)$$

### 2.2 Oncogenic Transformation and Population Compartments
Let $P_{\text{crisis}}(g)$ be the probability of crisis bypass:
$$P_{\text{trans}}(g) = \left(\text{Expression}^{1.4} \cdot 0.006\right) \cdot D_{\text{p53}} \cdot D_{\text{suicide}}$$
Where:
- $D_{\text{p53}} = 0.15$ if Super-p53 is active, else $1.00$.
- $D_{\text{suicide}} = 0.20$ if inducible suicide switch is armed, else $1.00$.

The cellular compartments evolve according to:
$$F_{\text{senescent}}(g) = \begin{cases} 
F_{\text{senescent}}(g-1) + (1.0 - F_{\text{senescent}}(g-1)) \cdot 0.08 & \text{if } L(g) \le 3,800 \\
F_{\text{senescent}}(g-1) \cdot (1.0 - \text{SenolyticRate}) & \text{if } L(g) > 3,800 
\end{cases}$$
$$F_{\text{transformed}}(g) = F_{\text{transformed}}(g-1) + \left(1.0 - F_{\text{senescent}}(g) - F_{\text{transformed}}(g-1)\right) \cdot P_{\text{trans}}(g)$$
$$F_{\text{proliferating}}(g) = \max\left(0.0, 1.0 - F_{\text{senescent}}(g) - F_{\text{transformed}}(g)\right)$$

---

## 3. Structural Roadmap

```text
src/aging/
├── guide.rs          // Real sgRNA design, PAM scanning, poly-run & CFD scoring
├── simulation.rs     // Multi-generation cell population & trajectory simulation
├── constraints.rs    // Biological safety constraints & auto-remediation
├── validation.rs     // Pre-flight design verification & audit
├── real_genome.rs    // Real chromosome & locus loaders
├── fidelity.rs       // High-fidelity biophysical kinetics
├── risk.rs           // Multi-pathway oncogenic hazard modeling
├── report.rs         // Comprehensive report with simulation trajectory & audit
├── pathways.rs       // 14 canonical longevity gene registry
├── design.rs         // Design Engine 2.0 with constraint enforcement
├── telomere.rs       // Core telomere states
├── telomerase.rs     // Core telomerase kinetics
└── clock.rs          // Biological age calculation
```
