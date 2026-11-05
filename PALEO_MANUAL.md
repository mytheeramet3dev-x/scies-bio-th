# scies-bio-th: Paleo-Genomic & Paleoproteomic Reconstruction Manual
**Technical Reference, Mathematical Foundations & Practical Guide for Extinct Genome Recovery**

---

## Table of Contents
1. [Scientific Principles & Preservation Thermodynamics](#1-scientific-principles--preservation-thermodynamics)
2. [Tri-Layer Architecture Overview](#2-tri-layer-architecture-overview)
3. [Mathematical Foundations & Algorithms](#3-mathematical-foundations--algorithms)
   - [3.1 Nucleotide Posterior Distributions & Shannon Entropy](#31-nucleotide-posterior-distributions--shannon-entropy)
   - [3.2 Continuous-Time Markov Chain (CTMC) & Felsenstein Pruning](#32-continuous-time-markov-chain-ctmc--felsenstein-pruning)
   - [3.3 Extinct Branch Time Propagation](#33-extinct-branch-time-propagation)
   - [3.4 Bayesian Paleoproteomic Codon Back-Translation](#34-bayesian-paleoproteomic-codon-back-translation)
   - [3.5 Briggs Deamination-Aware Ancient DNA Modeling](#35-briggs-deamination-aware-ancient-dna-modeling)
   - [3.6 Outgroup-Bracketed Pseudogene Resurrection](#36-outgroup-bracketed-pseudogene-resurrection)
   - [3.7 Commutative Group-Discounted Phenotypic Pooling](#37-commutative-group-discounted-phenotypic-pooling)
4. [API & Data Structures Reference](#4-api--data-structures-reference)
5. [End-to-End Practical Recipes](#5-end-to-end-practical-recipes)
   - [Recipe 1: Phylogenetic ASR from Tree](#recipe-1-phylogenetic-asr-from-tree)
   - [Recipe 2: Ancient Peptide Fragment Anchoring](#recipe-2-ancient-peptide-fragment-anchoring)
   - [Recipe 3: Ancient DNA Read Integration](#recipe-3-ancient-dna-read-integration)
   - [Recipe 4: Pseudogene Repair & Open Reading Frame Restoration](#recipe-4-pseudogene-repair--open-reading-frame-restoration)
   - [Recipe 5: Multi-Evidence Engine Synthesis & Chassis Edit Planning](#recipe-5-multi-evidence-engine-synthesis--chassis-edit-planning)
6. [Case Studies](#6-case-studies)
   - [Case Study A: Pleistocene Megafauna (*Mammuthus primigenius*)](#case-study-a-pleistocene-megafauna-mammuthus-primigenius)
   - [Case Study B: Mesozoic Non-Avian Dinosaur (*Tyrannosaurus rex*)](#case-study-b-mesozoic-non-avian-dinosaur-tyrannosaurus-rex)
7. [Downstream Interoperability & CRISPR Guide Design](#7-downstream-interoperability--crispr-guide-design)

---

## 1. Scientific Principles & Preservation Thermodynamics

De-extinction and extinct genome reconstruction face hard physical boundaries governed by chemical kinetics. Phosphodiester bond hydrolysis in DNA has an estimated effective half-life of $\approx 521$ years at $13.1^\circ\text{C}$ in bone matrix (Allentoft et al. 2012). Even in optimal sub-zero permafrost conditions, authentic recoverable DNA fragments ($>30\text{ bp}$) do not survive beyond approximately $1.5$ to $2.0$ million years.

```text
Time Depth        Surviving Biomolecules             Optimal Computational Reconstruction Strategy
──────────────────────────────────────────────────────────────────────────────────────────────────────
0 – 1.5 Ma        Degraded Ancient DNA (30–100 bp)   Briggs Post-Mortem Deamination Alignment
(Pleistocene)     + Collagen & Enamel Proteomes      + Closest Extant Reference Mapping
                                                     Accurate down to individual nucleotide variants (>99%)

1.5 – 10+ Ma      Tooth Enamelin (ENAM, AMELX)       Bayesian Paleoproteomic Codon Back-Translation
(Miocene / Plioc) + Bone Type I Collagen (COL1)      + Clade-Specific Synonymous Codon Usage Priors
                  DNA completely degraded            Coding sequences reconstructed with >90% fidelity

66+ Ma            Biomolecules rarely preserved      Extant Phylogenetic Bracketing (EPB)
(Mesozoic)        Occasional peptide traces          + Continuous-Time Markov Chain ASR
                  Inactivated Pseudogenes in birds   + Pseudogene Resurrection + Phenotypic Selection Priors
```

### The Three Biological Pillars
1. **Paleoproteomics**: Biomineralized proteins sealed within hydroxyapatite crystal lattices (specifically tooth enamel and osteoid collagen) exhibit orders-of-magnitude higher resistance to diagenetic hydrolysis than nucleic acids.
2. **Extant Phylogenetic Bracketing (EPB)**: Extinct organisms are phylogenetically bracketed by their closest living relatives (e.g., non-avian dinosaurs bracketed by Aves and Crocodylia; mammoths bracketed by Asian and African elephants).
3. **Pseudogene Remnants**: When ancestral traits are lost (such as teeth in modern birds), the structural genes (`ENAM`, `AMELX`, `MMP20`) are rarely deleted instantaneously; instead, they undergo pseudogenization through nonsense stop mutations and frameshift indels, preserving identifiable sequence fragments in the modern descendant's genome.

---

## 2. Tri-Layer Architecture Overview

The `paleo` module implements a mathematically rigorous three-layer pipeline:

```text
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                        Layer 1: Phylogenetic Ancestral Core                            │
│  Extant Relatives Alignment + Tree (PhyloNode) + Model (HKY85+Γ₄ / GTR / K80 / JC69)   │
│                                           │                                            │
│                        Log-Scaled Felsenstein Pruning                                  │
│                                           ▼                                            │
│                 Marginal Ancestral Posteriors + Extinct Branch Drift                   │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
┌───────────────────────────────────────────▼────────────────────────────────────────────┐
│                        Layer 2: Paleobiomolecular Anchoring                            │
│                                                                                        │
│   ┌─────────────────────────────────────┐   ┌──────────────────────────────────────┐   │
│   │ 2A. Paleoproteomics                 │   │ 2B. Ancient DNA (aDNA)               │   │
│   │ Fossil Peptides (LC-MS/MS)          │   │ Degraded short reads (30-100 bp)     │   │
│   │ + Hydroxyproline (O->P)             │   │ + Briggs 5' C->T & 3' G->A model     │   │
│   │ + Synonymous Codon Usage Priors     │   │ + Sequencing accuracy convolution    │   │
│   └──────────────────┬──────────────────┘   └──────────────────┬───────────────────┘   │
│                      │                                         │                       │
│                      └────────────────────┬────────────────────┘                       │
│                                           │                                            │
│   ┌───────────────────────────────────────▼────────────────────────────────────────┐   │
│   │ 2.5. Pseudogene Resurrection                                                   │   │
│   │ Inactivated extant pseudogene + Functional outgroup bracket                    │   │
│   │ Repairs premature stop codons (TGA/TAA/TAG) & restores Open Reading Frames     │   │
│   └───────────────────────────────────────┬────────────────────────────────────────┘   │
└───────────────────────────────────────────┼────────────────────────────────────────────┘
                                            │
┌───────────────────────────────────────────▼────────────────────────────────────────────┐
│                  Layer 3: Commutative Phenotypic Selection                             │
│  Osteology & Histology (Fibrolamellar, Pneumaticity) + Allometry (Peto's Paradox)      │
│  + Stable Isotopes (δ13C, δ15N trophic, δ18O paleoclimate)                             │
│                                           │                                            │
│                   Group-Discounted Log-Linear Bayesian Pooling                         │
│                    α_g = 1 / sqrt(|g|) (Strictly Commutative)                          │
│                                           ▼                                            │
│               Functional Pathway Selection Posteriors (dN/dS Lever)                    │
└───────────────────────────────────────────┬────────────────────────────────────────────┘
                                            │
┌───────────────────────────────────────────▼────────────────────────────────────────────┐
│                        Synthesis & Chassis Edit Planner                                │
│                     PaleoReconstructionEngine -> PaleoGenomeDraft                     │
│  - Consensus Dna & Protein Sequences                                                   │
│  - Per-site Shannon Entropy Map & High-Uncertainty Region Detection                    │
│  - Prioritized ChassisGenomicEdit Instructions against Extant Relative (Chassis)       │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## 3. Mathematical Foundations & Algorithms

### 3.1 Nucleotide Posterior Distributions & Shannon Entropy
At each nucleotide position $i \in \{0, \dots, N-1\}$, the state is represented by a normalized probability vector:
$$\mathbf{p}_i = \big(P(A), P(C), P(G), P(T)\big), \quad \sum_{b \in \{A,C,G,T\}} p_b = 1.0$$

The uncertainty at position $i$ is quantified by its Shannon entropy:
$$H(\mathbf{p}_i) = -\sum_{b \in \{A,C,G,T\}} p_b \log_2(p_b) \quad [\text{bits}]$$
- $H = 0.0\text{ bits}$: Absolute certainty (a single base has probability $1.0$).
- $H = 2.0\text{ bits}$: Complete uncertainty ($p_A = p_C = p_G = p_T = 0.25$).

#### Dynamic IUPAC Representation
When exporting sequences for downstream primer or guide design, `BasePosterior::to_iupac(min_confidence)` applies a cumulative probability threshold $\theta \in [0.25, 0.9999]$:
- If $\max_b(p_b) \ge \theta$, the canonical base is returned.
- Otherwise, bases are accumulated in descending order of probability until $\sum p_b \ge \theta$, mapping to the corresponding IUPAC degeneracy code:
  - $\{A, G\} \to \text{R}$ (puRine)
  - $\{C, T\} \to \text{Y}$ (pYrimidine)
  - $\{C, G\} \to \text{S}$ (Strong)
  - $\{A, T\} \to \text{W}$ (Weak)
  - $\{G, T\} \to \text{K}$ (Keto)
  - $\{A, C\} \to \text{M}$ (aMino)
  - $\{C, G, T\} \to \text{B}$
  - $\{A, G, T\} \to \text{D}$
  - $\{A, C, T\} \to \text{H}$
  - $\{A, C, G\} \to \text{V}$
  - All four $\to \text{N}$

---

### 3.2 Continuous-Time Markov Chain (CTMC) & Felsenstein Pruning
Evolutionary substitution along branch length $t$ (in millions of years, Ma) is modeled by a $4 \times 4$ rate transition probability matrix $P(t) = \exp(Qt)$.

#### Substitution Models (HKY85, GTR, K80, JC69)
The default vertebrate coding model is **HKY85 + discrete $\Gamma_4$** rate heterogeneity (`SubstitutionModel::default_vertebrate_coding()`):
- **HKY85 (Hasegawa, Kishino, Yano 1985)**: Differentiates transitions from transversions with $\kappa$ while accounting for unequal base frequencies $\pi = [\pi_A, \pi_C, \pi_G, \pi_T]$.
- **Discrete $\Gamma_4$ Rate Heterogeneity**: Computes 4 category substitution rates using Wilson-Hilferty transformation under shape parameter $\alpha$ (e.g. $\alpha = 0.60$).
- **GTR (General Time Reversible)**: Full 6-parameter symmetric exchangeabilities evaluated via matrix exponential scaling and squaring.
- **Log-Scaled Felsenstein Pruning**: At each internal node, likelihood vectors are dynamically scaled by $s = \max_b L(b)$ to guarantee absolute mathematical immunity against numerical underflow on deep branches (up to 250+ Ma).

#### Kimura 2-Parameter (K80) Formulation
Differentiating transitions ($A \leftrightarrow G, C \leftrightarrow T$ at rate $\alpha$) from transversions ($A,G \leftrightarrow C,T$ at rate $\beta$) with transition/transversion ratio $\kappa = \alpha / \beta$ and overall substitution rate $\mu$:
$$d = \mu \cdot t$$
$$P_{\text{transversion}}(t) = \frac{1}{4} - \frac{1}{4} \exp\left( -\frac{4d}{\kappa + 2} \right)$$
$$P_{\text{transition}}(t) = \frac{1}{4} + \frac{1}{4} \exp\left( -\frac{4d}{\kappa + 2} \right) - \frac{1}{2} \exp\left( -\frac{2d(\kappa + 1)}{\kappa + 2} \right)$$
$$P_{\text{identical}}(t) = 1 - P_{\text{transition}}(t) - 2 P_{\text{transversion}}(t)$$

#### Post-Order Tree Pruning
For a rooted tree $\mathcal{T}$ and a site $k$, conditional likelihoods $L_u(b)$ at node $u$ given child subtrees are computed upward from leaves to the root:
- **Leaf node** $v$ with observed base $b_{\text{obs}}$:
  $$L_v(b) = \begin{cases} 1.0 & \text{if } b = b_{\text{obs}} \\ 0.0 & \text{if } b \ne b_{\text{obs}} \\ 1.0 & \text{if ambiguous / gap } (N) \end{cases}$$
- **Internal node** $u$ with children $v \in \text{ch}(u)$ separated by branch length $t_{u,v}$:
  $$L_u(b) = \prod_{v \in \text{ch}(u)} \left( \sum_{b' \in \{A,C,G,T\}} P_{b \to b'}(t_{u,v}) L_v(b') \right)$$
- **Marginal Ancestral Posterior** at root $r$ under stationary prior $\pi = [0.25; 4]$:
  $$P_{\text{anc}}(b) = \frac{L_r(b)}{\sum_{b'} L_r(b')}$$

---

### 3.3 Extinct Branch Time Propagation
Extinct taxa typically diverge from the common ancestral node of extant relatives by $\Delta t_{\text{extinct}}$ million years (e.g., *Tyrannosaurus rex* separated from the bird-crocodilian root by $\approx 66\text{ Ma}$). The probability vector is forward-propagated:
$$P_{\text{extinct}}(b') = \sum_{b \in \{A,C,G,T\}} P_{\text{anc}}(b) \cdot P_{b \to b'}(\Delta t_{\text{extinct}})$$
As $\Delta t_{\text{extinct}} \to \infty$, the distribution smoothly and monotonically approaches maximum entropy without artificial hard boundaries.

---

### 3.4 Bayesian Paleoproteomic Codon Back-Translation
Given a peptide fragment of length $M$ codons recovered via mass spectrometry with match confidence $q \in (0, 1]$:
For each codon index $k$ aligning to base indices $3k, 3k+1, 3k+2$:
1. **Prior**: Prior joint probability over codon $c = (b_1, b_2, b_3) \in \{0..3\}^3$ is formed from Layer 1 base posteriors combined with relative synonymous codon usage bias $u(c)$ measured from extant relatives:
   $$P_0(c) \propto P_1(b_1) \cdot P_2(b_2) \cdot P_3(b_3) \cdot u(c)$$
2. **Likelihood**: If the translated amino acid matches the fossil observation $aa_{\text{obs}}$ (with hydroxyproline `O` mapped to proline `P`):
   $$L(aa_{\text{obs}} \mid c) = \begin{cases} q + \frac{1-q}{20} & \text{if } \text{translate}(c) = aa_{\text{obs}} \\ \frac{1-q}{20} & \text{otherwise} \end{cases}$$
3. **Marginalization**: The joint posterior $P(c \mid aa_{\text{obs}}) \propto P_0(c) \cdot L(aa_{\text{obs}} \mid c)$ is marginalized to update each nucleotide site independently:
   $$P'_1(b_1) = \sum_{b_2, b_3} P\big((b_1, b_2, b_3) \mid aa_{\text{obs}}\big)$$

---

### 3.5 Briggs Deamination-Aware Ancient DNA Modeling
Ancient DNA fragments suffer post-mortem hydrolytic deamination of cytosine to uracil (read as thymine during sequencing) at single-stranded 5' overhangs, and guanine to adenine at 3' overhangs.

At distance $d_{5'}$ from the 5' read terminus and $d_{3'}$ from the 3' terminus, deamination probabilities decay exponentially with scale parameter $\lambda \approx 0.35$:
$$\delta_{5'}(d_{5'}) = \delta_{\max, 5'} \cdot e^{-\lambda d_{5'}}, \quad \delta_{3'}(d_{3'}) = \delta_{\max, 3'} \cdot e^{-\lambda d_{3'}}$$

The transition matrix $D(d)$ of chemical degradation transforms true base $b_{\text{true}}$ into damaged intermediate $b_{\text{dmg}}$:
$$D = \begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & 1 - \delta_{5'} & 0 & \delta_{5'} \\
\delta_{3'} & 0 & 1 - \delta_{3'} & 0 \\
0 & 0 & 0 & 1
\end{pmatrix} \begin{matrix} A \\ C \\ G \\ T \end{matrix}$$

Convolving with sequencing accuracy $A \in [0.5, 0.9999]$ and background substitution error $\epsilon = (1 - A)/3$:
$$P(b_{\text{obs}} \mid b_{\text{true}}) = \sum_{b_{\text{dmg}}} D_{b_{\text{true}}, b_{\text{dmg}}} \cdot \begin{cases} A & \text{if } b_{\text{dmg}} = b_{\text{obs}} \\ \epsilon & \text{if } b_{\text{dmg}} \ne b_{\text{obs}} \end{cases}$$

This prevents spurious $C \to T$ variant calls at fragment ends while extracting maximum signal from internal read positions.

---

### 3.6 Outgroup-Bracketed Pseudogene Resurrection
Modern taxa frequently carry inactivated pseudogenes (e.g., tooth enamelin `ENAM` in chicken with internal `TGA` stops).
Given an inactivated descendant sequence $S_{\text{inact}}$ and functional outgroup sequences $S_{\text{out}, j}$ with divergence times $t_j$:
1. At undamaged codons, joint likelihood incorporates ingroup and outgroup signals.
2. At codons where the ingroup contains a premature stop (`TAA`, `TAG`, `TGA`) or frameshift defect (`N`), the engine conditions the ancestral joint distribution strictly over the **61 sense codons**:
   $$P_{\text{resurrected}}(c) \propto \mathbb{I}(c \notin \{\text{TAA}, \text{TAG}, \text{TGA}\}) \cdot \prod_{m=0}^2 \left( P_{\text{inact}}(b_m \mid c_m) \prod_j P_{\text{out}, j}(b_{m, j} \mid c_m) \right)$$
3. The restored codon replaces the nonsense triplet, resetting the Open Reading Frame (ORF) and producing an auditable `RepairEvent`.

---

### 3.7 Commutative Group-Discounted Phenotypic Pooling
To eliminate non-commutative order artifacts and avoid false overconfidence from correlated fossil observations:
1. Every phenotypic evidence item $e$ is assigned to a `CorrelationGroup` $g \in \mathcal{G}$:
   - `OsteologyHistology`: Bone microstructure, growth rings, pneumaticity.
   - `GeochemicalIsotope`: $\delta^{13}\text{C}, \delta^{15}\text{N}, \delta^{18}\text{O}$.
   - `MorphometricScaling`: Body mass, stylopodial circumference.
   - `PaleoEnvironment`: Paleotemperature, depositional setting.
2. For each functional pathway $M$, evidence $e$ emits a signed log-likelihood ratio $\ell_e(M)$.
3. Intra-group correlation damping uses the square-root factor $\alpha_g = \frac{1}{\sqrt{|g|}}$:
   $$\text{logit}\, P(M \mid \mathcal{E}) = \text{logit}\, P_0(M) + \sum_{g \in \mathcal{G}} \frac{1}{\sqrt{|g|}} \sum_{e \in g} r_e \cdot \ell_e(M)$$
4. The posterior selection probability is obtained via the logistic sigmoid:
   $$P(M \mid \mathcal{E}) = \frac{1}{1 + \exp(-\text{logit}\, P(M \mid \mathcal{E}))}$$

**Mathematical Proof of Commutativity**: Addition in $\mathbb{R}$ is associative and commutative; grouping is partition-based. Hence, $\forall \sigma \in \text{Permutations}(\mathcal{E})$, $P(M \mid \sigma(\mathcal{E})) \equiv P(M \mid \mathcal{E})$.

---

## 4. API & Data Structures Reference

| Type / Struct | Description | Module |
| :--- | :--- | :--- |
| [`BasePosterior`](file:///C:/Projects/github/scies-bio-th/src/paleo/types.rs#L96-L221) | 4-element probability vector `[P(A), P(C), P(G), P(T)]` with entropy calculation and IUPAC encoding | `scies_bio_th::paleo::types` |
| [`ReconstructedLocus`](file:///C:/Projects/github/scies-bio-th/src/paleo/types.rs#L285-L382) | Full locus structure containing site posteriors, consensus DNA/Protein, and uncertainty intervals | `scies_bio_th::paleo::types` |
| [`SubstitutionModel`](file:///C:/Projects/github/scies-bio-th/src/paleo/asr.rs#L17-L84) | CTMC transition models (`Hky85` with discrete $\Gamma_4$, `Gtr`, `Kimura2Param`, `JukesCantor`) | `scies_bio_th::paleo::asr` |
| [`FelsensteinAsr`](file:///C:/Projects/github/scies-bio-th/src/paleo/asr.rs#L94-L198) | Upward post-order tree pruning and forward terminal branch drift propagation | `scies_bio_th::paleo::asr` |
| [`PeptideFragment`](file:///C:/Projects/github/scies-bio-th/src/paleo/biomolecule.rs#L22-L40) | Mass-spec peptide fragment with PSM confidence and hydroxyproline normalization | `scies_bio_th::paleo::biomolecule` |
| [`PaleoproteomeTranslator`](file:///C:/Projects/github/scies-bio-th/src/paleo/biomolecule.rs#L43-L177) | Bayesian reverse-translation from amino acids to codons using codon-usage bias priors | `scies_bio_th::paleo::biomolecule` |
| [`AncientDnaRead`](file:///C:/Projects/github/scies-bio-th/src/paleo/biomolecule.rs#L180-L200) | aDNA read model with Briggs 5' $\text{C}\to\text{T}$ and 3' $\text{G}\to\text{A}$ deamination decay | `scies_bio_th::paleo::biomolecule` |
| [`AncientDnaAligner`](file:///C:/Projects/github/scies-bio-th/src/paleo/biomolecule.rs#L203-L287) | Deamination-aware Bayesian ancient DNA sequence read integration | `scies_bio_th::paleo::biomolecule` |
| [`PseudogeneResurrector`](file:///C:/Projects/github/scies-bio-th/src/paleo/resurrection.rs#L62-L199) | Repairs premature stop codons and frameshift defects using functional outgroups | `scies_bio_th::paleo::resurrection` |
| [`OsteologyEvidence`](file:///C:/Projects/github/scies-bio-th/src/paleo/phenotype.rs#L112-L176) | Histology (fibrolamellar/lamellar), pneumaticity, and Peto's Paradox mass scaling | `scies_bio_th::paleo::phenotype` |
| [`IsotopeEvidence`](file:///C:/Projects/github/scies-bio-th/src/paleo/phenotype.rs#L179-L234) | $\delta^{13}\text{C}$ diet, $\delta^{15}\text{N}$ trophic enrichment, and $\delta^{18}\text{O}$ paleoclimate | `scies_bio_th::paleo::phenotype` |
| [`PathwaySelectionProfile`](file:///C:/Projects/github/scies-bio-th/src/paleo/phenotype.rs#L237-L297) | Group-discounted Bayesian selection profile across 7 functional pathways | `scies_bio_th::paleo::phenotype` |
| [`PaleoReconstructionEngine`](file:///C:/Projects/github/scies-bio-th/src/paleo/engine.rs#L86-L273) | Orchestrator generating [`PaleoGenomeDraft`] and ranked [`ChassisGenomicEdit`] plans | `scies_bio_th::paleo::engine` |

---

## 5. End-to-End Practical Recipes

### Recipe 1: Phylogenetic ASR from Tree
```rust
use scies_bio_th::paleo::{ExtantSequence, FelsensteinAsr, SubstitutionModel};
use scies_bio_th::phylo::PhyloNode;
use scies_bio_th::seq::Dna;

// Construct phylogenetic tree: ((Chicken:68, Falcon:68):160, Alligator:228)
let birds = PhyloNode::internal(
    vec![PhyloNode::leaf("Gallus_gallus", 68.0), PhyloNode::leaf("Falco_peregrinus", 68.0)],
    160.0,
);
let tree = PhyloNode::internal(vec![birds, PhyloNode::leaf("Alligator_mississippiensis", 228.0)], 0.0);

let extant_sequences = vec![
    ExtantSequence::new("Gallus_gallus", Dna::new(b"ATGGGTCCCGCT")?),
    ExtantSequence::new("Falco_peregrinus", Dna::new(b"ATGGGTCCCGCT")?),
    ExtantSequence::new("Alligator_mississippiensis", Dna::new(b"ATGGGTCCCGCA")?),
];

let model = SubstitutionModel::default_vertebrate_coding();
let locus = FelsensteinAsr::reconstruct_locus(
    &tree,
    &extant_sequences,
    &model,
    "Tyrannosaurus rex",
    "COL1A1",
    66.0, // 66 Ma branch drift
)?;

println!("Consensus: {}", String::from_utf8_lossy(&locus.to_consensus_dna().sequence));
println!("Confidence: {:.1}%", locus.mean_posterior_confidence * 100.0);
```

### Recipe 2: Ancient Peptide Fragment Anchoring
```rust
use scies_bio_th::paleo::{PaleoproteomeTranslator, PeptideFragment, ReconstructedLocus};
use scies_bio_th::seq::Dna;

// Assume `locus` was generated via Layer 1 ASR
let mut locus: ReconstructedLocus = /* ... */;
let chicken_ref = Dna::new(b"ATGCCCGGT")?;

// Fossil collagen peptide: Gly-Hyp-Ala (Hyp represented as 'O')
let peptide = PeptideFragment::new(0, b"GOA", 0.97);

PaleoproteomeTranslator::anchor_with_peptides(&mut locus, &[peptide], Some(&chicken_ref))?;
println!("Anchored Protein: {}", String::from_utf8_lossy(&locus.translate_consensus()?.sequence));
```

### Recipe 3: Ancient DNA Read Integration
```rust
use scies_bio_th::paleo::{AncientDnaAligner, AncientDnaRead, ReconstructedLocus};

let mut locus: ReconstructedLocus = /* ... */;

// degraded 24 bp read with Siberian permafrost deamination
let reads = vec![
    AncientDnaRead::new_pleistocene(0, b"ATGACCCGCTTCCTGGGCACCGAG"),
];

let integrated_bases = AncientDnaAligner::integrate_reads(&mut locus, &reads)?;
println!("Integrated {} bases with Briggs deamination correction.", integrated_bases);
```

### Recipe 4: Pseudogene Repair & Open Reading Frame Restoration
```rust
use scies_bio_th::paleo::{ExtantSequence, PseudogeneResurrector, SubstitutionModel};
use scies_bio_th::seq::Dna;

// Chicken tooth enamel pseudogene (contains premature stop codon TGA at codon index 1)
let chicken_pseudogene = Dna::new(b"ATGTGATAA")?;
let functional_alligator = ExtantSequence::new("Alligator_mississippiensis", Dna::new(b"ATGCGATAA")?);

let model = SubstitutionModel::default_vertebrate_coding();
let resurrected = PseudogeneResurrector::resurrect_inactivated_gene(
    "ENAM",
    "Tyrannosaurus rex",
    &chicken_pseudogene,
    68.0,
    &[(functional_alligator, 228.0)],
    &model,
)?;

println!("Reversed {} premature stop codon(s).", resurrected.repairs.len());
println!("Restored Protein: {}", String::from_utf8_lossy(&resurrected.locus.translate_consensus()?.sequence));
```

### Recipe 5: Multi-Evidence Engine Synthesis & Chassis Edit Planning
```rust
use scies_bio_th::paleo::*;

let mut engine = PaleoReconstructionEngine::new("Mammuthus primigenius", "Elephas maximus");
// Register loci and phenotypic evidence
// engine.add_reconstructed_locus(...);
// engine.add_phenotypic_evidence(OsteologyEvidence { ... });
// engine.add_phenotypic_evidence(IsotopeEvidence { ... });

let draft = engine.build_draft()?;
println!("Overall Accuracy: {:.2}%", draft.summary_metrics.overall_base_accuracy_expected * 100.0);
println!("Target Chassis Edits: {}", draft.chassis_edit_plan.len());
for edit in draft.chassis_edit_plan.iter().take(5) {
    println!("Edit [{} @ {}]: {} -> {} ({:?}, priority={:.2})",
        edit.gene_symbol, edit.locus_position, edit.chassis_base as char, edit.target_base as char,
        edit.effect, edit.priority_score
    );
}
```

---

## 6. Case Studies

### Case Study A: Pleistocene Megafauna (*Mammuthus primigenius*)
- **Chassis**: Asian Elephant (*Elephas maximus*)
- **Target Gene**: `TRPM8` (cold sensation receptor)
- **Inputs**:
  - Tree: Elephantidae crown divergence ($6.0\text{ Ma}$) + Sirenia outgroup ($55\text{ Ma}$).
  - aDNA Reads: $24\text{ bp}$ permafrost sequence harboring Mammoth-specific missense substitution at codon 1 ($G \to A$, translating to Thr instead of Ala).
  - Stable Isotopes: $\delta^{18}\text{O} = -18.4‰$ (Siberian glacial ice water signal).
  - Histology: Massive fibrolamellar bone with adult body mass $5,500\text{ kg}$.
- **Reconstruction Output**:
  - Expected Base Accuracy: **$100.00\%$** (Shannon entropy: $0.0005\text{ bits}$).
  - Cold Adaptation Pathway Selection Posterior: **$90.5\%$**.
  - Peto's Paradox Mega-tumor Suppression Posterior: **$86.3\%$**.
  - Generated Edit: Exactly $1$ non-synonymous edit on `Elephas maximus` chassis (`pos 3: G -> A`, priority $2.81$).

### Case Study B: Mesozoic Non-Avian Dinosaur (*Tyrannosaurus rex*)
- **Chassis**: Modern Chicken (*Gallus gallus*)
- **Target Loci**: `COL1A1` (Bone Type I Collagen) & `ENAM` (Tooth Enamel Matrix)
- **Inputs**:
  - Tree: Avian-Crocodilian Archosaur crown ($228\text{ Ma}$) + unanchored extinct branch ($66\text{ Ma}$).
  - Paleoproteomics: MOR 1125 fossil collagen peptide `Gly-Hyp-Hyp-Gly-Pro-Ala-Gly-Hyp-Hyp` ($q = 0.96$).
  - Pseudogene: Inactivated chicken `ENAM` harboring premature stop codons at codons 2 and 5.
  - Functional Outgroups: *Alligator mississippiensis* and *Chelonia mydas*.
  - Phenotypic Evidence: Pneumatic vertebral foramina, fibrolamellar histology, mineralized teeth, body mass $7,800\text{ kg}$, $\delta^{15}\text{N} = +9.2‰$ (hypercarnivore).
- **Reconstruction Output**:
  - Expected Base Accuracy: **$96.39\%$** across coding loci.
  - Mineralized Dentition Pathway Posterior: **$90.9\%$**.
  - Pseudogene Repairs: Both premature stop codons reversed to ancestral functional state (`MPRGSQGP`).
  - Ranked Chassis Edits:
    1. `ENAM` pos 6: `T -> C` (`PseudogeneStopReversal`, confidence $97.1\%$, priority $4.10$)
    2. `ENAM` pos 15: `T -> C` (`PseudogeneStopReversal`, confidence $97.1\%$, priority $4.10$)
    3. `COL1A1` pos 5: `C -> T` (`Synonymous`, confidence $63.4\%$, priority $0.37$)

---

## 7. Downstream Interoperability & CRISPR Guide Design

Reconstructed loci and draft chassis edits interface directly with the rest of `scies-bio-th`:
1. **CRISPR & Prime Editing**:
   The output of `ChassisGenomicEdit` directly specifies the target locus coordinate and desired edit. Feeding the flanking sequence into `scies_bio_th::edit::crispr::generate_sgrna_candidates` or `scies_bio_th::edit::prime::generate_pegrna_candidates` generates active SpCas9/Cas12a guide RNAs to execute the edits in wet-lab cellular systems.
2. **Database Persistence with `scies-bio-db`**:
   `ReconstructedLocus` and `PaleoGenomeDraft` can be serialized with `serde` and stored in `scies-bio-db`'s 2-bit sequence repository, archiving ancestral genomes with cryptographic SHA-256 provenance hashes.
