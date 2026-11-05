# scies-bio-th — Paleo Genomic & Proteomic Reconstruction Design
**v1.0-paleo · Mechanistic Ancestral Reconstruction, Paleoproteomic Back-Translation & Phenotypic Selection**

**Status**: Implemented (v1.0)  
**Feature Flag**: `paleo`  
**Core Philosophy**: Pure Rust, zero external dependencies, commutative Bayesian inference, sequence-level output (`Dna`, `Protein`, and `ChassisGenomicEdit` diffs)

---

## 1. Scientific Foundation & Realistic Recovery Tiers

Biomolecular preservation follows strict thermodynamic and kinetic decay laws:

1. **Tier 1 — Late Pleistocene to ~1.5 Ma (Ancient DNA / aDNA)**:
   - Short DNA fragments (30–100 bp) survive in permafrost, tooth cementum, or petrous bone.
   - Primary artifact: **Post-mortem cytosine deamination** ($\text{C} \to \text{U}$, read as $\text{T}$ at 5' overhangs; $\text{G} \to \text{A}$ at 3' overhangs).
   - **Recovery**: >99% nucleotide accuracy when anchored against a close extant relative (e.g., *Mammuthus primigenius* against *Elephas maximus*).

2. **Tier 2 — Early Pleistocene / Miocene 1.5–10+ Ma (Paleoproteomics)**:
   - DNA is completely degraded, but biomineralized proteins (**tooth enamelin `ENAM`, amelogenin `AMELX`, and bone Type I collagen `COL1A1`/`COL1A2`**) survive for millions of years.
   - **Recovery**: Mass spectrometry yields peptide sequences. Combining deterministic 1st/2nd codon positions with phylogenetic codon-usage priors from extant relatives recovers coding DNA sequences with >90% accuracy.

3. **Tier 3 — Deep Time / Mesozoic 66+ Ma (Extant Phylogenetic Bracketing + Pseudogene Resurrection)**:
   - Relies on the **Extant Phylogenetic Bracket** (e.g., Birds + Crocodilians + Turtles for non-avian dinosaurs) via Continuous-Time Markov Chain (CTMC) **Ancestral Sequence Reconstruction (ASR)**.
   - **Pseudogene Resurrection**: Lineages that lost ancestral traits (e.g., birds losing teeth) retain inactivated pseudogene remnants (`ENAM`, `MMP20`) in their genomes. By aligning avian pseudogenes against functional crocodilian outgroups, frameshifts and premature stop codons (`TGA`/`TAA`/`TAG`) are repaired to their ancestral functional state.
   - **Phenotypic Priors**: Osteology, histology, body mass (Peto's Paradox), and stable isotopes ($\delta^{13}\text{C}, \delta^{15}\text{N}, \delta^{18}\text{O}$) modulate pathway-level selection ($\omega = dN/dS$) and prioritize functional edits on the extant chassis genome.

---

## 2. Tri-Layer Architecture

```text
scies_bio_th::paleo
├── mod.rs                // Public facade & re-exports
├── types.rs              // BasePosterior, ReconstructedLocus, LocusProvenance, PaleoError
├── asr.rs                // Layer 1: SubstitutionModel (JC69, K80), Felsenstein pruning ASR
├── biomolecule.rs        // Layer 2: PaleoproteomeTranslator & AncientDnaAligner (Briggs model)
├── resurrection.rs       // Layer 2.5: PseudogeneResurrector (premature stop & frameshift repair)
├── phenotype.rs          // Layer 3: Commutative group-discounted Bayesian pathway selection
└── engine.rs             // PaleoReconstructionEngine, ChassisGenomicEdit, PaleoGenomeDraft
```

---

## 3. Mathematical Specifications

### 3.1 Per-Site Nucleotide Posterior (`BasePosterior`)
Every nucleotide site $i$ maintains a normalized probability vector over $\{A, C, G, T\}$:
$$\mathbf{p}_i = (p_A, p_C, p_G, p_T), \quad \sum_{b \in \{A,C,G,T\}} p_b = 1, \quad H(\mathbf{p}_i) = -\sum_{b} p_b \log_2(p_b) \in [0, 2]\text{ bits}$$

### 3.2 Layer 1: Felsenstein Pruning & Extinct Branch Propagation (`asr`)
Given a rooted phylogenetic tree $\mathcal{T}$ of extant relatives and a substitution rate matrix $Q$ (JC69 or Kimura 2-Parameter K80 with transition/transversion ratio $\kappa$):
1. Upward post-order traversal computes conditional subtree likelihoods $L_u(b)$ at every internal node $u$:
   $$L_u(b) = \prod_{v \in \text{children}(u)} \left( \sum_{b' \in \{A,C,G,T\}} P_{b \to b'}(t_{u,v}) L_v(b') \right)$$
2. At the common ancestor node $r$, marginal posterior probabilities are:
   $$P_{\text{anc}}(b) = \frac{\pi_b L_r(b)}{\sum_{b'} \pi_{b'} L_r(b')}$$
3. Forward propagation along the extinct lineage branch of length $\Delta t_{\text{extinct}}$:
   $$P_{\text{extinct}}(b') = \sum_{b \in \{A,C,G,T\}} P_{\text{anc}}(b) P_{b \to b'}(\Delta t_{\text{extinct}})$$

### 3.3 Layer 2A: Bayesian Paleoproteomic Back-Translation (`biomolecule`)
Given an observed fossil amino acid $aa$ at codon position $k$ with mass-spectrometry confidence $q \in (0, 1]$:
- Prior over the 64 codons $c = (b_1, b_2, b_3)$ is formed from Layer 1 posteriors weighted by extant codon usage bias $u(c)$:
  $$P_0(c) \propto P_1(b_1) P_2(b_2) P_3(b_3) \cdot u(c)$$
- Likelihood of observing $aa$:
  $$L(aa \mid c) = \begin{cases} q + \frac{1-q}{20} & \text{if } \text{translate}(c) = aa \\ \frac{1-q}{20} & \text{otherwise} \end{cases}$$
- Marginalizing the posterior $P(c \mid aa) \propto P_0(c) L(aa \mid c)$ over each codon position yields the updated per-base posteriors $P'_1(b_1), P'_2(b_2), P'_3(b_3)$.

### 3.4 Layer 2B: Briggs Deamination-Aware Ancient DNA Model (`biomolecule`)
For an ancient DNA read aligned at position $i$ at distance $d_{5'}$ from the 5' end and $d_{3'}$ from the 3' end:
- Post-mortem deamination probabilities decay exponentially from read termini:
  $$\delta_{5'}(d_{5'}) = \delta_{0,5'} e^{-\lambda d_{5'}}, \quad \delta_{3'}(d_{3'}) = \delta_{0,3'} e^{-\lambda d_{3'}}$$
- If the true ancestral base is $\text{C}$, the probability of observing $\text{T}$ in the read is $(1-\epsilon)\delta_{5'} + \frac{\epsilon}{3}$, preventing false-positive $\text{C} \to \text{T}$ mutations at fragment ends.

### 3.5 Layer 3: Commutative Group-Discounted Phenotypic Pooling (`phenotype`)
To guarantee order-invariance (commutativity) and prevent overconfidence from correlated fossil measurements, evidences are partitioned into `CorrelationGroup` buckets $g \in \mathcal{G}$. With group discount $\alpha_g = 1/\sqrt{|g|}$:
$$\text{logit}\, P(M \mid \mathcal{E}) = \text{logit}\, P_0(M) + \sum_{g \in \mathcal{G}} \frac{1}{\sqrt{|g|}} \sum_{e \in g} r_e \cdot \ell_e(M)$$
This guarantees:
- **Order invariance**: Adding evidence in any order produces the exact same posterior.
- **No hard-cut fragility**: Conflicting evidence shifts posterior log-odds smoothly without zeroing out valid hypotheses.
- **Correlation damping**: Five osteological measurements on the same femur are discounted by $1/\sqrt{5}$.
