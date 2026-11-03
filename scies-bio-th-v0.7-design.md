# scies-bio-th — v0.7 Design Document
**Computational Precision & Structural Biology**

---

## 1. Context and Architectural Stance

`scies-bio-th` spans from fundamental 1D sequence analysis and precision CRISPR engineering (v0.1–v0.6) to **3D macromolecular structures, biological network systems, and statistical rigor (v0.7)**.

While maintaining our unwavering principle of **Zero External Dependencies** in the default build (pure standard library `std`), v0.7 establishes modular feature boundaries:
- `structure`: 3D coordinates, PDB/mmCIF streaming parsing, DSSP-like secondary structure, Kabsch RMSD & TM-score, and contact maps.
- `network`: Biological graphs, Protein-Protein Interaction (PPI) networks, KEGG-style pathway DAGs, and topological centrality.
- `stats`: Non-parametric hypothesis testing (Mann-Whitney U, Fisher's Exact, Chi-square), multiple testing corrections (Bonferroni, Benjamini-Hochberg FDR), and vector distance metrics.
- `docking`: Grid-based binding pocket detection and physics-based scoring (Van der Waals Lennard-Jones 6-12 + Coulombic electrostatics).
- `epigenome`: CpG island detection (Gardiner-Garden & Frommer), beta/M-value methylation transformations, ATAC-seq chromatin accessibility, and histone modification states.
- `experimental`: Isolates speculative and synthetic biology pipelines (e.g. `ivg`) from the core computation engine.

---

## 2. Module Tree

```text
scies_bio_th::
├── (v0.1–v0.6) seq, io, align, codon, edit, genome, annotation, variant,
│               motif, hmm, phylo, popgen, genomics, transcriptomics,
│               proteomics, metagenomics, aging
│
├── structure          — Protein 3D structure primitives
│   ├── pdb            — PDB and mmCIF streaming coordinate parser
│   ├── secondary      — DSSP-like secondary structure assignment (H/E/C)
│   ├── rmsd           — Kabsch optimal superposition RMSD and TM-score
│   └── contact        — Cα-Cα distance matrix and contact maps
│
├── network            — Biological network analysis
│   ├── graph          — Generic directed/undirected biological graph
│   ├── ppi            — Protein-Protein Interaction networks and hubs
│   ├── pathway        — KEGG-style pathway DAG and topological sorting
│   └── centrality     — Degree, Betweenness (Brandes), and Closeness centrality
│
├── stats              — Statistical methods for biological data
│   ├── test           — Mann-Whitney U, Fisher's Exact (2x2), Chi-square
│   ├── correction     — Bonferroni and Benjamini-Hochberg (BH) FDR
│   └── dist           — Euclidean, Manhattan, and Cosine distance
│
├── docking            — Molecular docking scoring (heuristic)
│   ├── grid           — 3D grid cavity and binding pocket detection
│   └── score          — Van der Waals (LJ 6-12) + distance-dependent dielectric electrostatics
│
├── epigenome          — Epigenomics support
│   ├── methylation    — CpG island detection and Beta/M-value conversion
│   ├── chromatin      — ATAC-seq chromatin accessibility index
│   └── histone        — Histone modification states and chromatin segmentation
│
└── [experimental] ivg — In vitro gametogenesis (moved behind feature gate)
```

---

## 3. Mathematical Specifications

### 3.1 Kabsch Optimal Superposition RMSD and TM-Score
Given two coordinate sets $P, Q \in \mathbb{R}^{N \times 3}$:
1. Center both around their respective centroids: $P_c = P - \bar{P}, \ Q_c = Q - \bar{Q}$.
2. Compute cross-covariance matrix $H = P_c^T Q_c$.
3. Compute rotation matrix $R$ minimizing $\| P_c R - Q_c \|_F$.
4. Root Mean Square Deviation:
   $$\text{RMSD} = \sqrt{\frac{1}{N} \sum_{i=1}^N \| P_{c,i} R - Q_{c,i} \|^2}$$
5. Template Modeling Score (TM-score):
   $$\text{TM-score} = \frac{1}{L_{\text{target}}} \sum_{i=1}^{L_{\text{aligned}}} \frac{1}{1 + \left(\frac{d_i}{d_0(L_{\text{target}})}\right)^2}, \quad d_0(L) = 1.24 \sqrt[3]{L - 15} - 1.8$$

### 3.2 Doench 2016 CFD Matrix for CRISPR Off-Target Prediction
Calculates cutting frequency determination using position-dependent mismatch weight matrix $W(\text{pos}, r_i, d_i)$:
$$\text{CFD} = \prod_{i=1}^{20} W(i, \text{sgRNA}_i, \text{DNA}_i) \times 100.0$$

### 3.3 Banded Needleman-Wunsch Alignment
For sequence lengths $M, N$ with band radius $k$:
Only dynamic programming cells $(i, j)$ where $|i - j| \le k$ are evaluated, reducing complexity from $O(M \cdot N)$ to $O(k \cdot \min(M, N))$, essential for chromosomes $\ge 10\text{k bp}$.

### 3.4 CpG Island Gardiner-Garden & Frommer Criteria
A genomic window $[i, i + W]$ is classified as a CpG island if:
1. Window length $W \ge 200\text{ bp}$.
2. $\text{GC content} = \frac{G + C}{W} \ge 0.50$.
3. $\text{Obs/Exp CpG} = \frac{\text{Count}(\text{CpG}) \cdot W}{(C \cdot G)} \ge 0.60$.
