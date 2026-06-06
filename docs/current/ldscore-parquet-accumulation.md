# LD-score Computation from R² Parquets — Mathematical Specification

This document defines, in math, **how LD scores are computed from a parquet R²
reference panel** by the streaming accumulator. It is the conceptual companion to
`docs/current/parquet-r2-format-and-read-pipeline.md` (the artifact format and the
read pipeline) and the design spec
`docs/superpowers/specs/2026-06-06-ldscore-parquet-pair-streaming-design.md`.

Every symbol is defined in §1 before it is used. Read §1 first.

---

## 1. Notation

Indices are **0-based** and refer to the *retained* SNP universe (the analysis
SNPs after intersection/restriction), in **genomic (position-sorted) order**.

| Symbol | Type | Meaning |
|---|---|---|
| $m$ | scalar | number of retained SNPs on the chromosome |
| $n_a$ | scalar | number of annotation columns (partitioned annotations **plus** the one regression-weight column) |
| $k,\ l$ | index in $\{0,\dots,m-1\}$ | generic retained-SNP indices (a "row" SNP $k$ and a "neighbour" SNP $l$) |
| $i,\ j$ | index in $\{0,\dots,m-1\}$ | the two endpoints of one *stored* pair, always with $i<j$ |
| $A$ | matrix $\mathbb{R}^{m\times n_a}$ | annotation matrix; $A_{l,a}$ = value of annotation $a$ at SNP $l$ |
| $A_{l,\cdot}$ | row vector $\mathbb{R}^{1\times n_a}$ | the whole annotation row for SNP $l$ |
| $r^2(k,l)$ | scalar in $(-\,\cdot,1]$ | **unbiased** squared correlation between SNPs $k$ and $l$ (may be slightly negative; see §6.4) |
| $\mathrm{win}(k)$ | set $\subseteq\{0,\dots,m-1\}$ | the LD-window neighbours of SNP $k$ (SNPs within the configured window of $k$) |
| $R$ | matrix $\mathbb{R}^{m\times m}$ | symmetric within-window LD matrix (defined in §3); **never materialized** |
| $C$ | matrix $\mathbb{R}^{m\times n_a}$ | the accumulator (`cor_sum` in code); on completion $C=RA$ |
| $\ell_{k,a}$ | scalar | the LD score of SNP $k$ for annotation $a$ (an entry of $C$) |
| $\mathrm{block\_left}[j]$ | index | smallest retained index inside the LD window of SNP $j$ (the window's left edge) |
| $P$ | set of triples | the stored, decoded, window-filtered pair set (defined in §4) |
| $v_{ij}$ | scalar | the final R² value of the stored pair $(i,j)$, $v_{ij}=r^2(i,j)$ (in vector form, $\mathbf{v}$ is a chunk's vector of these) |
| $\mathrm{remap}[\cdot]$ | array | maps a panel (build-order) index to its retained index, or $-1$ if the SNP is not retained |
| $\mathrm{IDX}_1,\mathrm{IDX}_2$ | on-disk int32 | a stored pair's two **panel** indices, with $\mathrm{IDX}_1<\mathrm{IDX}_2$ |
| $q$ | on-disk int16 | the quantized R² of a pair as stored on disk |
| $s$ | scalar | dequantization scale, $s=32767$ (`ldsc:r2_scale`) |
| $N$ | scalar | reference sample size (only used for the *raw*→unbiased correction of external panels; §6.4) |

Convention: $\mathrel{+}=$ means "add into" (accumulate); $\leftarrow$ means
assignment.

---

## 2. Purpose / what we compute

For each retained SNP $k$ and annotation $a$, the (partitioned) **LD score** is the
annotation-weighted sum of that SNP's within-window squared correlations,
**including the SNP itself** (self-correlation $1$):

$$
\ell_{k,a} \;=\; \underbrace{1\cdot A_{k,a}}_{\text{self }(l=k)} \;+\; \sum_{l\,\in\,\mathrm{win}(k)} r^2(k,l)\, A_{l,a}.
$$

The classic scalar LD score $\ell^2(k)=1+\sum_{l\in\mathrm{win}(k)} r^2(k,l)$ is the
special case $A_{\cdot,a}=\mathbf{1}$ (an all-ones annotation column).

---

## 3. Matrix form

Define the symmetric within-window LD matrix $R\in\mathbb{R}^{m\times m}$:

$$
R_{kl} \;=\;
\begin{cases}
1 & k=l \quad\text{(unit diagonal: self-correlation)}\\[2pt]
r^2(k,l) & l\in\mathrm{win}(k),\ k\neq l\\[2pt]
0 & \text{otherwise.}
\end{cases}
$$

Then the LD scores are exactly a matrix product:

$$
\boxed{\,C \;=\; R\,A\,}
\qquad\Longleftrightarrow\qquad
\ell_{k,a} \;=\; \sum_{l=0}^{m-1} R_{kl}\,A_{l,a}.
$$

$R$ is symmetric ($R_{kl}=R_{lk}$) because $r^2(k,l)=r^2(l,k)$, and it is **banded**:
$R_{kl}=0$ once $l$ leaves $k$'s window. $R$ is **never built**; it is represented
implicitly by its diagonal ($=1$) plus the stored off-diagonal pairs (§4), and the
product $RA$ is formed by streaming those pairs (§5).

---

## 4. What the parquet provides

The panel stores each **unordered** within-(build-)window pair **once**, in the
upper triangle. After the read-side per-pair transform (decode), a stored pair is a
triple $(i,j,v_{ij})$ where, for an on-disk row $(\mathrm{IDX}_1,\mathrm{IDX}_2,q)$:
$$
i=\mathrm{remap}[\mathrm{IDX}_1],\qquad
j=\mathrm{remap}[\mathrm{IDX}_2],\qquad
v_{ij}=\underbrace{q/s}_{\text{int16}\to\text{float}}\ \ (\text{unbiased already; §6.4}).
$$

Two filters define the consumed set $P$:

1. **Retention.** Drop pairs whose either endpoint is not in the analysis universe
   ($\mathrm{remap}=-1$).
2. **LD window.** The *ldscore* window (`--ld-wind-*`) may be **tighter** than the
   panel's build window, so a stored pair counts only if both endpoints are mutual
   window neighbours. For $i<j$ this is the single condition $i\ge\mathrm{block\_left}[j]$:

$$
P \;=\; \bigl\{\,(i,j,v_{ij})\ :\ i<j,\ \ i\ge\mathrm{block\_left}[j],\ \ \text{both retained}\,\bigr\}.
$$

Equivalence to $R$: for $i<j$,
$\ (i,j)\in\text{(unfiltered }P) \wedge i\ge\mathrm{block\_left}[j] \iff j\in\mathrm{win}(i) \iff R_{ij}=R_{ji}=v_{ij}.$
Diagonal entries $R_{kk}=1$ are **not** stored; they are added in §5 Step 1.

---

## 5. The streaming algorithm

Compute $C=RA$ by visiting every stored pair once. Because $R$ is symmetric, each
off-diagonal pair $(i,j)$ contributes to **two** rows of the output.

**Step 1 — seed the diagonal.** Since $R_{kk}=1$,

$$
C \;\leftarrow\; A
\qquad\bigl(\text{i.e. } C_{k,\cdot}\leftarrow A_{k,\cdot}\ \text{for every } k\bigr).
$$

**Step 2 — stream the off-diagonal pairs.** For each $(i,j,v_{ij})\in P$:

$$
C_{i,\cdot} \mathrel{+}= v_{ij}\,A_{j,\cdot},
\qquad\qquad
C_{j,\cdot} \mathrel{+}= v_{ij}\,A_{i,\cdot}.
$$

The first update is the term $R_{ij}A_{j,\cdot}$ in row $i$; the second is the
symmetric term $R_{ji}A_{i,\cdot}$ in row $j$.

**Step 3 — output.** After all pairs:

$$
C = RA,\qquad \ell_{k,a}=C_{k,a}.
$$

Columns $1..(n_a-1)$ are the reference/partitioned LD scores; the appended column is
the regression-universe weight $w_\ell$. Both come out of the **same pass**, because
the regression weight is just one more column of $A$.

### Vectorized form (per chunk)

Pairs arrive in chunks (one decoded row group). For a chunk with index vectors
$\mathbf{i},\mathbf{j}$ and value vector $\mathbf{v}$ (after applying the §4
filters), Step 2 is the unbuffered scatter-add

$$
C[\mathbf{i}] \mathrel{+}= \mathbf{v}\odot A[\mathbf{j}],
\qquad
C[\mathbf{j}] \mathrel{+}= \mathbf{v}\odot A[\mathbf{i}],
$$

where $A[\mathbf{j}]$ gathers the rows of $A$ at indices $\mathbf{j}$, and
$\mathbf{v}=(v_{ij})$ is the chunk's value vector whose entry scales each gathered
row (broadcast over the $n_a$ columns). **The scatter must be unbuffered**
(`np.add.at`): within a chunk
$\mathbf{i}$ contains repeats (one left SNP, many right neighbours), and a buffered
`+=` keeps only one of them. See §6.3.

---

## 6. Correctness notes

### 6.1 Symmetry is handled by the two updates
Each unordered pair is stored once but corresponds to two equal entries $R_{ij}=R_{ji}$.
Step 2's two scatter-adds place both. No pair is double-counted because $P$ holds
each unordered pair exactly once (parquet pair-uniqueness + injective remap; see the
format doc §2.2/§3.1).

### 6.2 The diagonal is counted exactly once
Every SNP's self term $R_{kk}A_{k,\cdot}=A_{k,\cdot}$ is added exactly once by the
Step 1 seed $C\leftarrow A$, and never by Step 2 (stored pairs have $i<j$, so they
are strictly off-diagonal).

### 6.3 Duplicate left indices
Within one chunk the left index $i$ repeats across a SNP's many neighbours. The
accumulation is therefore a true scatter-**add**, not an assignment: $C_i$ receives
the sum of all $v_{ij}\,A_{j,\cdot}$ for that $i$. A buffered $C[\mathbf{i}]\mathrel{+}=\dots$
would drop all but one contribution per repeated index, undercounting LD scores.

### 6.4 The stored R² is already unbiased
`build-ref-panel` writes the **unbiased** estimate $r^2$ (`_unbiased_r2_array`,
tagged `ldsc:r2_bias=unbiased`), so the decoded $v_{ij}=q/s$ is the final value and
no correction is applied. The raw→unbiased map
$v_{ij}\leftarrow v_{ij}-\tfrac{1-v_{ij}}{N-2}$ is applied **only** for external
*raw*-R² panels (`--r2-bias-mode raw --r2-sample-size N`), once per pair at decode.
Either way the $v_{ij}$ entering §5 is the final R². Possible small negative
$v_{ij}$ (≈15% of pairs) are preserved, not floored.

### 6.5 Numerical precision
$C$ accumulates in `float64` from `float32` products. The summation order differs
from the legacy genotype-block path, so results are **not bit-identical** (≈1e-6
drift) — well within LD-score tolerances (int16 quantization alone already perturbs
per-SNP LD scores by ≤~2e-3).

---

## 7. Complexity

Let $\mathrm{nnz}=|P|$ be the number of consumed pairs.

$$
\text{time} = O(\mathrm{nnz}\cdot n_a)\quad(\text{each pair touched once per annotation}),
\qquad
\text{memory} = \underbrace{O(m\,n_a)}_{C} + \underbrace{O(\text{one row group})}_{\text{stream}}.
$$

This is the optimal arithmetic for $C=RA$ given $R$ as an edge list, and the peak
memory is independent of the LD window width (the window only filters $P$). There is
no block tiling, no dense $R$ tile, and no decoded-row-group cache.
