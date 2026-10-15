# LD-score Computation from R² Parquets — Mathematical Specification

Last updated on: 2026-10-15

This document defines, in math, **how LD scores are computed from a parquet R²
reference panel** by the streaming accumulator. It is the conceptual companion to
[the artifact format and read pipeline](parquet-r2-format-and-read-pipeline.md).
The [original pair-streaming specification](../specs/2026-06-06-ldscore-parquet-pair-streaming-design.md) records the earlier full-reference CSR design; the current selected-output realization is described below.

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
| $\mathrm{win}(k)$ | set $\subseteq\{0,\dots,m-1\}$ | the other SNPs within the configured LD window of $k$, excluding $k$ itself |
| $R$ | matrix $\mathbb{R}^{m\times m}$ | symmetric within-window LD matrix (defined in §3); **never materialized** as a dense matrix |
| $U$ | sparse matrix $\mathbb{R}^{m\times m}$ | strict upper triangle of the stored pairs: $U_{ij}=v_{ij}$ for $(i,j)\in P$, $0$ elsewhere; gives $R=I+U+U^{\top}$ (§5) |
| $C$ | matrix $\mathbb{R}^{m\times n_a}$ | conceptual full-reference result $C=RA$; production retains only selected output rows |
| $\ell_{k,a}$ | scalar | the LD score of SNP $k$ for annotation $a$ (an entry of $C$) |
| $\mathrm{block\_left}[j]$ | index | smallest retained index inside the LD window of SNP $j$ (the window's left edge) |
| $P$ | set of triples | the stored, decoded, window-filtered pair set (defined in §4) |
| $v_{ij}$ | scalar | the final R² value of the stored pair $(i,j)$, $v_{ij}=r^2(i,j)$ (in vector form, $\mathbf{v}$ is a chunk's vector of these) |
| $\mathrm{remap}[\cdot]$ | array | maps a panel (build-order) index to its retained index, or $-1$ if the SNP is not retained |
| $\mathrm{IDX}_1,\mathrm{IDX}_2$ | on-disk int32 | a stored pair's two **panel** indices, with $\mathrm{IDX}_1<\mathrm{IDX}_2$ |
| $q$ | on-disk int16 | the quantized R² of a pair as stored on disk |
| $s$ | scalar | dequantization scale, $s=32767$ (`ldsc:r2_scale`) |
| $N$ | scalar | reference sample size (only used for the *raw*→unbiased correction of external panels; §6.4) |
| $K$ | scalar | maximum decoded pairs buffered before projection (§5) |
| $O$ | ordered retained-row indices | output SNP rows, a subset of the contributor universe |
| $C_{\mathrm{out}}$ | matrix $\mathbb{R}^{\lvert O\rvert\times n_a}$ | production float64 accumulator, equal to $(RA)[O,:]$ |

Convention: $\mathrel{+}=$ means "add into" (accumulate); $\leftarrow$ means
assignment; $\mathbf{v}\odot M$ is the row-wise scaling of matrix $M$ by vector
$\mathbf{v}$ (Hadamard product broadcast over columns, i.e.
$(\mathbf{v}\odot M)_{t,a}=\mathbf{v}_t\,M_{t,a}$, equivalently
$\operatorname{diag}(\mathbf{v})\,M$).

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

### Sparse realization (chunked CSR)

The conceptual full product can be expressed as two sparse matrix–dense products (SpMM) instead of element-wise scatter-adds.
With $R=I+U+U^{\top}$ (§1),

$$
C \;=\; RA \;=\; \underbrace{A}_{\text{diagonal}} \;+\; \underbrace{U A}_{i\text{-side}} \;+\; \underbrace{U^{\top}A}_{j\text{-side}}.
$$

$U$ is stored in **compressed sparse row (CSR)** form: a `data` array of nonzero
values (row by row), an `indices` array of their column indices, and an `indptr`
array of length $m{+}1$ where row $k$'s entries are `data[indptr[k]:indptr[k+1]]`.

**Numeric example** (a panel with $m=4$, pairs $(0,1){=}0.4,\ (0,2){=}0.2,\
(1,2){=}0.6,\ (2,3){=}0.5$, single all-ones annotation $A=[1,1,1,1]^{\top}$):

$$
U=\begin{bmatrix}0&0.4&0.2&0\\0&0&0.6&0\\0&0&0&0.5\\0&0&0&0\end{bmatrix},
\qquad
\begin{aligned}
\texttt{data}&=[\,0.4,\ 0.2,\ 0.6,\ 0.5\,]\\
\texttt{indices}&=[\,1,\ \ \ 2,\ \ \ 2,\ \ \ 3\,]\\
\texttt{indptr}&=[\,0,\ \ \ 2,\ \ \ 3,\ \ \ 4,\ \ \ 4\,]
\end{aligned}
$$

$$
UA=\begin{bmatrix}0.6\\0.6\\0.5\\0\end{bmatrix}\!\text{(row sums)},\quad
U^{\top}A=\begin{bmatrix}0\\0.4\\0.8\\0.5\end{bmatrix}\!\text{(col sums)},\quad
C=A+UA+U^{\top}A=\begin{bmatrix}1.6\\2.0\\2.3\\1.5\end{bmatrix}.
$$

This upper-triangle representation explains the two directed contributions.
Production combines the required orientations into a rectangular selected-output operator, described below; it does not allocate the full $U$ for this example as a separate workflow step.

#### Chunking

Materializing $U$ for a whole chromosome would hold all $|P|$ pairs in memory at
once. Instead the stored pairs are partitioned into chunks of $K$ pairs
($U=\sum_c U_c$), and each chunk's two SpMMs are accumulated:

$$
C \;\leftarrow\; A \;+\; \sum_c\bigl(U_c A + U_c^{\top}A\bigr).
$$

This identity holds in exact arithmetic by linearity of the matrix product,
$\bigl(\sum_c U_c\bigr)A=\sum_c U_c A$: no pair is dropped or double-counted because
each $(i,j)$ lands in exactly one chunk. Splitting the example after the first two
pairs:

$$
U_1=\begin{bmatrix}0&0.4&0.2&0\\0&0&0&0\\0&0&0&0\\0&0&0&0\end{bmatrix},\quad
U_2=\begin{bmatrix}0&0&0&0\\0&0&0.6&0\\0&0&0&0.5\\0&0&0&0\end{bmatrix},\quad
U_1+U_2=U,
$$

$$
U_1 A=\begin{bmatrix}0.6\\0\\0\\0\end{bmatrix},\quad
U_2 A=\begin{bmatrix}0\\0.6\\0.5\\0\end{bmatrix},\quad
U_1 A+U_2 A=\begin{bmatrix}0.6\\0.6\\0.5\\0\end{bmatrix}=UA.
$$

#### Current selected-output implementation

The workflow needs $C_{\mathrm{out}}=(RA)[O,:]$, not every row of $C$.
Each output still sums over all retained reference contributors.
For chunk $c$, let $W_c=(U_c+U_c^\top)[O,:]$; then

$$
C_{\mathrm{out}} = A[O,:] + \sum_c W_c A.
$$

[`ProjectionAccumulator.add_pairs()`](../../src/ldsc/_kernel/ldscore_projection.py) builds the needed directed edges only for output targets, compacts the involved target/contributor indices, and uses a rectangular float64 CSR for $W_c$.
If $O=[0,2]$ in the example, the selected scores are $[1.6,2.3]^\top$.
SNPs 1 and 3 still contribute; dropping them before multiplication would give the wrong answer.
The diagonal is seeded once on output rows by `add_diagonal()`.

[`ld_score_streaming_from_r2_reader()`](../../src/ldsc/_kernel/ldscore.py) splits decoded record batches so at most $K=16{,}000{,}000$ pairs are buffered before projection, with a final partial chunk.
$K$ is independent of Parquet row-group size and PLINK `snp_batch_size`; it is not clamped below by $m$.
The pair buffer, concatenated arrays, two-direction expansion, CSR and sorting workspace can coexist, so $K$ is not a byte/RSS cap.
The reader uses record batches of at most 65,536 rows even for oversized row groups.

`ProjectionAccumulator._accumulate()` tiles annotation columns and contributor/output rows.
It reads a float64 contributor tile once and reuses it across output blocks.
The default 16 MiB read budget targets one operand tile, with at least one row, rather than the complete working set.
Sparse slices, gathered accumulator entries, products, metadata, and conversion buffers have separate lifetimes.
The main score array spans output rows and the workflow's active annotation batch, including shared baseline and weight columns.
Each later execution batch traverses the reference again.

---

## 6. Correctness notes

### 6.1 Symmetry is handled by the two updates
Each unordered pair is stored once but corresponds to two equal entries $R_{ij}=R_{ji}$.
Step 2's two updates place both contributions in the conceptual full result.
Production retains the orientations whose target is in $O$.
No pair is double-counted because $P$ holds
each unordered pair exactly once (parquet pair-uniqueness + injective remap; see the
format doc §2.2/§3.1).

### 6.2 The diagonal is counted exactly once
Every SNP's self term $R_{kk}A_{k,\cdot}=A_{k,\cdot}$ is added exactly once by the
Step 1 seed $C\leftarrow A$, and never by Step 2 (stored pairs have $i<j$, so they
are strictly off-diagonal).

### 6.3 Repeated indices within a chunk
A SNP's left index $i$ repeats across its many neighbours, so several stored entries
share the same CSR row. The SpMM sums them by construction — row $i$ of $UA$ is the
sum over that row's entries, $\sum_{(i,j)\in P_c} v_{ij}A_{j,\cdot}$ — so repeated
indices accumulate correctly with no element-wise scatter and no risk of a buffered
overwrite. (This is what an earlier `np.add.at` scatter did element by element and
far more slowly.)

### 6.4 The stored R² is already unbiased
`build-r2-panel` writes the **unbiased** estimate $r^2$ (`_unbiased_r2_array`,
tagged `ldsc:r2_bias=unbiased`), so the decoded $v_{ij}=q/s$ is the final value and
no correction is applied. The raw→unbiased map
$v_{ij}\leftarrow v_{ij}-\tfrac{1-v_{ij}}{N-2}$ is applied **only** for external
*raw*-R² panels (those declaring `ldsc:r2_bias=raw` with `ldsc:n_samples=N` in
parquet schema metadata), once per pair at decode.
Either way the $v_{ij}$ entering §5 is the final R². Negative decoded $v_{ij}$ values are preserved, not floored.

### 6.5 Numerical precision
The score accumulator, sparse operator values, and annotation operands are float64.
This avoids accumulating many contributions in float32; public LD-score values are narrowed only at the established output boundary.
Annotation bit-packing is lossless, while source continuous annotations and decoded R² retain their existing float32 precision before promotion.
Int16 panel quantization is a separate lossy operation; promoting decoded values cannot restore the lost information.

Linearity proves the mathematical result, not bitwise equality between floating-point evaluation orders.
Changing $K$, tile layout, query width, or backend can reorder sums and alter roundoff.
Validate against independent float64 expectations with the established numerical tolerances, and compare transport-only changes at matched layouts.
There is no universal per-SNP LD-score error bound independent of pair count and annotation weights.
See [`test_ldscore_projection.py`](../../tests/test_ldscore_projection.py), [`test_ldscore_streaming.py`](../../tests/test_ldscore_streaming.py), and the [transport parity audit](../audits/2026-10-15_direct-streaming.md).

---

## 7. Complexity

Let $\mathrm{nnz}=|P|$ be the number of consumed pairs and $n_a$ the active annotation width including the optional weight.
The leading multiplication work is $O(\mathrm{nnz}\,n_a)$, with fewer directed updates when only some targets are output rows.
Reading, remapping, filtering, sorting/CSR construction, compression, and repeated execution batches add costs outside that arithmetic expression.

The main score payload is $8\,|O|\,n_a$ bytes.
Pair/operator workspace scales with the chunk bound $K$, and sidecars/maps with retained SNP count; annotation and output tiles introduce additional allocations.
There is no dense chromosome-wide $R$ or decoded-row-group cache.
This bounds particular components independently of total stored pair count, but it is not a formula for total process or job RSS.
The read-time window can change the retained sparse operator even though every stored pair is scanned.
Worker concurrency, genotype/index state, filesystem cache, and staged final outputs must be accounted for separately.
See [memory and CPU efficiency design](annotation-memory-design.md#precision-universes-and-resource-accounting).
