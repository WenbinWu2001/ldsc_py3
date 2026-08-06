# Exact gene LD-score index: mathematical algorithm

Last updated on: 2026-08-06

This document gives the input-to-output mathematical specification for the
exact gene-list index used by `ldsc ldscore`. It describes the online indexed
algorithm in full and the offline `build-gene-ldscore-index` construction more
concisely. The index is an exact factorization of the matching direct PLINK
calculation, not an approximation and not a sum of separately calculated
per-gene LD scores.

## 1. Scope and chromosome-local notation

All definitions below are chromosome-local unless a sum over chromosomes is
shown. Rows use the PLINK-authoritative genomic order after the baseline/PLINK
identifier-key inner intersection and genotype filtering.

| Symbol | Shape | Meaning |
| --- | --- | --- |
| $n$ | scalar | selected PLINK reference individuals |
| $m$ | scalar | retained LD-reference SNPs |
| $r$ | scalar | persisted regression/output SNPs, with $r\le m$ |
| $b$ | scalar | supplied baseline annotation columns |
| $g$ | scalar | embedded catalog genes |
| $t$ | scalar | disjoint gene atoms |
| $X$ | $n\times m$ | centered and variance-normalized PLINK genotypes |
| $A$ | $m\times b$ | supplied baseline annotations on retained reference SNPs |
| $R$ | $m\times m$ | within-window adjusted-$r^2$ LD operator |
| $p$ | $m\times 1$ | Boolean persisted-row/regression-SNP selector |
| $P$ | $r\times m$ | row-selection matrix induced by $p$ |
| $G$ | $g\times t$ | Boolean gene-to-atom membership |
| $H$ | $m\times t$ | Boolean reference-SNP-to-atom membership |
| $a_q$ | $g\times 1$ | Boolean catalog-gene selector for focal query $q$ |
| $z_q$ | $t\times 1$ | Boolean atom selector for focal query $q$ |
| $c$ | $m\times 1$ | Boolean common-SNP mask |
| $Y$ | $r\times t$ | stored gene LD-score operator $PRH$ |

The four SNP universes must not be conflated:

1. The **source universe** is the raw baseline and PLINK input.
2. The **LD-reference universe** of size $m$ is their identifier-key inner
   intersection after genotype usability checks, optional individual selection,
   and optional retained-reference MAF filtering.
3. The **common universe** is the subset selected by
   $c_i=\mathbf{1}[\operatorname{MAF}_i\ge\tau]$, where
   $\tau=\texttt{common_maf_min}$ (default 0.05).
4. The **persisted regression/output universe** of size $r$ is selected from the
   reference universe by the bundled HapMap3 or custom regression-SNP keys and
   then by the configured regression-region subtraction.

LD-score contributions, annotation counts, and annotation overlaps use the
LD-reference universe. Only persisted rows are written to the output tables.
Consequently, a nonpersisted SNP—including an MHC or centromeric SNP removed by
the regression-region policy—can contribute LD, counts, and overlaps.

## 2. The direct quantity being factorized

For retained SNPs $i$ and $j$, let

$$
\rho_{ij}=\frac{1}{n}X_{\cdot i}^{\mathsf T}X_{\cdot j}.
$$

Let $\nu=n-2$ when $n>2$ and $\nu=n$ otherwise, matching the PLINK kernel. The
adjusted squared-correlation contribution is

$$
\widetilde r_{ij}^{\,2}
=\rho_{ij}^{2}-\frac{1-\rho_{ij}^{2}}{\nu}.
$$

Define the windowed LD operator

$$
R_{ij}=
\begin{cases}
\widetilde r_{ij}^{\,2}, & i\text{ and }j\text{ are within the configured cM window},\\
0, & \text{otherwise}.
\end{cases}
$$

The standardized self-correlation makes $R_{ii}=1$. Off-diagonal adjusted
values may be negative; neither construction nor indexed assembly clamps them.

For any reference-SNP annotation vector $q\in\mathbb{R}^{m}$, its LD score on
the persisted rows is

$$
\ell(q)=PRq.
$$

The builder evaluates the same PLINK block kernel used by direct `ldscore`, so
$R$ is conceptual and is not materialized as a dense chromosome-wide matrix.

## 3. Concise offline index construction

### 3.1 Align the scientific inputs

For each chromosome, the builder:

1. reads and validates the baseline annotation shards;
2. prepares the selected PLINK individuals and usable genotype rows;
3. drops every duplicate effective-identity group independently in the mutable
   baseline and PLINK sources, then inner-joins by `SNP` for `rsid` or normalized
   `(CHR, POS)` for `chr_pos`;
4. orders the result canonically by PLINK-authoritative `CHR`, `POS`, and `SNP`
   while retaining the original BIM/BED column permutation for genotype reads; and
5. forms $A$, $X$, the common mask $c$, and the persisted selector $p$.

Baseline-only and PLINK-only SNPs are dropped. A duplicate group is removed in
full, warned, and recorded; no biologically ambiguous representative is chosen.
PLINK supplies the persisted `CHR`, `POS`, `SNP`, `A1`, and `A2`, while `SNP`
is a passive label rather than identity in `chr_pos` mode.

### 3.2 Precompute the fixed baseline and weight columns

The index stores the fixed baseline LD-score block

$$
L_A=PRA\in\mathbb{R}^{r\times b}
$$

and the regression-weight LD score

$$
w=PRp\in\mathbb{R}^{r}.
$$

They are evaluated together. The builder appends $p$ after the baseline
columns and makes one PLINK-kernel call,

$$
R[A\;p]=[RA\;Rp],
$$

then selects the persisted rows and splits the result into $L_A$ and $w$.
Thus baseline and regression-weight construction share one genotype-correlation
traversal without changing either numerical definition.

These columns, together with the persisted SNP identity rows, become
`baseline_rows.parquet`. Online gene-list runs reuse them without reading the
source baseline, PLINK, or regression-SNP inputs.

### 3.3 Convert padded genes into disjoint atoms

After optional gene-region filtering, gene $u$ with zero-based half-open
transcribed interval $[s_u,e_u)$ and padding $h$ has the projected interval

$$
I_u=[\max(0,s_u-h),\ e_u+h).
$$

The sorted union of all interval boundaries partitions the covered genome into
maximal nonempty half-open atoms. Within each atom, the set of covering genes is
constant. The builder records this relation as the Boolean matrix $G$.

Each retained reference SNP belongs to at most one atom. Its membership is
recorded in $H$, so every row of $H$ has either one nonzero or is all zero.

### 3.4 Store the exact LD-score operator

The central indexed operator is

$$
\boxed{Y=PRH}.
$$

Construction evaluates $RH$ in bounded atom-column blocks, selects the
persisted rows with $P$, and stores the result as a float64 CSR matrix in
`ldscore_operator.npz`. Atom blocking changes memory use only; linearity makes
it numerically the same operation as evaluating all columns together under the
same PLINK kernel.

These repeated atom-block kernel calls are distinct from the fixed
$R[A\;p]$ pass above and are intentional. They bound the dense working matrix
by the configured atom batch size; the builder does not fuse or materialize all
columns of $H$ in the baseline/weight pass.

### 3.5 Store count and overlap sufficient statistics

Let $\mathbf{1}$ be an $m$-vector of ones and let
$C=\operatorname{diag}(c)$. The builder stores

$$
\begin{aligned}
d      &=H^{\mathsf T}\mathbf{1},
&d_c   &=H^{\mathsf T}c,\\
D      &=A^{\mathsf T}H,
&D_c   &=A^{\mathsf T}CH,\\
M_A    &=A^{\mathsf T}\mathbf{1},
&M_{A,c}&=A^{\mathsf T}c,\\
O_A    &=A^{\mathsf T}A,
&O_{A,c}&=A^{\mathsf T}CA.
\end{aligned}
$$

Here $d,d_c\in\mathbb{N}^{t}$ count all/common reference SNPs in each atom;
$D,D_c\in\mathbb{R}^{b\times t}$ contain baseline-to-atom overlaps;
$M_A,M_{A,c}\in\mathbb{R}^{b}$ are baseline annotation counts; and
$O_A,O_{A,c}\in\mathbb{R}^{b\times b}$ are baseline overlap Gram matrices.
The corresponding reference-universe totals are $m$ and $\mathbf{1}^{\mathsf T}c$.

`atom_statistics.npz` stores $d,d_c,D,D_c$;
`baseline_statistics.npz` stores $M_A,M_{A,c},O_A,O_{A,c}$ and the two
reference-universe totals.

## 4. Full online indexed algorithm

Suppose the user supplies $Q$ focal gene-list files and optionally a fixed
control gene list.

### 4.1 Validate and resolve the inputs

The workflow loads one explicit complete index and validates its metadata,
chromosome coverage, shapes, sparse formats, column order, and `index_id`
binding. It then resolves each input token against `gene_catalog.parquet`.
Versioned Ensembl identifiers, aliases, duplicate entries, unresolved entries,
and genes excluded when the index was built follow the catalog-resolution
contract.

For a usable focal list $q$, resolution produces the Boolean catalog selector
$a_q$. Let $a_0$ denote the optional control selector.

### 4.2 Form Boolean unions of atoms

For each chromosome and focal query, compute

$$
z_q=\mathbf{1}[G^{\mathsf T}a_q>0].
$$

The comparison is elementwise. It implements a Boolean union, so duplicated,
overlapping, nested, or alias-selected genes do not multiply an atom's value.
The control selector is formed identically:

$$
z_0=\mathbf{1}[G^{\mathsf T}a_0>0].
$$

The implied SNP-level binary query annotation is

$$
q_q=Hz_q.
$$

Because atoms are disjoint, $q_q$ is exactly the direct padded-gene union on the
retained reference-SNP grid.

### 4.3 Assemble focal and control LD scores

Stack the focal selectors as

$$
Z=[z_1,\ldots,z_Q]\in\{0,1\}^{t\times Q}.
$$

All focal LD-score columns for the chromosome are obtained by one sparse matrix
product:

$$
\boxed{L_Q=YZ=PRHZ=PR[q_1,\ldots,q_Q]}.
$$

When enabled, the control column is

$$
\ell_0=Yz_0=PRq_0.
$$

Thus indexed and direct computation differ only in factorization and timing:
the expensive genotype/LD operation $PRH$ was performed once by the builder.
Online execution performs no PLINK, pairwise-$R^2$, or LD-window calculation.

### 4.4 Assemble annotation counts

For each focal query,

$$
M_q=d^{\mathsf T}z_q,
\qquad
M_{q,c}=d_c^{\mathsf T}z_q.
$$

The control counts are $d^{\mathsf T}z_0$ and $d_c^{\mathsf T}z_0$. Supplied
baseline counts come directly from $M_A$ and $M_{A,c}$. These quantities use the
LD-reference universe, not merely the $r$ persisted rows.

### 4.5 Assemble annotation overlaps

For two columns $u$ and $v$ of a generic reference-SNP annotation matrix $F$,
LDSC overlap is the Gram-matrix entry

$$
O_{uv}=\sum_{i=1}^{m} F_{iu}F_{iv}.
$$

For binary annotations this is the number of reference SNPs belonging to both
annotations. For continuous baseline annotations it is a weighted product and
need not be an integer.

The stored statistics give the required entries without reconstructing dense
query annotations:

$$
\begin{aligned}
\text{baseline--baseline:}\quad & O_A=A^{\mathsf T}A,\\
\text{baseline--query }q:\quad & Dz_q=A^{\mathsf T}Hz_q,\\
\text{query }q\text{ self-overlap:}\quad & d^{\mathsf T}z_q,\\
\text{baseline--control:}\quad & Dz_0,\\
\text{control self-overlap:}\quad & d^{\mathsf T}z_0,\\
\text{control--query }q:\quad & d^{\mathsf T}(z_0\land z_q).
\end{aligned}
$$

Replacing $d,D,O_A$ by $d_c,D_c,O_{A,c}$ gives the common-SNP overlaps.
The query self-overlap equals its count because the query annotation is binary.

The output overlap contract stores:

- rows for supplied baseline annotations and `gene_control` when enabled;
- columns for those same fixed annotations followed by every focal query; and
- a separate self-overlap entry for every focal query.

Focal query--query off-diagonal overlaps are intentionally not computed or
stored. `partitioned-h2` fits the fixed baseline/control block plus one focal
query per model, so no model requires an overlap between two different focal
queries. Their absence must not be interpreted as zero overlap.

### 4.6 Aggregate chromosomes

All score tables preserve the index's chromosome and persisted-row order. Counts,
overlap blocks, query diagonals, and reference-universe totals are additive:

$$
L_Q=\operatorname{row\_concat}_k L_Q^{(k)},
\qquad
M_q=\sum_k M_q^{(k)},
\qquad
O=\sum_k O^{(k)}.
$$

Chromosome-local atom identifiers never need to agree across chromosomes.

### 4.7 Prune unusable focal queries

After chromosome aggregation, a resolved focal query is retained only if

$$
M_q>0
$$

and its assembled values in $L_Q$ have nonzero variance across the persisted
regression rows. A zero-hit query receives status `zero_annotation_snps`; a
constant LD-score column receives status `zero_variance_ld_scores`. The workflow
removes a skipped focal query consistently from $L_Q$, its count record, its
fixed-row overlap column, and its query self-overlap. The fixed supplied
baseline and optional control are not focal-query pruning candidates.

If every focal query is skipped, the workflow writes query diagnostics but no
scientific LD-score tables and exits with an error.

### 4.8 Publish the canonical output

After query-status validation, the ordinary LD-score writer publishes:

| Artifact | Mathematical content |
| --- | --- |
| `ldscore.baseline.parquet` | persisted SNP identities, $w$, $L_A$, and optional $\ell_0$ |
| `ldscore.query.parquet` | persisted SNP identities and $L_Q$ |
| `ldscore.overlap.parquet` | all/common fixed-row overlap blocks and focal-query self-overlaps |
| `metadata.json` | annotation counts, SNP-universe totals/policies, columns, `index_id`, and provenance |
| `diagnostics/` | query resolution/status and workflow audit records |

Scientific matrix products and stored operator values use float64. The canonical
Parquet writer narrows LD-score columns to the public float32 storage dtype.
Negative adjusted-$r^2$ contributions are preserved, and no epsilon pruning or
post-assembly clamping is applied.

## 5. Why online assembly is fast

For $Q$ retained focal queries, the online numerical work is dominated by the
CSR--dense product $YZ$ and the atom-statistic products. In conventional sparse
notation their approximate arithmetic costs are

$$
O\!\left(\operatorname{nnz}(Y)Q\right)
\quad\text{and}\quad
O(btQ),
$$

respectively. Gene-to-atom rows are also sparse, so forming each Boolean union
touches the selected genes' stored memberships rather than the $m$-SNP genotype
grid. No online term depends on $n$ or requires genotype correlations.

The speedup comes from moving the expensive construction of $R$'s action on the
atom basis, $PRH$, offline and reusing it for every future gene list. The tradeoff
is index storage for $Y$ and its sufficient statistics, plus dense output memory
proportional to $rQ$. The current online implementation assembles all requested
query columns together and does not batch them.

## 6. Exactness identity

The complete equivalence argument for a focal list is

$$
\begin{aligned}
z_q &= \mathbf{1}[G^{\mathsf T}a_q>0],\\
q_q &= Hz_q,\\
\ell_q^{\mathrm{direct}} &= PRq_q,\\
\ell_q^{\mathrm{indexed}} &= Yz_q=(PRH)z_q=PR(Hz_q),\\
\therefore\quad
\ell_q^{\mathrm{indexed}}&=\ell_q^{\mathrm{direct}}.
\end{aligned}
$$

The equality is structural. Floating-point accumulation order can differ from a
separate direct run, so persisted float32 scores and downstream results are
compared under the validation tolerance rather than by requiring bitwise
identity. Annotation membership, identities, row/column order, binary counts,
and query statuses remain exact contracts. Continuous baseline overlap values
obey the same formulas and use the design specification's validated absolute
tolerance of $10^{-7}$.

## 7. File-to-symbol map

| Index component | Stored role |
| --- | --- |
| `gene_catalog.parquet` | catalog row identity and the information needed to form $a_q$ |
| `chromosomes/chrN/gene_to_atom.npz` | $G$ |
| `chromosomes/chrN/atoms.parquet` | genomic coordinates and ordering of the $t$ atoms |
| `chromosomes/chrN/ldscore_operator.npz` | $Y=PRH$ |
| `chromosomes/chrN/baseline_rows.parquet` | persisted identities, $w=PRp$, and $L_A=PRA$ |
| `chromosomes/chrN/atom_statistics.npz` | $d,d_c,D,D_c$ |
| `chromosomes/chrN/baseline_statistics.npz` | $M_A,M_{A,c},O_A,O_{A,c}$ and reference totals |
| root and chromosome `metadata.json` | shapes, ordering, scientific identity, policies, and provenance |

## 8. Related documentation and implementation

- [Exact gene LD-score index user guide](gene-ldscore-index.md)
- [Gene-list indexed `ldscore` tutorial](../wiki/main-functionalities/ldscore.md)
- [Detailed design specification](../specs/2026-08-03-exact-disjoint-atom-gene-ldscore-index-design.md)
- [General LD-score accumulation mathematics](ldscore-parquet-accumulation.md)
- Builder and online orchestration:
  [`src/ldsc/gene_ldscore_index.py`](../../src/ldsc/gene_ldscore_index.py)
- Atom construction and matrix kernels:
  [`src/ldsc/_kernel/gene_ldscore_index.py`](../../src/ldsc/_kernel/gene_ldscore_index.py)
