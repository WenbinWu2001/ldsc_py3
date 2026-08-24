# Conventional Versus Partitioned LDSC Heritability

Last updated on: 2026-08-24

Conventional one-column LDSC and multi-annotation partitioned LDSC target the
same genome-wide SNP heritability, but they do not generally produce identical
numerical estimates.
The two approaches make different assumptions about how heritability is
distributed across SNPs, fit different regression models, and can use different
estimation defaults.

## 1. Different models of per-SNP heritability

### Conventional one-column LDSC

Conventional LDSC fits the approximate model

```text
E[chi_j^2] = intercept + N_j * (h_g^2 / M) * l_j,
```

where `l_j` is the ordinary LD score of SNP `j`, `M` is the number of reference
SNPs, and `h_g^2 / M` is a single average per-SNP heritability parameter.

**The conventional model assumes that heritability is uniformly distributed
across the reference SNPs in expectation.**
This does not mean that every realized SNP effect has the same magnitude.
Rather, the working model assigns every SNP the same expected contribution to
heritability and uses one LD-score coefficient to summarize the genome.

### Multi-annotation partitioned LDSC

Partitioned LDSC fits

```text
E[chi_j^2] = intercept + N_j * sum_c tau_c * l(j, c),
```

where `l(j, c)` is the annotation-specific LD score and `tau_c` is the
conditional per-SNP heritability contribution associated with annotation `c`.
The model defines the expected per-SNP heritability as

```text
Var(beta_k) = sum_c a(k, c) * tau_c,
```

where `a(k, c)` records SNP `k`'s value for annotation `c`.

**The partitioned model allows per-SNP heritability to vary according to the
SNP's modeled importance, as defined by the included annotations.**
For example, SNPs in coding regions, promoters, enhancers, conserved regions,
or other functional categories can have different expected per-SNP
heritability from SNPs outside those categories.
The distribution is therefore not uniform in general.
The model can represent only the forms of importance encoded by the supplied
annotations, so annotation selection determines which departures from
uniformity it can capture.

The partitioned estimate of total heritability is

```text
h_g^2(partitioned) = sum_c M_c * tau_c,
```

where `M_c` is the reference-SNP count or sum associated with annotation `c`.
This total aggregates the annotation-specific fitted contributions.

## 2. Why the point estimates can differ

### 2.1 The one-column model is a constrained model

Suppose two annotations divide the genome into categories `A` and `B`.
The conventional model effectively imposes a common per-SNP contribution,

```text
tau_A = tau_B = tau.
```

The partitioned model allows

```text
tau_A != tau_B.
```

Even when the total LD score is the sum of the category-specific LD scores, the
slope from a regression on that sum does not generally equal the total formed
from several freely estimated slopes.
The conventional estimate comes from a constrained regression, whereas the
partitioned estimate comes from a multivariable regression and then aggregates
the fitted annotation coefficients.
Relaxing the common-effect constraint changes the weighted least-squares
projection and can change the estimated total.

The estimates should become similar when the uniform model is adequate, or when
the annotation-specific coefficients happen to be similar.
They can differ when functional categories have meaningfully different
expected per-SNP contributions.

### 2.2 The one-column model can average over annotation-specific architecture

If functional regions have higher per-SNP heritability than the rest of the
genome, the conventional model cannot represent this structure explicitly.
Its single coefficient becomes a weighted genome-wide average of the different
annotation-specific effects.

The partitioned model can assign different conditional contributions to the
included categories.
Adding these regressors changes how the model attributes variation in the GWAS
chi-square statistics, which can change both the individual coefficients and
their aggregate total.
This is the same general phenomenon that causes a coefficient to change when a
relevant predictor is added to a regression model.

The richer model is not automatically correct.
If important determinants of per-SNP heritability are absent from the
annotation set, the fitted coefficients and total can still reflect model
misspecification.

### 2.3 The regressions use different weighting behavior

LDSC uses weighted regression because GWAS chi-square statistics are
heteroscedastic and SNPs in LD are not independent.
The conventional and partitioned workflows do not necessarily give each SNP
the same influence.

In the current implementation, a single-annotation model uses the iterative
weighting path, whereas a multi-annotation model uses the legacy partitioned
LDSC weighting path.
Different weights produce different weighted least-squares projections, so the
estimated totals can differ even when both analyses begin with the same summary
statistics.

### 2.4 Default chi-square filtering differs

The current workflow applies different default outlier policies:

- A single-annotation model has no default global chi-square cap and uses the
  two-step estimator with a default cutoff of `30` when the intercept is free.
- A genuine multi-annotation model does not use the two-step estimator and
  applies the default cap `max(0.001 * N.max(), 80)` when `--chisq-max` is not
  specified.

Consequently, conventional and partitioned analyses can fit different sets of
regression SNPs unless the filtering choices are explicitly harmonized.
A small number of high-chi-square SNPs can have substantial leverage, so this
difference can affect the total heritability estimate.

### 2.5 Intercept estimation interacts with model specification

The LDSC intercept and the LD-score regressors jointly explain variation in the
GWAS chi-square statistics.
Replacing one LD-score column with several annotation-specific columns can
change the fitted intercept as well as the slopes.
Because total heritability is calculated from those slopes, the change can
propagate into the total estimate.

Both analyses should therefore use the same intercept policy for a controlled
comparison.
The intercept should either be estimated in both analyses or constrained to the
same value in both analyses.

### 2.6 Multi-annotation estimates can be less stable

Annotation-specific LD-score columns can be strongly correlated.
The multivariable regression must then separate contributions that produce
similar LD-score patterns.
This can lead to large coefficient standard errors, strong covariance among
coefficients, and sensitivity to small changes in SNP selection or weights.

The aggregate total can be more stable than the individual annotation
coefficients because correlated coefficient errors may partially cancel.
Nevertheless, the aggregate is not algebraically constrained to equal the
one-column estimate.

## 3. When should the estimates be close?

The conventional and partitioned estimates should generally be reasonably
close when:

- both analyses use the same GWAS summary statistics and regression SNPs;
- both use the same reference-SNP universe and `M` definition;
- chi-square filtering and intercept settings are harmonized;
- the annotations describe the same overall reference-SNP universe;
- the multi-annotation design is well-conditioned;
- sample size is large; and
- both the uniform and annotation-informed working models approximate the
  genetic architecture adequately.

The estimates can differ substantially when the trait has strong
annotation-specific architecture, the annotations are highly correlated, a few
large association statistics influence the fit, or the workflows use different
filtering or intercept settings.

## 4. How to compare the results

The appropriate comparison is

```text
h2.total_h2_obs
versus
partitioned_h2.total_h2_obs.
```

Do not compare conventional total heritability with `category_h2_obs`.
The latter is the conditional contribution associated with one annotation, not
the model's genome-wide total.

Compare the point estimates together with their standard errors and diagnostic
settings.
At minimum, verify the following quantities:

1. The post-filter number of regression SNPs.
2. The reference-SNP count kind, such as common versus all SNPs.
3. The effective chi-square cutoff.
4. Whether the intercept was estimated or constrained.
5. The LD-score annotations included in each model.

Exact equality is not a general LDSC identity.
It requires special constraints that make the multi-annotation model behave
like the one-column model, or it occurs by numerical coincidence.
Under compatible and adequately specified models, both estimates target the
same biological total and should approach that total as statistical information
increases.
