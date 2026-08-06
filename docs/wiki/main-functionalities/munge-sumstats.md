# Munge summary statistics

`ldsc munge-sumstats` converts a raw GWAS table into LDSC3-ready summary
statistics. See the current [munge-sumstats guide](../../current/munge-sumstats.md)
for the full workflow and output contract.

## Sample-size columns

Use one per-variant sample-size strategy:

- Direct N: `--N-col <column>`.
- Case/control counts: `--N-cas-col <cases> --N-con-col <controls>`.

The case and control flags are a required pair, and they cannot be combined
with `--N-col`. An explicit direct-N choice suppresses automatically inferred
case/control columns with a warning; an explicit case/control choice suppresses
an automatically inferred direct-N column with a warning.

### How the case/control pair becomes N

LDSC3 preserves the LDSC2 legacy case/control normalization. For variant
`i`, let `T_i = NCAS_i + NCON_i` and `p_i = NCAS_i / T_i`. Let `p_ref` be the
mean `p_i` among variants having the maximum `T_i`. The munger calculates:

```text
N_i = T_i * p_i / p_ref
```

This is **not** the commonly used harmonic effective sample size
`4 / (1/NCAS_i + 1/NCON_i)`. If the case fraction is constant, the legacy
calculation reduces to `NCAS_i + NCON_i`; if it varies, N is scaled relative
to the case fraction among the maximum-total-count variants.

Choose the strategy according to the scientific meaning of the source fields:

- Use `--N-col NEFF` when the data producer documents `NEFF` as the effective
  sample size corresponding to the reported association statistic—for example,
  after accounting for case/control imbalance, per-variant missingness,
  meta-analysis participation, or analysis weights—and those exact values are
  intended for LDSC.
- Use `--N-cas-col NCAS --N-con-col NCON` when the per-variant counts are the
  trusted inputs and LDSC2-compatible legacy case/control normalization is the
  intended rule.

`NEFF` definitions vary between producers. Check the study documentation and
compare `NEFF` with the count-derived values rather than assuming the two
strategies are interchangeable.

If the header automatically maps both strategies, such as `N + NCAS + NCON`,
LDSC3 stops rather than allowing one value to overwrite the other. Rerun with
either:

```bash
--N-col N
```

or:

```bash
--N-cas-col NCAS --N-con-col NCON
```

For `NEFF + NCAS + NCON`, use `--N-col NEFF` when effective N is the intended
quantity. The inferred `NCAS`/`NCON` columns are then suppressed automatically;
`--ignore` is unnecessary.
