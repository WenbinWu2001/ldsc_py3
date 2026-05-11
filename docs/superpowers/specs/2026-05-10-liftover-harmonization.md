# Liftover Conflict Harmonization — Design Proposal

## Context

Liftover collision handling is currently asymmetric across the two workflows
that perform coordinate liftover:

- **Sumstats munger** (`src/ldsc/sumstats_munger.py` + `_kernel/liftover.py`)
  has no public duplicate policy and unconditionally uses `drop-all` for both
  source and target duplicate `CHR/POS` groups. It logs counts and up to five
  examples at `INFO`. It writes no `dropped_snps/` sidecar.
- **Reference-panel builder** (`src/ldsc/ref_panel_builder.py`) exposes the
  legacy public flag `--duplicate-position-policy {error,drop-all}` with config
  field `ReferencePanelBuildConfig.duplicate_position_policy` (default
  `drop-all`). It writes `dropped_snps/chr{chrom}_dropped.tsv.gz` with schema
  `CHR, SNP, source_pos, target_pos, reason` containing **duplicate drops
  only**. Unmapped and cross-chromosome drops are logged but not in the
  sidecar.

The desired shared behavior is `drop-all` everywhere (every row in a
duplicate-coordinate cluster is dropped). The implementation already routes
both workflows through shared kernel helpers
(`duplicate_coordinate_drop_result`, `liftover_drop_report`,
`log_liftover_drop_report`) in `_kernel/liftover.py`, so the remaining work is
contract harmonization, not algorithmic.

This proposal collapses the policy surface, harmonizes audit output between
the two workflows, and quiets logs without losing debuggability.

---

## Current Behavior (Snapshot)

| Concern | Sumstats munger | Reference-panel builder |
| --- | --- | --- |
| Duplicate policy public knob | none (always `drop-all`) | `--duplicate-position-policy {error,drop-all}` |
| Source duplicate handling | drop-all before mapping | drop-all (or hard error) before mapping |
| Target duplicate handling | drop-all after mapping | drop-all (or hard error) after mapping |
| Sidecar for duplicates | none | `dropped_snps/chr{chrom}_dropped.tsv.gz`, duplicate-only |
| Sidecar for unmapped/cross-chrom | none | none (logs only) |
| All-dropped behavior | hard error (whole run) | skip chromosome; whole run errors only if zero chromosomes emit |
| Example logging | INFO with up to 5 examples | INFO/WARNING with up to 5 examples |
| HM3 quick liftover | yes | no |
| Liftover in `rsid` mode | hard error pre-IO | rejected when matching chain set |

---

## Recommended Behavior

### R1 — Single duplicate policy: `drop-all` everywhere

Drop the user-facing duplicate-position policy entirely. `drop-all` becomes the
only behavior for both workflows.

**Rationale.** The `error` mode produces a hard abort with a textual list of
collisions. The same information is available — and machine-readable — through
the `drop-all` sidecar plus log counts, so `error` adds no real diagnostic
power once sidecars cover all reasons (R2, R3). Eliminating the knob removes a
public contract divergence between the two workflows and matches the
already-locked sumstats behavior.

### R2 — One harmonized `dropped_snps/` sidecar shape, written by both workflows

The sidecar is **always written** by every run, even when zero rows were
dropped. Empty runs produce a gzipped TSV containing only the header line.
This eliminates the entire "stale sidecar from a prior run" problem class:
the path is always-owned and always-produced, so a re-run simply overwrites
the prior file in place and there is no deletion semantics to manage.

- **Sumstats** writes `{output_dir}/dropped_snps/dropped.tsv.gz`. One file
  per run; chromosome is a column in the file.
- **Reference-panel** writes `{output_dir}/dropped_snps/chr{chrom}_dropped.tsv.gz`
  for every chromosome processed by the current run. Chromosomes with zero
  drops produce a header-only file. Chromosomes *not* processed by the
  current run are not touched (see the "class-2 artifacts" note below).

Rationale for the always-written contract: the file's presence is a
positive attestation that the audit ran. Today's "no file = no drops"
convention is ambiguous (could mean the audit was skipped, or the run
crashed before writing, or there were truly no drops). An always-present
file with a known schema lets downstream consumers always
`pd.read_csv(sidecar, dtype=...)` without an `if path.exists()` branch.

**Schema (column-stable for ref-panel readers; dtype contract widened):**

| Column | Pandas dtype | Notes |
| --- | --- | --- |
| `CHR` | `string` | Nullable. NA only for sumstats `missing_coordinate` rows that lacked `CHR`. |
| `SNP` | `string` | Source label / rsID, preserved unchanged. Always present. |
| `source_pos` | `Int64` | Nullable. NA only for sumstats `missing_coordinate` rows that lacked `POS`. |
| `target_pos` | `Int64` | Nullable. NA for `missing_coordinate`, `source_duplicate`, `unmapped_liftover`, `cross_chromosome_liftover`. Populated for `target_collision`. |
| `reason` | `string` | One of the controlled vocabulary in R3. |

Nullable types are required because sumstats `missing_coordinate` rows can
have NA `CHR`, NA `POS`, or both (rows are dropped *because* their
coordinates are missing). The TSV is gzip-compressed; pandas writes empty
strings for nullable-NA cells with `na_rep=""`.

**Reader contract.** Consumers must pass an explicit `dtype=` map when
calling `pd.read_csv`:

```python
pd.read_csv(
    sidecar_path,
    sep="\t",
    compression="gzip",
    dtype={
        "CHR": "string",
        "SNP": "string",
        "source_pos": "Int64",
        "target_pos": "Int64",
        "reason": "string",
    },
)
```

Pandas cannot infer `string` / `Int64` extension types from an empty TSV
(header-only sidecars from clean runs), and integer inference from a
`source_pos` column with empty cells would default to `float64`. The
explicit dtype map is the only way to round-trip header-only sidecars
back to the canonical schema. The writer is responsible only for emitting
valid headers and proper NA encoding (`na_rep=""`); the consumer owns the
dtype map.

For ref-panel rows, the column **names** remain unchanged and `source_pos`
continues to be structurally non-null because BIM `BP` is non-null. The
unified dtype contract still widens both position columns to nullable `Int64`.
This is required because `target_pos` is intentionally empty for ref-panel
`source_duplicate`, `unmapped_liftover`, and `cross_chromosome_liftover` rows,
and clean runs may produce header-only sidecars. Readers that want a stable
parser across both workflows must use the explicit `Int64`/`string` dtype map
above rather than plain NumPy integer dtypes.

### R3 — Controlled vocabulary of drop reasons in sidecar

The shared vocabulary covers every reason that could remove a row inside
the liftover stage of either workflow:

- `missing_coordinate` — row lacked usable `CHR` or `POS` before mapping.
- `source_duplicate` — duplicate source `CHR/POS` group (drop-all).
- `unmapped_liftover` — chain returned no hit, or HM3 map had no key.
- `cross_chromosome_liftover` — chain produced only cross-chromosome hits.
- `target_collision` — multiple source rows mapped to the same target
  `CHR/POS` (drop-all).

**Per-workflow reason coverage** (intentionally asymmetric):

| Reason | Sumstats | Reference-panel |
| --- | --- | --- |
| `missing_coordinate` | yes | **never** — PLINK `.bim` `BP` is structurally non-null; malformed values raise during parsing before sidecar logic runs |
| `source_duplicate` | yes | yes (chr_pos mode only) |
| `unmapped_liftover` | yes (when liftover requested) | yes (when matching chain provided) |
| `cross_chromosome_liftover` | yes (chain method only) | yes (chain method only) |
| `target_collision` | yes (when liftover requested) | yes (chr_pos mode only) |

The vocabulary is shared so consumers can write one downstream parser; not
every reason fires in every workflow.

### R4 — Quiet logs by default; examples behind `DEBUG`

At `INFO`/`WARNING`, log only counts and a single pointer to the sidecar:

```text
Reference-panel liftover dropped 142 SNPs on chromosome 6
  (12 source_duplicate, 84 target_collision, 31 unmapped_liftover, 15 cross_chromosome_liftover).
  See 'out/dropped_snps/chr6_dropped.tsv.gz' for row-level provenance.
```

Per-reason example rows (currently `liftover_drop_report.examples`) move to a
single `DEBUG` line per reason, e.g.:

```text
DEBUG  liftover example dropped (target_collision): {SNP: rsX, CHR: 6, source_POS: 100, target_POS: 200}
```

This addresses the "examples at normal log levels makes logs noisy" concern
while keeping the data accessible via `--log-level DEBUG`.

### R5 — Keep HM3 quick liftover sumstats-only

Locked. No change. HM3 curated map is a sparse coordinate dictionary tailored
to HM3 SNPs and does not cover the dense PLINK universe of a reference panel.

### R6 — All-dropped contract stays workflow-specific

Locked. Sumstats raises if liftover removes every row (single output,
catastrophic). Ref-panel skips the chromosome and only fails the run when
zero chromosomes emit (per-chromosome outputs, partial recovery is the
documented behavior). Harmonization here would change the ref-panel
output contract.

### R7 — Per-workflow contracts retained for class-2 (conditionally-produced) artifacts; user owns directory hygiene

**Governing principle.** For class-2 artifacts, the workflow follows three
explicit rules:

1. **Don't touch what we don't recognize.** Files outside the current
   run's expected-output set are never deleted, modified, or moved by the
   workflow. This is true under `--overwrite=True` as well — the
   overwrite flag opts into replacing files this run *will produce*, not
   into deleting unrelated siblings.
2. **Document the limitation as a user contract.** The user is
   responsible for keeping the output directory internally consistent.
   The workflow's job is to inform; the user's job is to act.
3. **Inform loudly when divergence is possible.** When the workflow
   detects pre-existing class-2 artifacts that the current run will not
   touch, it should emit a `WARNING` naming the affected files and the
   silent-mix risk. The warning is informational — it must not abort
   the run, because chromosome-by-chromosome batch processing with
   consistent configuration is a legitimate workflow that produces the
   same on-disk pattern as the dangerous silent-mix case, and the
   workflow has no way to distinguish the two cases on its own.
   The warning behavior is implemented by the follow-up stale-class-2 design
   linked below.

This is the user-contract model used by `make`, `bazel`, and most data
pipeline tools: the workflow guarantees consistency for what it
*produces* in a single invocation, not for the directory's full state
across invocations. Cross-invocation consistency is the user's
responsibility.



Two distinct artifact classes exist in both workflows:

- **Class 1 — always-produced.** Workflow always writes this file regardless
  of run configuration. Examples: `sumstats.metadata.json`, `sumstats.log`,
  the new `dropped_snps/dropped.tsv.gz`, and the new ref-panel
  `dropped_snps/chr{chrom}_dropped.tsv.gz` (per R2, always written for every
  chromosome processed by the current run; empty header-only TSV when no
  liftover-stage drops occurred).
- **Class 2 — conditionally-produced.** The set of files written depends on
  per-run flags. Sumstats: `sumstats.parquet` vs `sumstats.sumstats.gz`
  depending on `--output-format`. Ref-panel: per-chromosome
  `chr{N}_r2.parquet` and `chr{N}_meta.tsv.gz` only for chromosomes
  processed by the current run. Per-build artifacts under `hg19/` vs
  `hg38/` only when a matching liftover chain is supplied.

**Workflow-log asymmetry.** Both workflows produce a `.log` file under their
parsed-CLI entry points (`sumstats.log` from `SumstatsMunger.run`;
`build-ref-panel.log` from `run_build_ref_panel_from_args`). Direct
programmatic invocation of `ReferencePanelBuilder.run()` leaves
`self._workflow_log_path = None` (`ref_panel_builder.py:170`), and the
`workflow_logging` context at `ref_panel_builder.py:239` therefore writes
no log file. So the "always-produced" framing for `.log` files is
**CLI/parsed-workflow scoped**, not "any invocation of the workflow object."
Direct API consumers who want a log can set `_workflow_log_path` themselves
(or use the parsed-CLI entry point). This pre-existing distinction is
unchanged by liftover harmonization; the `dropped_snps/` sidecar work does
not introduce any new logging path.

R2 unifies the new sidecar across both workflows by making it class-1 plus
the always-written-empty pattern. **Class-2 artifacts retain their existing
per-workflow contracts**:

| Workflow | Helper used | Stale-sibling behavior under `--overwrite` |
| --- | --- | --- |
| `munge-sumstats` | `preflight_output_artifact_family` + `remove_output_artifacts` | Auto-removed after successful write. Example: switching `--output-format` from `both` to `parquet` removes the prior `sumstats.sumstats.gz`. |
| `build-ref-panel` | `ensure_output_paths_available` only | **Not** removed. Stale per-chromosome or per-build siblings from prior runs persist on disk. Documented contract (post-harmonization wording): *"Use a fresh output directory when changing emitted builds, liftover/coordinate configuration, or chromosome scope"* (`ref_panel_builder.py:21-26`). The pre-harmonization wording quoted "duplicate-position policy" as one of the trigger conditions, but that public knob is removed by this work; the module docstring update is part of the implementation (Task 5). |

**Why the asymmetry is preserved.** The two contracts emerged from a real
cost asymmetry between the workflows. Sumstats produces small artifacts
quickly, so auto-cleanup of `output_format` siblings is convenient and
the cost of an accidental deletion is small. Ref-panel artifacts can take
hours of CPU and reach hundreds of MB per chromosome; auto-deleting "stale"
siblings under `--overwrite` would risk destroying expensive compute on a
mistyped invocation. The defensive answer was "don't touch what we don't
recognize, document the limitation, and let the user manage directory
hygiene." Both choices are defensible. Changing either one means
re-litigating that cost tradeoff and is out of scope for liftover
harmonization.

**Concrete example of the ref-panel "fresh dir" contract.** Run 1 builds
chr1-22 into `out/`. Run 2 reuses `out/` with `--plink-prefix only_chr6_7.@
--overwrite`. The Run-2 preflight list is constructed by
`_expected_ref_panel_output_paths` and contains only chr6/chr7 paths
(loop at `ref_panel_builder.py:152` iterates over the chromosomes resolved
from the current PLINK input, not all 22). The collision check passes for
those two; chr1-5 and chr8-22 files from Run 1 are never inspected. Result:
`out/hg19/` ends up holding chr6/chr7 from Run 2 and chr1-5,8-22 from Run 1
side-by-side, silently mixing two different builds. Downstream tools that
load the directory as a "complete reference panel" see 22 chromosomes
without any indication that they came from different inputs. The runtime
metadata sidecars (`chr*_meta.tsv.gz`) record the per-chromosome provenance,
but most consumers do not cross-check.

This is the silent-mix failure mode the existing module docstring and warning
message describe. Liftover harmonization does **not** clean it automatically.
Users running
`build-ref-panel` against the same directory more than once must continue
to either:

1. start from a fresh `--output-dir`, or
2. manually clean `out/{build}/chr*_r2.parquet`,
   `out/{build}/chr*_meta.tsv.gz`, and `out/dropped_snps/chr*_dropped.tsv.gz`
   for chromosomes outside the current run's scope.

A future spec may revisit this contract, but it is independent of the
liftover work.

**Stale class-2 sibling warning — implemented in follow-up.** To make
the user contract observable rather than implicit, ref-panel emits a single
`WARNING` log line per run when pre-existing class-2 artifacts in `output_dir`
are not in the current run's expected set. That behavior is implemented
separately in
[`docs/superpowers/specs/2026-05-11-ref-panel-stale-class2-warning.md`](./2026-05-11-ref-panel-stale-class2-warning.md)
with its own plan at
[`docs/superpowers/plans/2026-05-11-ref-panel-stale-class2-warning.md`](../plans/2026-05-11-ref-panel-stale-class2-warning.md).
The principles in this section (R7's three governing rules) are the
contract that warning enforces; the split kept the liftover change set focused
on liftover while still making the class-2 hygiene contract observable.

---

## Public Contract Changes

| Change | Workflow | Magnitude | Approval needed? |
| --- | --- | --- | --- |
| Remove `--duplicate-position-policy` CLI flag | ref-panel | breaking CLI removal | **YES — explicit user approval before implementation** |
| Remove `ReferencePanelBuildConfig.duplicate_position_policy` field | ref-panel | breaking Python API removal | **YES — explicit user approval** |
| `error` policy ceases to exist (always `drop-all`) | ref-panel | behavior change for users currently relying on hard-fail | YES — implied by above |
| Add `dropped_snps/dropped.tsv.gz` sidecar | sumstats | additive new artifact, always written (header-only when no drops) | minor — additive |
| Switch ref-panel `dropped_snps/chr{chrom}_dropped.tsv.gz` from conditional-write to always-write | ref-panel | contract change: file presence no longer signals "had drops"; consumers must read row count instead | YES — confirmed user direction (option R1) |
| Expand sidecar `reason` vocabulary (additive) | ref-panel | additive: column names unchanged; dtype contract widened per R2; more rows | minor — additive |
| Move per-row example log lines from `INFO` to `DEBUG` | both | log volume reduction, no API change | minor — log only |

No changes to: `MungeConfig`, `GlobalConfig`, sumstats `metadata.json`
sidecar shape, ref-panel runtime metadata TSV/parquet schemas, or any
liftover validation rule already locked in
`docs/current/liftover-harmonization-decisions.md`.

---

## Output / Logging Contract

### Sumstats (`munge-sumstats`)

```
{output_dir}/sumstats.parquet           # class 2 — produced if output_format ∈ {parquet, both}
{output_dir}/sumstats.sumstats.gz       # class 2 — produced if output_format ∈ {tsv.gz, both}
{output_dir}/sumstats.metadata.json     # class 1 — always produced
{output_dir}/sumstats.log               # class 1 (CLI/parsed entry) — always produced via SumstatsMunger.run
{output_dir}/dropped_snps/dropped.tsv.gz   # NEW class 1 — always produced (header-only when no drops)
```

### Reference-panel (`build-ref-panel`)

```
{output_dir}/{build}/chr{chrom}_r2.parquet                  # class 2 — per-chromosome, per-build
{output_dir}/{build}/chr{chrom}_meta.tsv.gz                 # class 2 — per-chromosome, per-build
{output_dir}/dropped_snps/chr{chrom}_dropped.tsv.gz         # NEW class 1 — always produced for each chromosome processed (header-only when no drops); broader reason coverage
{output_dir}/build-ref-panel.log (or .chr{chrom}.log)       # class 1 (CLI/parsed entry only) — produced when run_build_ref_panel_from_args sets _workflow_log_path; direct ReferencePanelBuilder.run() leaves it None and no log file is written
```

### Log behavior summary

- `DEBUG`: up to 5 example rows per reason (existing `examples` payload).
- `INFO`: per-stage counts, total dropped, breakdown by reason, sidecar path
  pointer.
- `WARNING`: large drop fractions or "all rows dropped on chromosome X"
  conditions (existing thresholds).
- Sidecar is the source of truth for row-level audit; logs do not duplicate
  it at default verbosity.

---

## Backward Compatibility Impact

**Breaks (require explicit user approval):**

1. CLI users invoking `build-ref-panel --duplicate-position-policy error`
   will see an "unrecognized argument" error.
2. Python callers passing
   `ReferencePanelBuildConfig(duplicate_position_policy="error")` will see a
   `TypeError` (unexpected keyword argument) once the field is removed.
3. Anyone relying on a build to *abort* on duplicates must instead inspect
   the `dropped_snps/` sidecar after the run to detect collisions.

**Non-breaking (additive or behind `DEBUG`):**

- New sumstats `dropped_snps/dropped.tsv.gz` artifact (always written, even
  on clean runs — header-only TSV when no rows dropped). Output preflight
  declares it unconditionally in both `produced_paths` and `owned_paths`.
- Ref-panel `dropped_snps/chr{chrom}_dropped.tsv.gz` sidecar switches from
  conditional-write to always-write (per spec R2 confirmation 2026-05-11),
  produced for every chromosome the run processes (header-only when no
  drops). Sidecar gains rows with new `reason` values for chain-mapping
  drops; reference-panel-applicable reasons are
  `source_duplicate`/`unmapped_liftover`/`cross_chromosome_liftover`/`target_collision`
  (no `missing_coordinate` — see R3 per-workflow coverage table).
  Consumers that previously used `path.exists()` to detect "had drops"
  must switch to checking the row count.
- Default log verbosity drops; previous example payloads remain fully
  available via `--log-level DEBUG`.

**Unchanged:**

- Sumstats metadata sidecar shape (`format`, `trait_name`,
  `config_snapshot`).
- Ref-panel runtime metadata TSV/parquet schemas.
- All locked liftover validation rules.
- `validate_config_compatibility` and downstream regression behavior.

---

## Test / Doc Updates Needed

### Tests

- Remove `tests/...test_ref_panel_builder...` cases that exercise
  `duplicate_position_policy="error"` and the corresponding CLI flag.
- Add a ref-panel test that asserts collisions populate the sidecar (with
  the new `reason` values) and the run continues without aborting.
- Add a sumstats test that asserts `dropped_snps/dropped.tsv.gz` is written
  with rows covering all five sumstats-applicable reasons:
  `missing_coordinate`, `source_duplicate`, `unmapped_liftover`,
  `cross_chromosome_liftover`, `target_collision`.
- Add a sumstats test that asserts the sidecar is written as a header-only
  TSV (no data rows) when the run drops nothing, matching the always-write
  contract in R2.
- Add a logging test that confirms example rows appear at `DEBUG` and not
  `INFO` after the threshold change (simple `caplog` check).
- Update preflight/overwrite tests for sumstats to include the new sidecar
  in the owned-artifact family.
- Snapshot test that sumstats sidecar JSON shape is unchanged (regression
  guard against accidental scope creep).

### Docs

- `docs/current/liftover-harmonization-decisions.md`: update the
  reference-panel section to remove `error` policy mention and document the
  expanded sidecar reason vocabulary plus the new sumstats sidecar.
- `docs/current/data_flow.md`, `docs/current/architecture.md`: reflect the
  new sumstats artifact in the munge stage diagram.
- `tutorials/munge_sumstats.md` (or equivalent): add a one-paragraph
  "auditing dropped SNPs" note pointing to the sidecar.
- `docs/superpowers/specs/2026-05-09-munge-sumstats-liftover.md`: append a
  follow-up addendum recording the policy collapse and sumstats sidecar
  addition.

### Code touchpoints (for the eventual implementation plan, not this design)

- `src/ldsc/config.py`: remove `duplicate_position_policy` field and its
  validation branch on `ReferencePanelBuildConfig`.
- `src/ldsc/ref_panel_builder.py`: remove the `--duplicate-position-policy`
  parser entry, drop the `policy=` plumbing through `_resolve_unique_snp_set`
  (always `drop-all`), broaden the dropped-frame to include
  `unmapped_liftover` and `cross_chromosome_liftover` rows that today only
  hit the log (NB: `missing_coordinate` does **not** apply to ref-panel —
  PLINK BIM `BP` is structurally non-null per R3 per-workflow coverage
  table), and downgrade the example-emitting log calls.
- `src/ldsc/sumstats_munger.py` + `src/ldsc/_kernel/liftover.py`: collect
  the existing per-reason `LiftoverDropReport` payloads in
  `apply_sumstats_liftover` into a single drop frame, write
  `dropped_snps/dropped.tsv.gz` from the workflow layer (matching the
  ref-panel writer style), include the path in preflight/overwrite, and
  switch example logging to `DEBUG`.
- `src/ldsc/_kernel/liftover.py::log_liftover_drop_report`: change the
  `INFO` line to count-only with a sidecar pointer placeholder; emit
  examples on a separate `DEBUG` line.

---

## Verification Plan (post-implementation)

1. `pytest tests/_kernel/test_liftover_module.py
   tests/_kernel/test_sumstats_munger_liftover_stage.py
   tests/test_munge_config_liftover.py
   tests/test_sumstats_munger_liftover.py
   tests/test_ref_panel_builder*.py -v`
2. CLI smoke:
   `ldsc build-ref-panel --help | grep -v duplicate-position-policy` (must
   show no match) and `ldsc munge-sumstats --help` (unchanged).
3. End-to-end: run `build-ref-panel` on a fixture with intentional
   collisions across the four ref-panel-applicable reason categories
   (`source_duplicate`, `unmapped_liftover`, `cross_chromosome_liftover`,
   `target_collision`); confirm `dropped_snps/chr*_dropped.tsv.gz`
   contains rows from each reason and the `INFO` log shows counts only.
   Confirm a chromosome with no drops produces a header-only sidecar.
4. End-to-end: run `munge-sumstats` with a chain liftover request on a
   fixture seeded with all five sumstats-applicable categories: missing
   coordinates, source duplicates, unmapped chain hits, cross-chromosome
   hits, and target collisions; confirm `dropped_snps/dropped.tsv.gz`
   contains rows from each reason. Confirm a clean run produces a
   header-only sidecar.
5. Re-run `--log-level DEBUG` on (3)/(4) and confirm example rows appear
   for each reason.

---

## Confirmed Decisions

User confirmed (2026-05-10):

1. **Remove `--duplicate-position-policy` outright.** Delete the CLI flag and
   the `ReferencePanelBuildConfig.duplicate_position_policy` field. `drop-all`
   becomes the only behavior. This is a breaking change for any caller
   currently passing `error`; downstream auditing happens through the
   sidecar.
2. **Sidecars cover the workflow-applicable reasons in both workflows.**
   Sumstats gains `dropped_snps/dropped.tsv.gz` covering all five
   reasons (`missing_coordinate`, `source_duplicate`,
   `unmapped_liftover`, `cross_chromosome_liftover`, `target_collision`).
   Ref-panel sidecar gains `unmapped_liftover` and
   `cross_chromosome_liftover` alongside its existing
   `source_duplicate`/`target_collision` rows. Ref-panel never
   produces `missing_coordinate` rows — PLINK BIM `BP` is structurally
   non-null. Schema columns unchanged; dtypes widened to nullable
   `string`/`Int64` so both workflows share one schema.
3. **Examples move from `INFO` to `DEBUG`.** Default logs become count-only
   with a sidecar pointer; row-level examples remain available via
   `--log-level DEBUG`.

User confirmed (2026-05-11):

4. **Always-write the new sidecar, in both workflows (option R1).** Replaces
   the earlier "conditionally-produced" wiring. The sidecar is class-1
   (always-owned, always-produced); empty runs write a header-only TSV.
   Eliminates stale-cleanup complexity for this artifact entirely. No
   `Path.unlink()` is required for the new sidecar in either workflow.
5. **Class-2 (conditionally-produced) artifacts retain their existing
   per-workflow contracts.** Sumstats keeps `output_format`-sibling
   auto-cleanup via `preflight_output_artifact_family` +
   `remove_output_artifacts`. Ref-panel keeps the documented "use a fresh
   output directory" contract (no auto-cleanup) for per-chromosome and
   per-build artifacts. The asymmetry is rooted in cost differences
   between the workflows and is out of scope for liftover harmonization.
   See R7 for the rationale and the silent-mix failure mode users must
   manage themselves.
6. **The stale-class-2 `WARNING` implementation ships as a separate
   spec.** R7's principles (governing rules and silent-mix failure mode)
   stay in this spec because they shape the per-workflow contracts the
   liftover harmonization preserves. The warning implementation — a directory
   scan, a log call, and focused tests — lives at
   `docs/superpowers/specs/2026-05-11-ref-panel-stale-class2-warning.md`
   to keep the liftover change set focused on liftover.
