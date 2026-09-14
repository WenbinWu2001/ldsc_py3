# Staged input preflight and progress audit

Last updated on: 2026-09-14

## Verified starting state

Local `ldsc_py3_restructured`, branch `restructure`, began clean at `547bf93`. This includes `67c2e57` (PLINK chromosome resolution), `6454dd8` (output retries), and the shared worker policy. Historical Longleaf revisions and job timing are background supplied by the user; this work does not inspect remote data, run scheduler operations, or deploy changes remotely.

## Confirmed order and late checks at the starting revision

| Workflow | Observed order | Consequence |
| --- | --- | --- |
| Direct LD scores | options/build inference → output gate → gene/restriction loading → `prepare_annotation_sources` → `validate_direct_scope` for queries → calculation → writer | Reference integrity can fail after annotation staging; baseline-only runs lack the query scope gate. |
| Gene-index construction | config/output transaction → catalog → restriction/map → PLINK inspection → annotation preparation → workers → publication/reload | PLINK ordering is already fixed; independent source families still fail separately, and annotation header defects in late files are discovered during staging. |
| Indexed LD scores | index root/catalog → each chromosome's required paths then its full operator validation → gene lists → assembly/writer | Late missing component paths follow earlier operator loads; the direct Python API defers output collisions to the writer. |
| R2 construction | PLINK inspection → source-build/state resolution → output gate/log → per-chromosome construction/publication | PLINK errors aggregate, but auxiliary paths are inspected separately; keep-file resolution occurs inside chromosome work. |
| Standalone annotate | output gate → catalog/gene lists/build inference → baseline preparation → BED preparation/projection → writer | Missing BED inputs and late baseline headers can follow expensive preparation. |
| Annotation preparation API | baseline/query resolution → preparation → catalog/gene/BED work | Query declarations can fail after baseline preparation. |

Sources: `run_ldscore_from_args`, `prepare_direct_annotations`, `inspect_direct_inputs`, `run_build_gene_ldscore_index_from_args`, `_run_gene_ldscore_index_build`, `_load_gene_ldscore_index`, `ReferencePanelBuilder._run`, `run_annotate`, and `build_annotation_shards` in `src/ldsc/`.

## Gate boundaries

1. **Options and output ownership:** retain configuration validation and the existing artifact-family/transaction checks. Output authorization is separate from evidence that inputs are valid. Worker options use `_parallelism._validate_threads`.
2. **Input declarations and headers:** resolve every independent requested object; inspect required companions and bounded headers before content scans. Collect repair records in the existing six-column input-issue format. Concrete paths and annotation widths are retained for preparation.
3. **Content validation:** inspect BIM/FAM and BED dimensions; validate immutable index payloads; stream annotations into bounded metadata/value scratch; validate alignment and identities using those scratch results. Full source contents cannot be inferred from names or headers.
4. **Derived coverage and support:** compare validated chromosome sets, resolved catalog selections, identity-cleaned reference intersections, and post-filter support. These checks require preceding content; moving them to path discovery would change scientific policy.
5. **Computation and publication:** preserve numerical routines, worker dispatch, output schemas, and overwrite/failure-marker lifecycles. Publication logs do not imply earlier private staging was a complete result.

Progress uses the existing logging routes and INFO level: phase start/completion/failure plus time-throttled chunk progress, current source/chromosome, completed counts, totals when known, and elapsed seconds. It does not install another console handler, estimate completion dates, hash inputs, or retain whole annotation matrices.

## Remaining command audit

The following inspection covers the other CLI owners; output preflight alone is never classified as input validation.

| Commands | Existing input boundary / content-dependent limitation |
| --- | --- |
| `h2`, `partitioned-h2`, `rg` | Shared runner configuration, summary-statistic loading, LD-score manifest/schema validation, then identity alignment and fitting. Trait contents and cross-trait identity require loading; independent trait failures can aggregate before fitting. |
| `munge-sumstats` | Build/liftover option validation and output gate precede raw-input inference, restriction and QC. Header profiles and source-build inference depend on sampled/content data; stdin-like or infer-only behavior must retain its no-write policy. |
| `convert-ldsc2-ldscores` | Explicit legacy family discovery and conversion diagnostics own chromosome coverage, companion counts, schemas, and cross-family equivalence. Duplicate compressed/plain representations require content comparison under the existing contract. |
| `quantile-h2` | Output gate precedes fitted-model loading and annotation/statistic preparation with alignment diagnostics. Requested targets and provenance depend on the saved model; bin/support checks require selected values. |
| `query-r2` | Output ownership and pair/panel declarations now precede pair loading. Stdin is consumed only at its content stage. Pair membership and values require panel content. |
| `convert-h2-scale`, `plot` | Single saved-result workflows derive output locations from validated result metadata and have existing artifact-family checks. Required fitted-result fields and plot compatibility depend on the saved result; no chromosome preparation loop is involved. |

## Implemented order and limits

- Direct and indexed LD scores retain the shared `_validate_threads` option gate and existing worker dispatch. Direct input/header checks precede gene resolution, reference/control content validation, annotation staging, derived scope, computation, and publication. Reference scope validation is reused after annotation preparation.
- Index construction checks baseline headers (including empty inputs, declared first-row chromosome mismatches, and differing declared shard columns), auxiliary declarations, and all selected PLINK companions before the existing catalog/content stages. Its reference/control content gate aggregates regression-list, map, and PLINK defects before annotation scans. The loader checks all chromosome component paths, small identity metadata, and readable Parquet footers before the first operator load.
- R2 construction aggregates companion, auxiliary path, and genetic-map header failures first. PLINK content determines output chromosome scope; source-build inference determines emitted build paths. The existing output-family gate precedes shared build-state preparation. Early repair tables retain the previous PLINK-diagnostic lifecycle; canonical file logging starts only when that output scope is known.
- Standalone annotate checks all baseline headers and BED/gene/catalog paths before baseline preparation. Its gene Gate A stops before baseline scans. The preparation-only Python API also resolves gene/interval inputs before baseline staging; ordinary BED skip policies remain unchanged.
- Regression, quantile, munging, and R2-query declarations now precede their input loading or numerical stages after their own output authorization. `rg` aggregates independent trait content failures before fitting. Legacy conversion, plotting, and scale conversion retain their specialized result/family validators; their audit did not justify another generic chromosome staging gate.

Each successful input gate clears its own earlier repair table, so a repaired retry cannot retain stale input-failure diagnostics; the real R2 retry test covers this without overwrite. Scientific files and failure markers remain under their existing workflow owners.

No new whole-file content hashes were added. Header reads are bounded; annotation widths, concrete BED paths, PLINK discovery members, and validated reference mappings are reused. Existing annotation scans still stage each source once and reuse metadata/value streams for alignment and identity passes. Whole-chromosome numerical data are not cached by the new gates.

Content-dependent limitations remain explicit: global identity/alignment and ordinary-glob chromosome coverage need content; genotype-filtered support needs prepared genotypes; map ordering and interpolation support need map contents and resolved reference coordinates; immutable CSR dimensions/values and existing integrity identities need operator payloads. Baseline-only R2 validation does not add a whole pair scan: pair integrity remains checked by the consuming reader, while query-scope runs retain their existing streamed pair validation before staging. Ignored genetic-map and opposite-direction-chain policies are unchanged. Progress is time-throttled at returned chunks/objects and cannot interrupt a blocking library call. Annotation scanning validates while staging, and R2 construction writes pair chunks as they are computed; combined phase labels describe these existing streamed boundaries without adding a second scan or a complete intermediate pair table.

## Verification record

Completed locally in the editable `ldsc3-dev` environment:

| Check | Result |
| --- | --- |
| `python -m pytest tests/test_staged_input_preflight.py tests/test_annotation_preparation.py -q` | 40 passed |
| `python -m pytest -q` | 1,799 passed, 132 subtests passed, one skip; 311 warnings; 138.38 seconds |
| `python -m unittest discover -s tests -p 'test*.py' -v` | 1,000 tests, OK with one skip |
| `ldsc --help`, `python -m ldsc --help` | Both exited successfully |
| `git diff --check` | Clean |

The focused tests intercept expensive annotation scans, index operator loads, and chromosome computation; verify aggregated errors across late chromosomes, headers, controls, and malformed manifests; check structured progress records; and prove reference discovery is reused. Real PLINK fixtures retain independently expected self-window LD scores, while the full suite includes direct/indexed numerical agreement, both identity modes, artifact validation, worker policy, and overwrite/failure-marker contracts. The R2 retry test repairs a missing companion, successfully builds both chromosomes, and verifies that the earlier input repair table is removed.

Pytest and unittest were run sequentially for final verification. Historical remote revisions and job observations were not used as claims about current deployment; no remote repository or scheduler operation was performed.
