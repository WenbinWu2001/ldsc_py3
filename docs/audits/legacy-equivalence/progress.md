# Progress Log — Legacy↔Refactored Equivalence Audit

## Session: 2026-06-14

### Phase 0: Setup
- **Status:** complete
- Actions taken:
  - Confirmed refactored repo on branch `restructure`; legacy at ldsc_py2_Bulik.
  - Found existing root planning files belong to a prior "Threshold Comparison Audit" — NOT clobbered.
  - Confirmed workspace root sits inside parent git repo `.../LDSC` → placed audit files in dedicated
    subdir `docs/audits/legacy-equivalence/` under the refactored repo instead.
  - Read design_map.md + restructured CLAUDE.md → refactored CLI = `ldsc` with 8 subcommands.
  - Scaffolded task_plan.md, findings.md, progress.md.
- Files created:
  - docs/audits/legacy-equivalence/{task_plan,findings,progress}.md

### Phase 1: Reference ingestion — legacy
- **Status:** complete
- Read legacy ldsc.py (49 flags + dispatch 582-661), make_annot.py (7), munge_sumstats.py (29). All -> findings.

### Phase 2: Reference ingestion — refactored
- **Status:** complete (wiki examples deferred to Phase 5)
- Read cli.py (8 subcommands), io-argument-inventory.md (full), and verified every flag against
  actual add_argument source in annotation_builder/ldscore_calculator/ref_panel_builder/r2_query/
  sumstats_munger(+_kernel)/regression_runner. Captured refactored inventory -> findings.

### Phase 3: Command coverage mapping
- **Status:** complete
- Mapping table in findings. Verified absence of cts-bin/per-allele/pq-exp/h2-cts/ref-ld-chr-cts/
  gene-set/merge-alleles in refactored via grep. read_cts is dead code; cell-type regime exists.

### Phase 4: Flag-level coverage + behavior
- **Status:** complete
- 7 per-command flag tables (A-G) + 7 behavioral discrepancies in findings. Grounded exclude-regions
  default (ldscore_calculator.py:930), common-maf-min 0.05, chunksize 5e6->1e6, print-* NONE.

### Phase 5: Wiki-example reproduction
- **Status:** complete
- Read 6 tutorial PDFs (H&GC, LD Score Estimation, Partitioned, Cell-type, Continuous Annotations; munge
  examples across them). All commands extracted + translated in findings. Confirmed CRITICAL legacy-LD-score
  format blocker (load_ldscore_from_dir:2465). Refined cts-bin (HIGH) vs continuous-annot-files (supported).
  Remaining wiki pages (Home/FAQ/LD File Formats/Summary Stats Format/Tests/What Data) are conceptual; noted.

### Phase 6: Severity assessment + report assembly
- **Status:** complete
- Reconciled severities (cts-bin CRITICAL->HIGH; 2 true CRITICALs = format blocker + exclude-regions default).
- Delivered report: docs/audits/legacy-equivalence/legacy-equivalence-audit.md (373 lines).
- Gap counts: CRITICAL 2, HIGH 7, MEDIUM ~12, LOW ~14. Flag coverage 65/85 (~76%). All 7 capabilities >= partial.

## 5-Question Reboot Check (final)
| Question | Answer |
|----------|--------|
| Where am I? | DONE — all 6 phases complete |
| Where am I going? | Nothing pending; report delivered |
| What's the goal? | Interface/behavior equivalence audit + severity report (NOT numerical) |
| What have I learned? | findings.md (full per-flag + per-tutorial detail, all source-grounded) |
| What have I done? | Read both codebases + 6 wiki tutorials; built coverage map, 7 flag tables, behavioral discrepancies, wiki translations; assembled report |

### Phase 6: Severity assessment + report assembly
- **Status:** pending

## Error Log
| Timestamp | Error | Attempt | Resolution |
|-----------|-------|---------|------------|
|           |       | 1       |            |

## 5-Question Reboot Check
| Question | Answer |
|----------|--------|
| Where am I? | Phase 1 (legacy reference ingestion) |
| Where am I going? | Phases 2-6: wiki+refactored ingestion, coverage map, flag diff, wiki translation, report |
| What's the goal? | Source-grounded interface/behavior equivalence audit + severity-rated report |
| What have I learned? | See findings.md |
| What have I done? | Phase 0 setup complete; see above |
