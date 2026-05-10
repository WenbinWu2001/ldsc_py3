# Duplicate-Position Policy Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `--duplicate-position-policy {error,drop-all}` to `build-ref-panel` so that SNPs sharing a `CHR:POS` key in any emitted build are either rejected loudly or dropped with a provenance sidecar.

**2026-05-10 update:** Liftover harmonization changed the final contract from
this original implementation plan. `duplicate_position_policy` now defaults to
`drop-all`, coordinate duplicate handling applies only in `chr_pos` mode,
source-only `rsid` builds log that the policy is not applicable, and matching
reference-panel chain liftover is rejected in `rsid` mode. The `error` policy
remains available explicitly.

**Architecture:** A new helper `_resolve_unique_snp_set()` in `ref_panel_builder.py` runs two detection passes (source-build duplicates, then target-build collisions) on the restriction-filtered SNP set and returns a cleaned `keep_snps` array plus a provenance DataFrame. `_build_chromosome()` calls this helper after liftover and before the per-build emit loop. Cross-build consistency is enforced because the same cleaned set feeds every emitted build.

**Tech Stack:** Python 3.11, pandas, numpy, unittest (existing test harness)

---

## File Map

| File | Change |
|---|---|
| `src/ldsc/config.py` | Add `duplicate_position_policy: str = "drop-all"` field + validation to `ReferencePanelBuildConfig` |
| `src/ldsc/ref_panel_builder.py` | Add CLI flag, extend `config_from_args()`, add `_resolve_unique_snp_set()`, add `_write_dropped_sidecar()`, wire into `_build_chromosome()` |
| `tests/test_ref_panel_builder.py` | Add `DuplicatePositionPolicyTest` unit-test class and two integration test methods |

---

## Task 1: Add `duplicate_position_policy` to `ReferencePanelBuildConfig`

**Files:**
- Modify: `src/ldsc/config.py:444-487`
- Test: `tests/test_ref_panel_builder.py`

- [ ] **Step 1: Write the failing test**

Add a new test class after the existing `ReferencePanelBuildConfigTest` class in `tests/test_ref_panel_builder.py`:

```python
class ReferencePanelBuildConfigDuplicatePolicyTest(unittest.TestCase):
    def _base_config(self, **kwargs):
        return ReferencePanelBuildConfig(
            plink_prefix="plink/panel.@",
            source_genome_build="hg19",
            output_dir="out",
            ld_wind_kb=1.0,
            **kwargs,
        )

    def test_default_policy_is_drop_all(self):
        config = self._base_config()
        self.assertEqual(config.duplicate_position_policy, "drop-all")

    def test_drop_all_policy_is_accepted(self):
        config = self._base_config(duplicate_position_policy="drop-all")
        self.assertEqual(config.duplicate_position_policy, "drop-all")

    def test_invalid_policy_raises(self):
        with self.assertRaisesRegex(ValueError, "duplicate_position_policy"):
            self._base_config(duplicate_position_policy="keep-one")
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/test_ref_panel_builder.py::ReferencePanelBuildConfigDuplicatePolicyTest -v
```

Expected: `FAILED` — `ReferencePanelBuildConfig` has no `duplicate_position_policy` field.

- [ ] **Step 3: Add the field to `ReferencePanelBuildConfig` in `src/ldsc/config.py`**

After `overwrite: bool = False` (currently the last field, line ~445), add:

```python
    duplicate_position_policy: str = "drop-all"
```

In `__post_init__` (after the `snp_batch_size` validation block), add:

```python
        if self.duplicate_position_policy not in {"error", "drop-all"}:
            raise ValueError(
                f"duplicate_position_policy must be 'error' or 'drop-all', "
                f"got {self.duplicate_position_policy!r}."
            )
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/test_ref_panel_builder.py::ReferencePanelBuildConfigDuplicatePolicyTest -v
```

Expected: 3 tests PASS.

- [ ] **Step 5: Run the full suite to check for regressions**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/test_ref_panel_builder.py -v
```

Expected: all existing tests still PASS.

- [ ] **Step 6: Commit**

```bash
cd ldsc_py3_Jerry && git add src/ldsc/config.py tests/test_ref_panel_builder.py && git commit -m "feat: add duplicate_position_policy field to ReferencePanelBuildConfig"
```

---

## Task 2: Extend CLI parser and `config_from_args()`

**Files:**
- Modify: `src/ldsc/ref_panel_builder.py` — `build_parser()` and `config_from_args()`
- Test: `tests/test_ref_panel_builder.py`

- [ ] **Step 1: Write the failing tests**

Add these methods to `ReferencePanelBuildConfigFromArgsTest` in `tests/test_ref_panel_builder.py`:

```python
    def test_build_parser_defaults_duplicate_position_policy_to_drop_all(self):
        parser = ref_panel_builder.build_parser()
        self.assertEqual(parser.get_default("duplicate_position_policy"), "drop-all")

    def test_build_parser_accepts_drop_all_policy(self):
        parser = ref_panel_builder.build_parser()
        args = parser.parse_args([
            "--plink-prefix", "plink/panel.@",
            "--output-dir", "out",
            "--ld-wind-kb", "1",
            "--duplicate-position-policy", "drop-all",
        ])
        self.assertEqual(args.duplicate_position_policy, "drop-all")

    def test_config_from_args_passes_policy_to_config(self):
        parser = ref_panel_builder.build_parser()
        args = parser.parse_args([
            "--plink-prefix", "plink/panel.@",
            "--output-dir", "out",
            "--ld-wind-kb", "1",
            "--duplicate-position-policy", "drop-all",
        ])
        build_config, _ = ref_panel_builder.config_from_args(args)
        self.assertEqual(build_config.duplicate_position_policy, "drop-all")
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/test_ref_panel_builder.py::ReferencePanelBuildConfigFromArgsTest::test_build_parser_defaults_duplicate_position_policy_to_drop_all tests/test_ref_panel_builder.py::ReferencePanelBuildConfigFromArgsTest::test_build_parser_accepts_drop_all_policy tests/test_ref_panel_builder.py::ReferencePanelBuildConfigFromArgsTest::test_config_from_args_passes_policy_to_config -v
```

Expected: `FAILED` — `--duplicate-position-policy` not recognised.

- [ ] **Step 3: Add the flag to `build_parser()` in `src/ldsc/ref_panel_builder.py`**

After the `--snp-batch-size` argument, add:

```python
    parser.add_argument(
        "--duplicate-position-policy",
        default="drop-all",
        choices=("error", "drop-all"),
        help=(
            "How to handle SNPs that share a CHR:POS key in any emitted build. "
            "'drop-all' drops every SNP in each colliding cluster and writes a "
            "provenance sidecar to {output_dir}/dropped_snps/chr{chrom}_dropped.tsv.gz (default). "
            "'error' aborts and reports all duplicate clusters."
        ),
    )
```

- [ ] **Step 4: Pass the policy in `config_from_args()` in `src/ldsc/ref_panel_builder.py`**

In `config_from_args()`, inside the `ReferencePanelBuildConfig(...)` constructor call, add:

```python
        duplicate_position_policy=args.duplicate_position_policy,
```

- [ ] **Step 5: Run tests to verify they pass**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/test_ref_panel_builder.py::ReferencePanelBuildConfigFromArgsTest -v
```

Expected: all tests in `ReferencePanelBuildConfigFromArgsTest` PASS.

- [ ] **Step 6: Commit**

```bash
cd ldsc_py3_Jerry && git add src/ldsc/ref_panel_builder.py tests/test_ref_panel_builder.py && git commit -m "feat: add --duplicate-position-policy CLI flag to build-ref-panel"
```

---

## Task 3: Implement `_resolve_unique_snp_set()`

**Files:**
- Modify: `src/ldsc/ref_panel_builder.py` — add new module-private function
- Test: `tests/test_ref_panel_builder.py` — add `DuplicatePositionPolicyTest` class

- [ ] **Step 1: Write the failing tests**

Add a new test class `DuplicatePositionPolicyTest` in `tests/test_ref_panel_builder.py`:

```python
class DuplicatePositionPolicyTest(unittest.TestCase):
    """Unit tests for _resolve_unique_snp_set()."""

    def _make_chrom_df(self, snp_ids, bps):
        """Build a minimal .bim-style DataFrame indexed by PLINK row indices."""
        return pd.DataFrame(
            {"SNP": snp_ids, "BP": bps},
            index=list(range(len(snp_ids))),
        )

    def test_no_duplicates_returns_keep_snps_unchanged(self):
        chrom_df = self._make_chrom_df(["rs1", "rs2", "rs3"], [100, 200, 300])
        keep = np.array([0, 1, 2])
        hg19 = {0: 100, 1: 200, 2: 300}
        hg38 = {0: 1000, 1: 2000, 2: 3000}

        cleaned, dropped = ref_panel_builder._resolve_unique_snp_set(
            "1", chrom_df, keep, hg19, hg38, "error"
        )

        np.testing.assert_array_equal(cleaned, keep)
        self.assertTrue(dropped.empty)

    def test_source_duplicate_error_policy_raises(self):
        chrom_df = self._make_chrom_df(["rs1", "rs2", "rs3"], [100, 100, 300])
        keep = np.array([0, 1, 2])

        with self.assertRaisesRegex(ValueError, "source-build duplicate"):
            ref_panel_builder._resolve_unique_snp_set(
                "1", chrom_df, keep, {0: 100, 1: 100, 2: 300}, {}, "error"
            )

    def test_source_duplicate_drop_all_removes_cluster(self):
        chrom_df = self._make_chrom_df(["rs1", "rs2", "rs3"], [100, 100, 300])
        keep = np.array([0, 1, 2])

        cleaned, dropped = ref_panel_builder._resolve_unique_snp_set(
            "1", chrom_df, keep, {0: 100, 1: 100, 2: 300}, {}, "drop-all"
        )

        np.testing.assert_array_equal(cleaned, np.array([2]))
        self.assertEqual(len(dropped), 2)
        self.assertTrue((dropped["reason"] == "source_duplicate").all())
        self.assertTrue(dropped["target_pos"].isna().all())
        self.assertSetEqual(set(dropped["SNP"]), {"rs1", "rs2"})

    def test_target_collision_error_policy_raises(self):
        chrom_df = self._make_chrom_df(["rs1", "rs2", "rs3"], [100, 200, 300])
        keep = np.array([0, 1, 2])
        hg38 = {0: 5000, 1: 5000, 2: 6000}  # rs1 and rs2 collide in hg38

        with self.assertRaisesRegex(ValueError, "target-build collision"):
            ref_panel_builder._resolve_unique_snp_set(
                "1", chrom_df, keep, {0: 100, 1: 200, 2: 300}, hg38, "error"
            )

    def test_target_collision_drop_all_removes_cluster_from_all_builds(self):
        chrom_df = self._make_chrom_df(["rs1", "rs2", "rs3"], [100, 200, 300])
        keep = np.array([0, 1, 2])
        hg19 = {0: 100, 1: 200, 2: 300}
        hg38 = {0: 5000, 1: 5000, 2: 6000}  # rs1 and rs2 collide in hg38

        cleaned, dropped = ref_panel_builder._resolve_unique_snp_set(
            "1", chrom_df, keep, hg19, hg38, "drop-all"
        )

        np.testing.assert_array_equal(cleaned, np.array([2]))
        self.assertEqual(len(dropped), 2)
        self.assertTrue((dropped["reason"] == "target_collision").all())
        self.assertSetEqual(set(dropped["SNP"]), {"rs1", "rs2"})
        # target_pos must be the colliding hg38 position, not NA
        self.assertFalse(dropped["target_pos"].isna().any())

    def test_target_collision_in_hg19_only_drops_both_snps(self):
        chrom_df = self._make_chrom_df(["rs1", "rs2", "rs3"], [100, 200, 300])
        keep = np.array([0, 1, 2])
        hg19 = {0: 5000, 1: 5000, 2: 6000}  # rs1 and rs2 collide in hg19
        hg38 = {}  # source-only build

        cleaned, dropped = ref_panel_builder._resolve_unique_snp_set(
            "1", chrom_df, keep, hg19, hg38, "drop-all"
        )

        np.testing.assert_array_equal(cleaned, np.array([2]))
        self.assertEqual(len(dropped), 2)
        self.assertTrue((dropped["reason"] == "target_collision").all())

    def test_collision_in_both_builds_deduplicates_provenance_rows(self):
        # rs1 and rs2 collide in both hg19 and hg38 — should appear once each
        chrom_df = self._make_chrom_df(["rs1", "rs2", "rs3"], [100, 200, 300])
        keep = np.array([0, 1, 2])
        hg19 = {0: 5000, 1: 5000, 2: 6000}
        hg38 = {0: 7000, 1: 7000, 2: 8000}

        cleaned, dropped = ref_panel_builder._resolve_unique_snp_set(
            "1", chrom_df, keep, hg19, hg38, "drop-all"
        )

        np.testing.assert_array_equal(cleaned, np.array([2]))
        self.assertEqual(len(dropped), 2)  # one row per SNP, not per build
        self.assertSetEqual(set(dropped["SNP"]), {"rs1", "rs2"})

    def test_all_snps_dropped_returns_empty_keep(self):
        chrom_df = self._make_chrom_df(["rs1", "rs2"], [100, 100])
        keep = np.array([0, 1])

        cleaned, dropped = ref_panel_builder._resolve_unique_snp_set(
            "1", chrom_df, keep, {0: 100, 1: 100}, {}, "drop-all"
        )

        self.assertEqual(len(cleaned), 0)
        self.assertEqual(len(dropped), 2)

    def test_provenance_columns_are_correct(self):
        chrom_df = self._make_chrom_df(["rs1", "rs2", "rs3"], [100, 100, 300])
        keep = np.array([0, 1, 2])

        _, dropped = ref_panel_builder._resolve_unique_snp_set(
            "1", chrom_df, keep, {0: 100, 1: 100, 2: 300}, {}, "drop-all"
        )

        self.assertListEqual(list(dropped.columns), ["CHR", "SNP", "source_pos", "target_pos", "reason"])
        self.assertEqual(dropped["CHR"].iloc[0], "1")
        self.assertEqual(dropped["source_pos"].dtype, int)
        # target_pos must be nullable Int64
        self.assertEqual(dropped["target_pos"].dtype.name, "Int64")
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/test_ref_panel_builder.py::DuplicatePositionPolicyTest -v
```

Expected: `AttributeError: module ... has no attribute '_resolve_unique_snp_set'`.

- [ ] **Step 3: Implement `_resolve_unique_snp_set()` in `src/ldsc/ref_panel_builder.py`**

Add the following function just before `build_parser()` (around line 809):

```python
def _resolve_unique_snp_set(
    chrom: str,
    chrom_df: pd.DataFrame,
    keep_snps: np.ndarray,
    hg19_lookup: dict[int, int],
    hg38_lookup: dict[int, int],
    policy: str,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Detect and resolve SNPs sharing a CHR:POS key in any emitted build.

    Pass 1 detects source-build duplicates (same BP in the .bim). Pass 2
    detects target-build collisions (same translated position in either lookup)
    against the SNPs surviving pass 1. Both passes apply the same policy.
    Each dropped SNP contributes one provenance row; a two-SNP cluster produces
    two rows.

    Parameters
    ----------
    chrom : str
        Chromosome label, used in error messages and provenance.
    chrom_df : pd.DataFrame
        Source .bim rows for this chromosome, indexed by PLINK row index.
        Must have ``SNP`` and ``BP`` columns.
    keep_snps : np.ndarray
        Restriction-filtered PLINK row indices to check.
    hg19_lookup : dict[int, int]
        Map from PLINK row index to hg19 position. Empty for hg38-source-only builds.
    hg38_lookup : dict[int, int]
        Map from PLINK row index to hg38 position. Empty for hg19-source-only builds.
    policy : str
        ``"error"`` raises ``ValueError`` listing all duplicate clusters.
        ``"drop-all"`` removes all SNPs in each cluster and returns provenance rows.

    Returns
    -------
    cleaned_keep_snps : np.ndarray
        Subset of ``keep_snps`` with all colliders removed (same dtype).
    provenance : pd.DataFrame
        One row per dropped SNP with columns
        ``CHR, SNP, source_pos, target_pos, reason``. Empty if nothing dropped.
    """
    keep_list = keep_snps.tolist()
    drop_set: set[int] = set()
    prov_rows: list[dict] = []

    # Pass 1: source-build duplicates
    source_pos = chrom_df.loc[keep_list, "BP"].to_numpy(dtype=int)
    dup_mask = pd.Series(source_pos).duplicated(keep=False).to_numpy()
    if dup_mask.any():
        if policy == "error":
            dup_positions = sorted(set(source_pos[dup_mask]))
            lines = []
            for pos in dup_positions:
                for idx, p in zip(keep_list, source_pos):
                    if p == pos:
                        lines.append(f"    CHR={chrom} source_pos={pos} SNP={chrom_df.loc[idx, 'SNP']}")
            raise ValueError(
                f"Source-build duplicate CHR:POS on chromosome {chrom}. "
                "Use --duplicate-position-policy=drop-all to drop colliding clusters.\n"
                + "\n".join(lines)
            )
        for i, idx in enumerate(keep_list):
            if dup_mask[i]:
                drop_set.add(int(idx))
                prov_rows.append({
                    "CHR": chrom,
                    "SNP": chrom_df.loc[idx, "SNP"],
                    "source_pos": int(source_pos[i]),
                    "target_pos": pd.NA,
                    "reason": "source_duplicate",
                })

    # Survivors of pass 1
    surviving = [idx for idx in keep_list if idx not in drop_set]

    # Pass 2: target-build collisions — check both non-empty lookups against surviving set
    for lookup in (hg19_lookup, hg38_lookup):
        if not lookup:
            continue
        tgt_pos = [lookup[idx] for idx in surviving]
        dup_mask2 = pd.Series(tgt_pos).duplicated(keep=False).to_numpy()
        if not dup_mask2.any():
            continue
        if policy == "error":
            dup_positions = sorted(set(np.array(tgt_pos)[dup_mask2]))
            lines = []
            for pos in dup_positions:
                for idx, p in zip(surviving, tgt_pos):
                    if p == pos:
                        lines.append(
                            f"    CHR={chrom} source_pos={chrom_df.loc[idx, 'BP']} "
                            f"target_pos={pos} SNP={chrom_df.loc[idx, 'SNP']}"
                        )
            raise ValueError(
                f"Target-build collision CHR:POS on chromosome {chrom}. "
                "Use --duplicate-position-policy=drop-all to drop colliding clusters.\n"
                + "\n".join(lines)
            )
        for i, idx in enumerate(surviving):
            if dup_mask2[i] and idx not in drop_set:
                drop_set.add(int(idx))
                prov_rows.append({
                    "CHR": chrom,
                    "SNP": chrom_df.loc[idx, "SNP"],
                    "source_pos": int(chrom_df.loc[idx, "BP"]),
                    "target_pos": int(tgt_pos[i]),
                    "reason": "target_collision",
                })

    cleaned = np.asarray([idx for idx in keep_snps if idx not in drop_set], dtype=keep_snps.dtype)

    if prov_rows:
        dropped_df = pd.DataFrame(prov_rows, columns=["CHR", "SNP", "source_pos", "target_pos", "reason"])
        dropped_df["target_pos"] = dropped_df["target_pos"].astype("Int64")
    else:
        dropped_df = pd.DataFrame(columns=["CHR", "SNP", "source_pos", "target_pos", "reason"])

    return cleaned, dropped_df
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/test_ref_panel_builder.py::DuplicatePositionPolicyTest -v
```

Expected: all 9 tests PASS.

- [ ] **Step 5: Run the full suite to check for regressions**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/test_ref_panel_builder.py -v
```

Expected: all tests PASS.

- [ ] **Step 6: Commit**

```bash
cd ldsc_py3_Jerry && git add src/ldsc/ref_panel_builder.py tests/test_ref_panel_builder.py && git commit -m "feat: implement _resolve_unique_snp_set() for duplicate-position detection"
```

---

## Task 4: Implement `_write_dropped_sidecar()` and wire into `_build_chromosome()`

**Files:**
- Modify: `src/ldsc/ref_panel_builder.py` — add `_write_dropped_sidecar()`, update `_build_chromosome()`
- Test: `tests/test_ref_panel_builder.py` — add integration tests

- [ ] **Step 1: Write the failing integration tests**

Add these two methods to `ReferencePanelBuilderRunTest` in `tests/test_ref_panel_builder.py`.

The tests use `_write_plink_prefix_rows` (already in the class) to create a `.bim` with
duplicate positions and mock `_resolve_mappable_snp_positions` so the builder never
reaches the `.bed` read.  The "drop-all" test additionally mocks the per-build emit
helpers so no real parquet is written.

```python
    def test_error_policy_raises_on_source_build_duplicates(self):
        """run() raises ValueError when .bim has duplicate CHR:POS and policy=error."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_plink_prefix_rows(tmpdir, "panel.1", [
                ("1", "rs1", 100),
                ("1", "rs2", 100),  # duplicate position
                ("1", "rs3", 200),
            ])
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
                duplicate_position_policy="error",
            )
            builder = ref_panel_builder.ReferencePanelBuilder(
                global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg19")
            )

            def fake_resolve(self_arg, *, build_state, chrom, source_build, chrom_df, keep_snps):
                keep = np.asarray(keep_snps, dtype=int)
                hg19 = {int(idx): int(chrom_df.loc[idx, "BP"]) for idx in keep}
                return keep, hg19, {}

            with mock.patch.object(
                ref_panel_builder.ReferencePanelBuilder,
                "_resolve_mappable_snp_positions",
                side_effect=fake_resolve,
            ):
                with self.assertRaisesRegex(ValueError, "source-build duplicate"):
                    builder.run(config)

    def test_drop_all_policy_writes_sidecar_under_dropped_snps(self):
        """run() writes chr{chrom}_dropped.tsv.gz under dropped_snps when duplicates are dropped."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            self._write_plink_prefix_rows(tmpdir, "panel.1", [
                ("1", "rs1", 100),
                ("1", "rs2", 100),  # duplicate of rs1
                ("1", "rs3", 200),
            ])
            config = ReferencePanelBuildConfig(
                plink_prefix=tmpdir / "panel.@",
                source_genome_build="hg19",
                output_dir=tmpdir / "out",
                ld_wind_kb=1.0,
                duplicate_position_policy="drop-all",
            )
            builder = ref_panel_builder.ReferencePanelBuilder(
                global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="hg19")
            )

            def fake_resolve(self_arg, *, build_state, chrom, source_build, chrom_df, keep_snps):
                keep = np.asarray(keep_snps, dtype=int)
                hg19 = {int(idx): int(chrom_df.loc[idx, "BP"]) for idx in keep}
                return keep, hg19, {}

            # Mock emit helpers so no real parquet/bed read occurs
            with mock.patch.object(
                ref_panel_builder.ReferencePanelBuilder,
                "_resolve_mappable_snp_positions",
                side_effect=fake_resolve,
            ), mock.patch(
                "ldsc.ref_panel_builder.kernel_builder.write_r2_parquet"
            ), mock.patch(
                "ldsc.ref_panel_builder.kernel_builder.write_runtime_metadata_sidecar"
            ), mock.patch(
                "ldsc.ref_panel_builder.kernel_ldscore.PlinkBEDFile"
            ) as mock_bed:
                mock_bed.return_value.kept_snps = [2]
                mock_bed.return_value.maf = np.array([0.3])
                mock_bed.return_value.m = 1
                mock_bed.return_value.n = 1
                mock_bed.return_value.nextSNPs = lambda: iter([np.array([0.0])])

                with self.assertLogs("LDSC.ref_panel_builder", level="WARNING") as log_ctx:
                    builder.run(config)

            sidecar = tmpdir / "out" / "dropped_snps" / "chr1_dropped.tsv.gz"
            self.assertTrue(sidecar.exists(), "sidecar must be written under dropped_snps, not under a build subdir")

            import gzip
            with gzip.open(sidecar, "rt") as fh:
                dropped_df = pd.read_csv(fh, sep="\t")
            self.assertEqual(len(dropped_df), 2)
            self.assertTrue((dropped_df["reason"] == "source_duplicate").all())
            self.assertTrue(
                any("chr1_dropped.tsv.gz" in line for line in log_ctx.output),
                "WARNING log must reference the sidecar path",
            )
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/test_ref_panel_builder.py::ReferencePanelBuilderRunTest::test_error_policy_raises_on_source_build_duplicates tests/test_ref_panel_builder.py::ReferencePanelBuilderRunTest::test_drop_all_policy_writes_sidecar_at_panel_root -v
```

Expected: `FAILED` — `_resolve_unique_snp_set` is not called yet from `_build_chromosome`.

- [ ] **Step 3: Add `_write_dropped_sidecar()` to `src/ldsc/ref_panel_builder.py`**

Add just before `build_parser()`:

```python
def _write_dropped_sidecar(dropped_df: pd.DataFrame, path: Path, chrom: str) -> None:
    """Write a provenance sidecar for dropped duplicate-position SNPs."""
    path.parent.mkdir(parents=True, exist_ok=True)
    dropped_df.to_csv(path, sep="\t", index=False, compression="gzip")
    n = len(dropped_df)
    n_source = int((dropped_df["reason"] == "source_duplicate").sum())
    n_target = int((dropped_df["reason"] == "target_collision").sum())
    LOGGER.warning(
        f"Dropped {n} SNPs on chromosome {chrom} due to duplicate positions "
        f"({n_source} source-build duplicates, {n_target} target-build collisions). "
        f"Provenance written to '{path}'."
    )
```

- [ ] **Step 4: Wire `_resolve_unique_snp_set()` into `_build_chromosome()`**

In `ReferencePanelBuilder._build_chromosome()`, after the call to `self._resolve_mappable_snp_positions(...)` and before the `keep_indivs` block (around line 393), add:

```python
        keep_snps, dropped_df = _resolve_unique_snp_set(
            chrom=chrom,
            chrom_df=chrom_df,
            keep_snps=keep_snps,
            hg19_lookup=hg19_lookup,
            hg38_lookup=hg38_lookup,
            policy=config.duplicate_position_policy,
        )
        if not dropped_df.empty:
            sidecar_path = Path(config.output_dir) / "dropped_snps" / f"chr{chrom}_dropped.tsv.gz"
            _write_dropped_sidecar(dropped_df, sidecar_path, chrom)
        if len(keep_snps) == 0:
            LOGGER.info(f"Skipping chromosome {chrom}: no SNPs remain after duplicate-position filtering.")
            return None
```

- [ ] **Step 5: Run integration tests to verify they pass**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/test_ref_panel_builder.py::ReferencePanelBuilderRunTest::test_error_policy_raises_on_source_build_duplicates tests/test_ref_panel_builder.py::ReferencePanelBuilderRunTest::test_drop_all_policy_writes_sidecar_at_panel_root -v
```

Expected: both tests PASS.

- [ ] **Step 6: Run the full suite to check for regressions**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/test_ref_panel_builder.py -v
```

Expected: all tests PASS.

- [ ] **Step 7: Commit**

```bash
cd ldsc_py3_Jerry && git add src/ldsc/ref_panel_builder.py tests/test_ref_panel_builder.py && git commit -m "feat: wire duplicate-position detection into _build_chromosome() with sidecar provenance"
```
