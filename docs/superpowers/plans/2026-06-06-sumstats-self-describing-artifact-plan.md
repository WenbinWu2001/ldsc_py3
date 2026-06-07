# Self-Describing Single-File Sumstats Artifact — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `munge-sumstats` write a self-describing single `.parquet` (identity in the Parquet footer, no `metadata.json`), make downstream loading tolerate metadata-less files, add a chr_pos genome-build guard in regression, and remove `schema_version` from artifact metadata package-wide.

**Architecture:** The sumstats identity payload moves from a root `metadata.json` sidecar into discrete `ldsc:*` Parquet footer keys, mirroring the existing reference-panel footer convention (`src/ldsc/_kernel/ref_panel.py`). `load_sumstats` reads the footer and no longer searches for a sidecar; absent metadata yields `config_snapshot=None` instead of an error. Regression adds a chr_pos-only build-consistency check using `infer_chr_pos_build`. `schema_version` is dropped from all writers and readers, leaving `artifact_type` as the sole identity guard.

**Tech Stack:** Python 3, pandas, pyarrow (Parquet footer key/value metadata), pytest.

**Spec:** [`docs/superpowers/specs/2026-06-06-sumstats-self-describing-artifact-design.md`](../specs/2026-06-06-sumstats-self-describing-artifact-design.md)

**Conventions for every task:**
- Run tests with the project env: `source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && pytest <args>`.
- Commit messages use Conventional Commits; no AI tool names in messages.
- TDD: write the failing test, see it fail, implement, see it pass, commit.

---

## File Structure

- `src/ldsc/_kernel/snp_identity.py` — remove `SCHEMA_VERSION`; drop the field from `identity_artifact_metadata`; drop the check from `validate_identity_artifact_metadata`. (Part 4)
- `src/ldsc/_kernel/ref_panel.py` — drop `ldsc:schema_version` from required keys / parsed dict / package-marker set / TSV-sidecar required set. (Part 4)
- `src/ldsc/outputs.py`, `src/ldsc/ref_panel_builder.py`, `src/ldsc/regression_runner.py`, `src/ldsc/annotation_builder.py` — drop `schema_version` from emitted metadata and the LD-score provenance check. (Part 4)
- `src/ldsc/sumstats_munger.py` — add footer write/read helpers; thread footer metadata into the parquet writer; remove `metadata.json` from `run`/`write_output`; rewrite `load_sumstats` metadata resolution. (Parts 1–2)
- `src/ldsc/regression_runner.py` — add `_regression_genome_build` guard, call it in `build_dataset` and `build_rg_dataset`. (Part 3)
- `docs/current/munge-sumstats.md`, `docs/current/artifact-metadata-field-inventory.md`, `README.md`, `design_map.md`, tutorials — documentation. (Part 5)
- Tests across `tests/test_snp_identity.py`, `tests/test_ref_panel.py`, `tests/test_output.py`, `tests/test_ldscore_workflow.py`, `tests/test_regression_workflow.py`, `tests/test_ref_panel_builder.py`, `tests/test_plink_io.py`, `tests/test_sumstats_munger.py`.

---

## Part 4 — Remove `schema_version` package-wide (do first)

Doing this first settles the metadata helpers that Parts 1–2 reuse.

### Task 1: Drop `schema_version` from the core identity helpers

**Files:**
- Modify: `src/ldsc/_kernel/snp_identity.py:35`, `:115-138`
- Test: `tests/test_snp_identity.py`

- [ ] **Step 1: Update the failing test**

In `tests/test_snp_identity.py`, find the assertions on `identity_artifact_metadata` / `validate_identity_artifact_metadata` and rewrite them to the new contract:

```python
def test_identity_artifact_metadata_has_no_schema_version():
    meta = identity_artifact_metadata(
        artifact_type="sumstats", snp_identifier="chr_pos_allele_aware", genome_build="hg19"
    )
    assert "schema_version" not in meta
    assert meta == {
        "artifact_type": "sumstats",
        "snp_identifier": "chr_pos_allele_aware",
        "genome_build": "hg19",
    }


def test_validate_identity_artifact_metadata_accepts_payload_without_schema_version():
    meta = {"artifact_type": "sumstats", "snp_identifier": "rsid_allele_aware", "genome_build": None}
    assert validate_identity_artifact_metadata(meta, expected_artifact_type="sumstats") == "rsid_allele_aware"


def test_validate_identity_artifact_metadata_rejects_wrong_artifact_type():
    meta = {"artifact_type": "ldscore", "snp_identifier": "rsid", "genome_build": None}
    with pytest.raises(LDSCInputError):
        validate_identity_artifact_metadata(meta, expected_artifact_type="sumstats")
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_snp_identity.py -k "schema_version or artifact_type" -v`
Expected: FAIL (current code still emits/checks `schema_version`).

- [ ] **Step 3: Implement**

In `src/ldsc/_kernel/snp_identity.py`:
- Delete line 35: `SCHEMA_VERSION = 1`.
- In `identity_artifact_metadata`, remove the `"schema_version": SCHEMA_VERSION,` entry from the returned dict.
- In `validate_identity_artifact_metadata`, replace the guard with an `artifact_type`-only check:

```python
def validate_identity_artifact_metadata(metadata: dict[str, object], *, expected_artifact_type: str) -> str:
    """Validate minimal identity metadata and return the normalized SNP identifier mode."""
    if metadata.get("artifact_type") != expected_artifact_type:
        actual_type = metadata.get("artifact_type")
        raise LDSCInputError(
            f"Could not read LDSC {expected_artifact_type} artifact metadata: expected "
            f"artifact_type={expected_artifact_type!r}, got artifact_type={actual_type!r}. "
            "Most likely the artifact was written by an older LDSC version or by another tool. "
            "Regenerate the artifact with the current LDSC package. "
            f"Other causes & fixes: {_ARTIFACT_METADATA_TROUBLESHOOTING}"
        )
    try:
        return normalize_snp_identifier_mode(str(metadata.get("snp_identifier")))
    except LDSCConfigError as exc:
        actual_mode = metadata.get("snp_identifier")
        raise LDSCInputError(
            f"Could not read LDSC {expected_artifact_type} artifact metadata: snp_identifier={actual_mode!r} is invalid. "
            "Most likely the artifact was hand-edited. Regenerate it with the current LDSC package."
        ) from exc
```

- [ ] **Step 4: Run to verify pass**

Run: `pytest tests/test_snp_identity.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/snp_identity.py tests/test_snp_identity.py
git commit -m "refactor(metadata): drop schema_version from identity helpers"
```

### Task 2: Drop `ldsc:schema_version` from the ref-panel reader

**Files:**
- Modify: `src/ldsc/_kernel/ref_panel.py:120-139`, `:151-157`, `:915-930`
- Test: `tests/test_ref_panel.py`

- [ ] **Step 1: Update the failing test**

In `tests/test_ref_panel.py`, adjust any assertion expecting a `schema_version` key in parsed footer/sidecar metadata to expect only `artifact_type`, `snp_identifier`, `genome_build`. Add:

```python
def test_ref_panel_footer_reads_without_schema_version(tmp_path):
    # build a minimal R2 parquet via the package writer used elsewhere in this module,
    # then assert the parsed identity dict has no schema_version and still validates.
    metadata = _read_r2_identity_schema_metadata(path, expected_artifact_type="ref_panel_r2")
    assert "schema_version" not in metadata
    assert metadata["artifact_type"] == "ref_panel_r2"
```

(Reuse the existing fixture/helper that the surrounding tests use to produce a package R2 parquet; do not invent a new builder.)

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_ref_panel.py -k "schema_version or identity" -v`
Expected: FAIL.

- [ ] **Step 3: Implement**

In `src/ldsc/_kernel/ref_panel.py`:
- In the footer reader (lines ~120-139): remove `b"ldsc:schema_version"` from `required_keys`, and remove the `"schema_version": int(raw[b"ldsc:schema_version"]...)` entry from the parsed `metadata` dict.
- In `_r2_path_has_ldsc_package_schema` (lines ~151-157): remove `b"ldsc:schema_version"` from `package_keys` (keep `b"ldsc:artifact_type"`, `b"ldsc:sorted_by_build"`, `b"ldsc:row_group_size"`).
- In the TSV-sidecar metadata reader (lines ~915-930): remove `"schema_version"` from the `required` set and the `"schema_version": int(raw["schema_version"])` entry.

- [ ] **Step 4: Run to verify pass**

Run: `pytest tests/test_ref_panel.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/ref_panel.py tests/test_ref_panel.py
git commit -m "refactor(metadata): drop schema_version from ref-panel readers"
```

### Task 3: Drop `schema_version` from result/annotation writers and the LD-score provenance check

**Files:**
- Modify: `src/ldsc/outputs.py:671`, `:872`, `:913`; `src/ldsc/ref_panel_builder.py:330`; `src/ldsc/regression_runner.py:1342`, `:2194-2196`; `src/ldsc/annotation_builder.py:1082`
- Test: `tests/test_output.py`, `tests/test_ldscore_workflow.py`, `tests/test_regression_workflow.py`, `tests/test_ref_panel_builder.py`, `tests/test_plink_io.py`

- [ ] **Step 1: Update the failing tests**

In each listed test file, remove assertions that a written `metadata.json` (or footer/sidecar) contains `schema_version`, and where a test asserts the full metadata dict, drop the `schema_version` key from the expected value. Add one representative assertion in `tests/test_output.py`:

```python
def test_h2_result_metadata_has_no_schema_version(tmp_path):
    # write an h2 result dir using the existing H2DirectoryWriter fixture in this file
    payload = json.loads((result_dir / "diagnostics" / "metadata.json").read_text())
    assert "schema_version" not in payload
    assert payload["artifact_type"] == "h2_result"
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_output.py tests/test_ldscore_workflow.py tests/test_regression_workflow.py tests/test_ref_panel_builder.py tests/test_plink_io.py -k schema_version -v`
Expected: FAIL.

- [ ] **Step 3: Implement**

- `src/ldsc/outputs.py`: delete the `"schema_version": 1,` line from the partitioned-h2 payload (line ~671), the rg pair payload (line ~872), and `payload["schema_version"] = 1` in `_result_metadata` (line ~913).
- `src/ldsc/ref_panel_builder.py`: delete `"schema_version": 1,` (line ~330).
- `src/ldsc/annotation_builder.py`: delete `"schema_version": 1,` (line ~1082).
- `src/ldsc/regression_runner.py`: delete `"schema_version": 1,` (line ~1342). Update the LD-score provenance check (lines ~2194-2196) to require `artifact_type` only:

```python
    if "artifact_type" not in metadata:
        raise LDSCInputError(
            "Regression could not read LD-score artifact provenance: metadata lacks `artifact_type`. "
            "Most likely the LD-score directory was written by an older LDSC version. "
            f"Regenerate it with the current `ldsc ldscore`. Other causes & fixes: {_REGRESSION_SCHEMA_DOC}"
        )
```

- [ ] **Step 4: Run to verify pass**

Run: `pytest tests/test_output.py tests/test_ldscore_workflow.py tests/test_regression_workflow.py tests/test_ref_panel_builder.py tests/test_plink_io.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/outputs.py src/ldsc/ref_panel_builder.py src/ldsc/annotation_builder.py src/ldsc/regression_runner.py tests/
git commit -m "refactor(metadata): drop schema_version from result and annotation metadata"
```

---

## Part 1 — Munge writes a self-describing single `.parquet`

### Task 4: Add the sumstats footer builder + thread it into the parquet writer

**Files:**
- Modify: `src/ldsc/sumstats_munger.py` (`_write_sumstats_parquet` ~1638, `_write_sumstats_outputs` ~1690)
- Test: `tests/test_sumstats_munger.py`

- [ ] **Step 1: Write the failing test**

```python
def test_sumstats_parquet_embeds_identity_footer(tmp_path):
    import pyarrow.parquet as pq
    from ldsc.config import GlobalConfig
    from ldsc.sumstats_munger import _sumstats_footer_metadata, _write_sumstats_outputs

    data = pd.DataFrame({"SNP": ["rs1", "rs2"], "CHR": ["1", "1"], "POS": [10, 20],
                         "A1": ["A", "C"], "A2": ["G", "T"], "Z": [0.1, 0.2], "N": [100, 100]})
    snapshot = GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg19")
    footer = _sumstats_footer_metadata(snapshot, trait_name="height")
    out = str(tmp_path / "sumstats.parquet")
    _write_sumstats_outputs(data, output_files={"parquet": out}, output_format="parquet",
                            footer_metadata=footer)

    raw = pq.read_schema(out).metadata
    assert raw[b"ldsc:artifact_type"] == b"sumstats"
    assert raw[b"ldsc:snp_identifier"] == b"chr_pos_allele_aware"
    assert raw[b"ldsc:genome_build"] == b"hg19"
    assert raw[b"ldsc:trait_name"] == b"height"
    assert b"ldsc:schema_version" not in raw


def test_sumstats_footer_encodes_none_genome_build_as_empty():
    from ldsc.config import GlobalConfig
    from ldsc.sumstats_munger import _sumstats_footer_metadata
    footer = _sumstats_footer_metadata(GlobalConfig(snp_identifier="rsid_allele_aware"), trait_name=None)
    assert footer[b"ldsc:genome_build"] == b""
    assert footer[b"ldsc:trait_name"] == b""
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_sumstats_munger.py -k footer -v`
Expected: FAIL (`_sumstats_footer_metadata` not defined; `_write_sumstats_outputs` has no `footer_metadata` arg).

- [ ] **Step 3: Implement**

In `src/ldsc/sumstats_munger.py`, add the builder (place it near `_write_sumstats_parquet`):

```python
def _sumstats_footer_metadata(config_snapshot: GlobalConfig, trait_name: str | None) -> dict[bytes, bytes]:
    """Return discrete ``ldsc:*`` Parquet footer keys for a curated sumstats artifact.

    Mirrors the reference-panel footer convention in ``_kernel/ref_panel.py``. A
    ``None`` ``genome_build`` (rsID-family) or ``trait_name`` is encoded as an
    empty string and decoded back to ``None`` on read.
    """
    identity = identity_artifact_metadata(
        artifact_type="sumstats",
        snp_identifier=config_snapshot.snp_identifier,
        genome_build=config_snapshot.genome_build,
    )
    footer = {
        f"ldsc:{key}".encode("utf-8"): ("" if value is None else str(value)).encode("utf-8")
        for key, value in identity.items()
    }
    footer[b"ldsc:trait_name"] = ("" if trait_name is None else str(trait_name)).encode("utf-8")
    return footer
```

Thread it through the writer. Change `_write_sumstats_parquet` to accept the footer and merge it into the schema metadata:

```python
def _write_sumstats_parquet(
    data: pd.DataFrame, path: str, footer_metadata: dict[bytes, bytes] | None = None
) -> list[dict[str, Any]]:
    ...
    frame = _prepare_sumstats_parquet_frame(data)
    schema = pa.Schema.from_pandas(frame, preserve_index=False)
    if footer_metadata:
        merged = {**(schema.metadata or {}), **footer_metadata}
        schema = schema.with_metadata(merged)
    ...
```

Change `_write_sumstats_outputs` to accept and forward it:

```python
def _write_sumstats_outputs(
    data: pd.DataFrame,
    *,
    output_files: dict[str, str],
    output_format: str,
    footer_metadata: dict[bytes, bytes] | None = None,
) -> tuple[str, list[dict[str, Any]]]:
    output_format = _normalize_output_format(output_format)
    parquet_row_groups: list[dict[str, Any]] = []
    if "tsv.gz" in output_files:
        _write_sumstats_tsv_gz(data, output_files["tsv.gz"])
    if "parquet" in output_files:
        parquet_row_groups = _write_sumstats_parquet(data, output_files["parquet"], footer_metadata)
    primary = output_files["parquet"] if output_format in {"parquet", "both"} else output_files["tsv.gz"]
    return primary, parquet_row_groups
```

- [ ] **Step 4: Run to verify pass**

Run: `pytest tests/test_sumstats_munger.py -k footer -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/sumstats_munger.py tests/test_sumstats_munger.py
git commit -m "feat(munge): embed identity metadata in sumstats parquet footer"
```

### Task 5: Remove `metadata.json` from `run` and `write_output`

**Files:**
- Modify: `src/ldsc/sumstats_munger.py` (`run` ~507-616, `write_output` ~685-705)
- Test: `tests/test_sumstats_munger.py`

- [ ] **Step 1: Write the failing test**

```python
def test_munge_run_writes_no_metadata_json(tmp_path, tiny_raw_sumstats_path):
    # Use the existing end-to-end munge fixture/helper in this test module.
    out = tmp_path / "munged"
    run_munge(raw=tiny_raw_sumstats_path, output_dir=str(out))  # existing helper
    assert (out / "sumstats.parquet").exists()
    assert not (out / "metadata.json").exists()
    assert (out / "diagnostics" / "sumstats.log").exists()
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_sumstats_munger.py -k no_metadata_json -v`
Expected: FAIL (`metadata.json` still written).

- [ ] **Step 3: Implement**

In `run`:
- Delete `metadata_path = output_dir / "metadata.json"` (line ~509).
- Remove `metadata_path` from `produced_paths` (line ~515) and `owned_paths` (line ~518).
- Build the footer and pass it to `_write_sumstats_outputs`:

```python
            footer_metadata = _sumstats_footer_metadata(table_config_snapshot, raw_sumstats_config.trait_name)
            primary_sumstats_file, parquet_row_groups = _write_sumstats_outputs(
                data,
                output_files=output_files,
                output_format=munge_config.output_format,
                footer_metadata=footer_metadata,
            )
```

  Note: `table_config_snapshot` is computed at line ~571 (`_effective_sumstats_config`). Move the `_write_sumstats_outputs` call to *after* that line so the footer carries the effective (post-liftover) build, or compute `table_config_snapshot` before writing. Keep the existing ordering of `coordinate_metadata`/`drop_frame` derivation; only the parquet write must follow `table_config_snapshot`.
- Delete the `_write_sumstats_metadata(metadata_path, ...)` call (lines ~587-592).
- Remove `"metadata_json": str(metadata_path),` from `run_output_paths` (line ~614) and from the `provenance` dict's `"metadata_path"` entry (line ~604).

In `write_output`:
- Delete `metadata_path = output_root / "metadata.json"` (line ~685).
- Remove `metadata_path` from `produced_paths` and the owned list (lines ~688, ~691).
- Replace the `_write_sumstats_metadata(...)` call (lines ~700-705) with footer threading:

```python
        primary_sumstats_file, _parquet_row_groups = _write_sumstats_outputs(
            sumstats.data,
            output_files=output_files,
            output_format=output_format,
            footer_metadata=_sumstats_footer_metadata(sumstats.config_snapshot, sumstats.trait_name),
        )
```

- Delete the now-unused `_write_sumstats_metadata` function and `_sumstats_metadata_path` if no remaining callers (verify with grep before deleting).

- [ ] **Step 4: Run to verify pass**

Run: `pytest tests/test_sumstats_munger.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/sumstats_munger.py tests/test_sumstats_munger.py
git commit -m "feat(munge): stop writing metadata.json sidecar"
```

---

## Part 2 — Loader: footer-or-nothing, metadata optional

### Task 6: Add the footer reader

**Files:**
- Modify: `src/ldsc/sumstats_munger.py` (near `_read_curated_sumstats_artifact`)
- Test: `tests/test_sumstats_munger.py`

- [ ] **Step 1: Write the failing test**

```python
def test_read_sumstats_parquet_footer_roundtrip(tmp_path):
    from ldsc.config import GlobalConfig
    from ldsc.sumstats_munger import _sumstats_footer_metadata, _write_sumstats_outputs, _read_sumstats_parquet_footer

    data = pd.DataFrame({"SNP": ["rs1"], "CHR": ["1"], "POS": [10], "A1": ["A"], "A2": ["G"], "Z": [0.1], "N": [100]})
    out = str(tmp_path / "sumstats.parquet")
    footer = _sumstats_footer_metadata(GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg19"), "height")
    _write_sumstats_outputs(data, output_files={"parquet": out}, output_format="parquet", footer_metadata=footer)

    meta = _read_sumstats_parquet_footer(out)
    assert meta == {"artifact_type": "sumstats", "snp_identifier": "chr_pos_allele_aware",
                    "genome_build": "hg19", "trait_name": "height"}


def test_read_sumstats_parquet_footer_absent_returns_none(tmp_path):
    from ldsc.sumstats_munger import _read_sumstats_parquet_footer
    out = str(tmp_path / "plain.parquet")
    pd.DataFrame({"SNP": ["rs1"], "Z": [0.1], "N": [100]}).to_parquet(out)
    assert _read_sumstats_parquet_footer(out) is None


def test_read_sumstats_parquet_footer_non_parquet_returns_none(tmp_path):
    from ldsc.sumstats_munger import _read_sumstats_parquet_footer
    out = str(tmp_path / "x.sumstats.gz")
    pd.DataFrame({"SNP": ["rs1"], "Z": [0.1], "N": [100]}).to_csv(out, sep="\t", index=False, compression="gzip")
    assert _read_sumstats_parquet_footer(out) is None
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_sumstats_munger.py -k footer -v`
Expected: FAIL (`_read_sumstats_parquet_footer` not defined).

- [ ] **Step 3: Implement**

```python
def _read_sumstats_parquet_footer(path: str) -> dict[str, Any] | None:
    """Return the embedded sumstats identity payload, or ``None`` when absent.

    Only Parquet artifacts carry a footer. ``.sumstats(.gz)`` and footer-less
    Parquet files return ``None``. Empty ``genome_build``/``trait_name`` footer
    values decode back to ``None``.
    """
    if not str(path).endswith(".parquet"):
        return None
    import pyarrow.parquet as pq

    raw = pq.read_schema(path).metadata or {}
    if b"ldsc:artifact_type" not in raw:
        return None
    genome_build = raw.get(b"ldsc:genome_build", b"").decode("utf-8") or None
    trait_name = raw.get(b"ldsc:trait_name", b"").decode("utf-8") or None
    return {
        "artifact_type": raw[b"ldsc:artifact_type"].decode("utf-8"),
        "snp_identifier": raw[b"ldsc:snp_identifier"].decode("utf-8"),
        "genome_build": genome_build,
        "trait_name": trait_name,
    }
```

- [ ] **Step 4: Run to verify pass**

Run: `pytest tests/test_sumstats_munger.py -k footer -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/sumstats_munger.py tests/test_sumstats_munger.py
git commit -m "feat(loader): read identity metadata from sumstats parquet footer"
```

### Task 7: Rewrite `load_sumstats` metadata resolution (footer-or-none)

**Files:**
- Modify: `src/ldsc/sumstats_munger.py` (`load_sumstats` ~338-415); remove sidecar helpers `_sumstats_metadata_path`/`_read_sumstats_metadata` if unused.
- Test: `tests/test_sumstats_munger.py`

- [ ] **Step 1: Write the failing tests**

```python
def test_load_sumstats_footer_only_parquet(tmp_path, munged_chr_pos_parquet_path):
    # munged_chr_pos_parquet_path: a parquet written by the new munge (footer, no metadata.json)
    table = load_sumstats(munged_chr_pos_parquet_path)
    assert table.config_snapshot is not None
    assert table.config_snapshot.snp_identifier == "chr_pos_allele_aware"


def test_load_sumstats_legacy_sumstats_gz_has_no_config(tmp_path):
    out = tmp_path / "legacy.sumstats.gz"
    pd.DataFrame({"SNP": ["rs1", "rs2"], "A1": ["A", "C"], "A2": ["G", "T"],
                  "Z": [0.1, 0.2], "N": [100, 100]}).to_csv(out, sep="\t", index=False, compression="gzip")
    table = load_sumstats(str(out))
    assert table.config_snapshot is None
    assert set(["SNP", "Z", "N"]).issubset(table.data.columns)


def test_load_sumstats_footerless_parquet_has_no_config(tmp_path):
    out = tmp_path / "plain.parquet"
    pd.DataFrame({"SNP": ["rs1"], "A1": ["A"], "A2": ["G"], "Z": [0.1], "N": [100]}).to_parquet(out)
    table = load_sumstats(str(out))
    assert table.config_snapshot is None
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_sumstats_munger.py -k load_sumstats -v`
Expected: FAIL (current loader raises on missing `metadata.json`).

- [ ] **Step 3: Implement**

Replace the metadata block in `load_sumstats` (the section from `metadata_path = _sumstats_metadata_path(resolved)` through the end of the identity-cleanup validation, lines ~352-403) with:

```python
    metadata = _read_sumstats_parquet_footer(resolved)
    if metadata is None:
        LOGGER.info(
            f"No embedded identity metadata in '{resolved}'. The SNP identifier mode will be "
            "inferred from the LD-score panel during regression."
        )
        table = SumstatsTable(
            data=df.reset_index(drop=True),
            has_alleles={"A1", "A2"}.issubset(df.columns),
            source_path=resolved,
            trait_name=_resolve_sumstats_trait_name(trait_name, None, resolved),
            provenance={},
            config_snapshot=None,
        )
        table.validate()
        return table
    config_snapshot = _global_config_from_sumstats_metadata(metadata, artifact_path=resolved)
    mode = normalize_snp_identifier_mode(config_snapshot.snp_identifier)
    # ... existing allele-aware check, clean_identity_artifact_table, effective_merge_key_series
    #     validation block stays unchanged from here ...
```

Keep the existing construction of the metadata-present `SumstatsTable` (lines ~404-415) but set `provenance={"metadata": metadata}` (drop the `metadata_path` key, since there is no sidecar path). After editing, grep for `_sumstats_metadata_path` and `_read_sumstats_metadata`; delete them if they now have zero callers.

- [ ] **Step 4: Run to verify pass**

Run: `pytest tests/test_sumstats_munger.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/sumstats_munger.py tests/test_sumstats_munger.py
git commit -m "feat(loader): make sumstats metadata optional, drop sidecar lookup"
```

---

## Part 3 — Regression genome-build guard (chr_pos only)

### Task 8: Add `_regression_genome_build` and wire it into h2 + rg

**Files:**
- Modify: `src/ldsc/regression_runner.py` (new helper near `_resolve_h2_identity` ~949; calls in `build_dataset` ~292 and `build_rg_dataset` ~435)
- Test: `tests/test_regression_workflow.py`

- [ ] **Step 1: Write the failing tests**

```python
def test_regression_build_guard_errors_on_mismatch():
    from ldsc.regression_runner import _regression_genome_build
    from ldsc.config import GlobalConfig
    sumstats = _make_sumstats_table(snp_identifier="chr_pos_allele_aware", genome_build="hg19")  # helper in this file
    ldscore = _make_ldscore_result(snp_identifier="chr_pos_allele_aware", genome_build="hg38")
    with pytest.raises(LDSCInputError, match="liftover"):
        _regression_genome_build(sumstats, ldscore, "chr_pos_allele_aware", context="sumstats")


def test_regression_build_guard_skips_rsid():
    from ldsc.regression_runner import _regression_genome_build
    sumstats = _make_sumstats_table(snp_identifier=None, genome_build=None, config_snapshot=None)
    ldscore = _make_ldscore_result(snp_identifier="rsid_allele_aware", genome_build=None)
    # rsID family: returns without inference or error even though sumstats has no build
    _regression_genome_build(sumstats, ldscore, "rsid_allele_aware", context="sumstats")


def test_regression_build_guard_infers_when_metadata_absent(chr_pos_sumstats_table_no_metadata, hg19_ldscore_result):
    from ldsc.regression_runner import _regression_genome_build
    # chr_pos sumstats whose coordinates are hg19; panel is hg19 -> passes via inference
    _regression_genome_build(chr_pos_sumstats_table_no_metadata, hg19_ldscore_result, "chr_pos", context="sumstats")
```

(Reuse or add small local helpers `_make_sumstats_table` / `_make_ldscore_result` consistent with existing fixtures in `tests/test_regression_workflow.py`.)

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_regression_workflow.py -k build_guard -v`
Expected: FAIL (`_regression_genome_build` not defined).

- [ ] **Step 3: Implement**

Add near the identity resolvers in `src/ldsc/regression_runner.py` (import `infer_chr_pos_build` and `identity_mode_family` at top if not already):

```python
def _regression_genome_build(
    sumstats_table: SumstatsTable,
    ldscore_result: LDScoreResult,
    identifier_mode: str,
    *,
    context: str,
) -> None:
    """Verify a chr_pos sumstats shares the LD-score panel's genome build.

    rsID-family runs skip this entirely. When the sumstats carries no build
    metadata, the build is inferred from its coordinates; an inconclusive
    inference or a build mismatch is a hard error directing the user to liftover.
    """
    if identity_mode_family(identifier_mode) != "chr_pos":
        return
    panel_build = getattr(ldscore_result.config_snapshot, "genome_build", None)
    snapshot = sumstats_table.config_snapshot
    sumstats_build = getattr(snapshot, "genome_build", None) if snapshot is not None else None
    if sumstats_build is None:
        try:
            inference = infer_chr_pos_build(sumstats_table.data, context=context)
        except ValueError as exc:
            raise LDSCInputError(
                f"Regression could not verify the genome build of sumstats '{context}': its coordinates "
                "did not match the reference closely enough to infer a build, and it carries no build metadata. "
                "Most likely it is a legacy or third-party file. Re-munge it with the current `ldsc munge-sumstats` "
                "so it records its build, or liftover it to match the LD-score panel, then re-run."
            ) from exc
        sumstats_build = inference.genome_build
        LOGGER.info(f"Inferred genome build '{sumstats_build}' for sumstats '{context}'.")
    if panel_build is not None and sumstats_build != panel_build:
        raise LDSCInputError(
            f"Regression genome-build mismatch: sumstats '{context}' is {sumstats_build!r} but the LD-score "
            f"panel is {panel_build!r}. Coordinate merges require the same build. Liftover the sumstats to "
            f"{panel_build!r} and re-run."
        )
    LOGGER.info(f"Regression genome build for '{context}': {sumstats_build!r} (panel {panel_build!r}).")
```

In `build_dataset`, right after `identifier_mode` is resolved (after line ~292, before the merge branches at ~328):

```python
        _regression_genome_build(
            sumstats_table, ldscore_result, identifier_mode,
            context=sumstats_table.source_path or sumstats_table.trait_name or "sumstats",
        )
```

In `build_rg_dataset`, after its `identifier_mode` is resolved (~435), call it once per sumstats table:

```python
        for _table in (sumstats_table_1, sumstats_table_2):
            _regression_genome_build(
                _table, ldscore_result, identifier_mode,
                context=_table.source_path or _table.trait_name or "sumstats",
            )
```

- [ ] **Step 4: Run to verify pass**

Run: `pytest tests/test_regression_workflow.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/regression_runner.py tests/test_regression_workflow.py
git commit -m "feat(regression): add chr_pos genome-build guard with liftover guidance"
```

---

## Part 5 — Documentation

### Task 9: Rewrite the munge-sumstats module doc

**Files:**
- Modify: `docs/current/munge-sumstats.md`

- [ ] **Step 1: Update the Output Artifacts tree**

Replace the tree at lines ~75-77 (remove `metadata.json`):

```
<output_dir>/
  sumstats.parquet
  sumstats.sumstats.gz        # only with --output-format tsv.gz or both
  diagnostics/
    sumstats.log
    dropped_snps/dropped.tsv.gz
```

- [ ] **Step 2: Rewrite the "Root Metadata Schema" section (lines ~93-112) as embedded footer metadata**

Replace the entire section with a clearly-titled description of the included metadata:

```markdown
## Embedded Identity Metadata

The curated `sumstats.parquet` is self-describing: its downstream identity payload
is written into the Parquet **footer** as discrete `ldsc:*` key/value entries, so
the data file alone is sufficient for regression. No `metadata.json` sidecar is
written.

| Footer key | Meaning |
|---|---|
| `ldsc:artifact_type` | Always `sumstats`. The sole identity guard. |
| `ldsc:snp_identifier` | Identity mode used when writing the artifact. |
| `ldsc:genome_build` | Final output build in coordinate-family modes; empty for rsID-family modes (decoded as `null`). |
| `ldsc:trait_name` | Trait label from the raw sumstats config when supplied; empty otherwise. |

There is no `schema_version`. `.sumstats.gz` is a plain TSV and carries no embedded
metadata; in `both` mode only the `.parquet` is self-describing. When a loaded file
has no footer metadata (legacy `.sumstats.gz` or a footer-less parquet), regression
infers the identifier mode from the LD-score panel and logs it.
```

- [ ] **Step 3: Fix scattered references**

- Line ~31 (workflow step 9): change "Write fixed outputs, root metadata, and diagnostic audit files." to "Write the self-describing `sumstats.parquet` (identity in its footer) and diagnostic audit files."
- Lines ~67-68: change "`metadata.json` and in-memory `SumstatsTable.config_snapshot` store `genome_build=null`..." to "The Parquet footer `ldsc:genome_build` key and in-memory `SumstatsTable.config_snapshot` store an empty/`null` build for rsID-family artifacts."
- Lines ~108-109: update the downstream-compatibility note to describe the new chr_pos genome-build guard (build from footer or inferred from coordinates; mismatch → liftover error).

- [ ] **Step 4: Commit**

```bash
git add docs/current/munge-sumstats.md
git commit -m "docs: document self-describing sumstats parquet footer metadata"
```

### Task 10: Update the artifact metadata field inventory

**Files:**
- Modify: `docs/current/artifact-metadata-field-inventory.md`

- [ ] **Step 1: Update the sumstats section (lines ~48-72)**

- Remove `metadata.json` from the `munge-sumstats` output tree (line ~52).
- Replace the field table: drop the `schema_version` and `files` rows; reframe the remaining rows as footer keys (`ldsc:artifact_type`, `ldsc:snp_identifier`, `ldsc:genome_build`, `ldsc:trait_name`); state that the parquet footer is the downstream-required identity source and no `metadata.json` is written.

- [ ] **Step 2: Update the boundary rule and per-artifact tables**

- Lines ~3-24: update wording that sumstats identity lives in `sumstats/metadata.json` to "the `sumstats.parquet` footer"; state `artifact_type` (not `schema_version` + `artifact_type`) is the identity marker.
- Remove every `schema_version` row across the per-artifact tables: lines ~67, ~121, ~136 (`ldsc:schema_version`), ~159, ~186, ~217, and any others (grep `schema_version` in this file and delete each row).

- [ ] **Step 3: Verify no stray `schema_version` remains**

Run: `grep -n "schema_version" docs/current/artifact-metadata-field-inventory.md`
Expected: no output.

- [ ] **Step 4: Commit**

```bash
git add docs/current/artifact-metadata-field-inventory.md
git commit -m "docs: drop schema_version and metadata.json from sumstats metadata inventory"
```

### Task 11: Update README, design_map, tutorials, and sweep stray references

**Files:**
- Modify: `README.md` (lines ~101-109), `design_map.md`, `tutorials/munge-sumstats.ipynb`, `tutorials/heritability-estimates.md`, `tutorials/cross-trait-genetic-correlation.md`, plus sumstats references in `docs/current/{data-flow,architecture,class-and-features,code-structure,path-specification,snp-identifier-genome-build-defaults}.md`.

- [ ] **Step 1: README**

Rewrite lines ~101-109 to state munge writes a self-describing `sumstats.parquet` (footer metadata, no `sumstats.metadata.json` sidecar). Remove the "writes a thin `sumstats.metadata.json` beside the selected output" sentence.

- [ ] **Step 2: design_map.md**

Add a row mapping the new design doc to its implementation:

```markdown
| `docs/superpowers/specs/2026-06-06-sumstats-self-describing-artifact-design.md` | sumstats parquet footer metadata + `load_sumstats` footer reader (`sumstats_munger.py`), chr_pos genome-build guard (`regression_runner._regression_genome_build`), package-wide `schema_version` removal |
```

- [ ] **Step 3: Sweep sumstats-specific `metadata.json` references**

For each listed `docs/current/*.md` and tutorial, grep for `metadata.json` / `sumstats.metadata.json` / `schema_version` and update only the references that describe the **sumstats** artifact (leave ldscore/h2/rg `metadata.json` references intact — those artifacts keep their sidecars). For sumstats, point to the parquet footer.

Run per file, e.g.:
`grep -nE "metadata\.json|schema_version" docs/current/data-flow.md`
and edit sumstats-related hits.

- [ ] **Step 4: Verify**

Run: `grep -rn "sumstats.metadata.json\|sumstats/metadata.json" docs/ README.md tutorials/`
Expected: no output (all sumstats sidecar references removed).

- [ ] **Step 5: Commit**

```bash
git add README.md design_map.md docs/current/ tutorials/
git commit -m "docs: update README, design map, and tutorials for self-describing sumstats"
```

---

## Final verification

### Task 12: Full suite + grep gates

- [ ] **Step 1: Run the whole test suite**

Run: `source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && pytest`
Expected: all green.

- [ ] **Step 2: Confirm `schema_version` is gone from source**

Run: `grep -rn "schema_version\|SCHEMA_VERSION" src/ldsc/`
Expected: no output.

- [ ] **Step 3: Confirm no sumstats `metadata.json` writer remains**

Run: `grep -rn "_write_sumstats_metadata\|metadata.json" src/ldsc/sumstats_munger.py`
Expected: no output.

- [ ] **Step 4: Update lessons.md only if a recurring mistake surfaced during execution** (per CLAUDE.md). Otherwise skip.
```
