"""Repeated identical SNP-restriction loads should parse and collapse the file once.

``ldsc ldscore`` loads the same packaged HM3 map twice per run (once for the
reference-panel restriction, once for the regression restriction). Reusing the
memoized result avoids reparsing and re-normalizing the genome-wide map, while
still returning equal keys. The cache must invalidate when the file content,
identifier mode, or genome build changes.
"""

import ldsc._kernel.identifiers as identifiers
from ldsc._kernel.identifiers import (
    clear_snp_restriction_keys_cache,
    read_snp_restriction_keys,
)


def _write_restriction(path, rows):
    path.write_text("CHR\tPOS\tA1\tA2\n" + "".join(f"{c}\t{p}\t{a1}\t{a2}\n" for c, p, a1, a2 in rows))


def _count_parses(monkeypatch):
    calls = {"n": 0}
    original = identifiers._read_snp_restriction_table

    def counting(*args, **kwargs):
        calls["n"] += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(identifiers, "_read_snp_restriction_table", counting)
    return calls


def test_repeated_identical_load_parses_file_once(tmp_path, monkeypatch):
    clear_snp_restriction_keys_cache()
    path = tmp_path / "restr.tsv"
    _write_restriction(path, [("1", "100", "A", "C"), ("1", "200", "G", "T")])
    calls = _count_parses(monkeypatch)

    first = read_snp_restriction_keys(path, "chr_pos_allele_aware", genome_build="hg19")
    second = read_snp_restriction_keys(path, "chr_pos_allele_aware", genome_build="hg19")

    assert calls["n"] == 1, "second identical load must hit the cache, not reparse"
    assert first.keys == second.keys
    assert first.match_kind == second.match_kind
    assert first.n_input_rows == second.n_input_rows
    assert first.n_retained_keys == second.n_retained_keys


def test_different_mode_recomputes(tmp_path, monkeypatch):
    clear_snp_restriction_keys_cache()
    path = tmp_path / "restr.tsv"
    _write_restriction(path, [("1", "100", "A", "C")])
    calls = _count_parses(monkeypatch)

    read_snp_restriction_keys(path, "chr_pos_allele_aware", genome_build="hg19")
    read_snp_restriction_keys(path, "chr_pos", genome_build="hg19")

    assert calls["n"] == 2, "a different identifier mode must not reuse the cached result"


def test_modified_file_recomputes(tmp_path, monkeypatch):
    clear_snp_restriction_keys_cache()
    path = tmp_path / "restr.tsv"
    _write_restriction(path, [("1", "100", "A", "C")])
    calls = _count_parses(monkeypatch)

    first = read_snp_restriction_keys(path, "chr_pos_allele_aware", genome_build="hg19")
    _write_restriction(path, [("1", "100", "A", "C"), ("2", "300", "A", "G")])
    second = read_snp_restriction_keys(path, "chr_pos_allele_aware", genome_build="hg19")

    assert calls["n"] == 2, "a changed file must be reparsed"
    assert second.keys != first.keys
    assert "2:300:A:G" in second.keys
