"""Tests for `rustscenic.pipeline._load_tfs` species-shortcut handling.

Regression test for the v0.4.0 bug: passing ``tfs="hs"`` (or "mm", or any
of the aliases ``data.tfs()`` accepts) was treated as a filesystem path,
producing ``FileNotFoundError: [Errno 2] No such file or directory: 'hs'``.
The README documents the species shortcut as the default zero-config path,
so a user following the docs hit a hard crash before any compute ran.

The fix is in ``rustscenic.pipeline._load_tfs``: detect species aliases
before falling through to ``Path.read_text``.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from rustscenic.pipeline import _load_tfs


def test_load_tfs_none_defaults_to_hs():
    """``tfs=None`` returns the bundled human list."""
    out = _load_tfs(None)
    assert isinstance(out, list)
    assert len(out) > 1000  # 1,839 in the v0.4 bundle
    assert all(isinstance(t, str) for t in out)


@pytest.mark.parametrize(
    "alias,expected_min",
    [
        ("hs", 1000),
        ("mm", 1000),
        ("human", 1000),
        ("mouse", 1000),
        ("hg38", 1000),
        ("mm10", 1000),
        ("HS", 1000),  # case-insensitivity
        ("Mouse", 1000),
        ("homo_sapiens", 1000),
        ("mus_musculus", 1000),
    ],
)
def test_load_tfs_species_shortcut(alias, expected_min):
    """Every species alias ``data.tfs`` accepts must also work via pipeline.run.

    Regression for the FileNotFoundError observed on v0.4.0 PyPI install.
    """
    out = _load_tfs(alias)
    assert isinstance(out, list)
    assert len(out) >= expected_min, f"alias {alias!r} returned only {len(out)} TFs"
    assert all(isinstance(t, str) for t in out)


def test_load_tfs_hs_and_mm_differ():
    """``hs`` and ``mm`` must not collapse to the same list, and contain
    known canonical TFs (guards against an accidental file swap)."""
    hs = _load_tfs("hs")
    mm = _load_tfs("mm")
    assert hs != mm
    # Known canonical symbols: SPI1 / PAX5 are core PBMC lineage TFs in HGNC,
    # Pax6 / Sox2 are core neural TFs in MGI. Spot-check that the file shipped
    # with the wheel actually contains them (catches a truncated or swapped
    # bundle).
    for sym in ("SPI1", "PAX5", "GATA3"):
        assert sym in hs, f"expected {sym!r} in hs TF list, got {len(hs)} TFs"
    for sym in ("Pax6", "Sox2", "Neurod2"):
        assert sym in mm, f"expected {sym!r} in mm TF list, got {len(mm)} TFs"


def test_load_tfs_path_with_alias_name_routes_to_shortcut(tmp_path):
    """A ``Path`` whose string is a species alias hits the shortcut, not the
    file reader. Documents that the alias check applies to both ``str`` and
    ``Path`` inputs (audit follow-up after the v0.4.0 -> v0.4.1 fix)."""
    out = _load_tfs(Path("hs"))
    assert isinstance(out, list)
    assert len(out) > 1000
    assert "SPI1" in out


def test_load_tfs_path(tmp_path):
    """A real file path still works (the original code path)."""
    p = tmp_path / "tfs.txt"
    p.write_text("FOO\nBAR\n  BAZ  \n\n")
    out = _load_tfs(str(p))
    assert out == ["FOO", "BAR", "BAZ"]


def test_load_tfs_path_object(tmp_path):
    """``Path`` instances are read as files (not species shortcuts)."""
    p = tmp_path / "tfs.txt"
    p.write_text("FOO\nBAR\n")
    out = _load_tfs(p)
    assert out == ["FOO", "BAR"]


def test_load_tfs_list():
    """Passing a list returns a copy of the list."""
    inp = ["TF1", "TF2", "TF3"]
    out = _load_tfs(inp)
    assert out == inp


def test_load_tfs_unknown_string_still_treated_as_path(tmp_path):
    """Unknown strings keep the original path-reading semantics."""
    p = tmp_path / "missing.txt"  # does not exist
    with pytest.raises(FileNotFoundError):
        _load_tfs(str(p))
