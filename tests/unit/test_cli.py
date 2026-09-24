"""CLI and plane-file I/O."""

from __future__ import annotations

from pathlib import Path

from vd3d.cli import format_result, main, resolve_planes, warn_if_not_general_position
from vd3d.events.samples import planes_four_through_123, planes_through_123
from vd3d.geometry.plane import Plane
from vd3d.io import dump_planes, format_plane, load_planes
from vd3d.sweep import vertical_decomposition_3d


def test_load_planes_with_and_without_ids(tmp_path: Path):
    path = tmp_path / "planes.txt"
    path.write_text(
        "# comment\n"
        "1 0 1 -4\n"
        "7 0 1 1 -5\n"
        "1 1 1 -6\n",
        encoding="utf-8",
    )
    planes = load_planes(path)
    assert [p.id for p in planes] == [1, 7, 2]
    assert planes[0].a == 1 and planes[0].d == -4
    assert planes[1].id == 7 and planes[1].b == 1


def test_load_planes_accepts_fractions(tmp_path: Path):
    path = tmp_path / "frac.txt"
    path.write_text("1 0 1 -3/2\n", encoding="utf-8")
    planes = load_planes(path)
    assert planes[0].d.numerator == -3
    assert planes[0].d.denominator == 2


def test_format_plane_equation():
    planes = planes_through_123()
    assert format_plane(planes[0]) == "x+z-4=0"
    assert format_plane(planes[1]) == "y+z-5=0"
    assert format_plane(planes[2]) == "x+y+z-6=0"
    assert format_plane(Plane(id=9, a=2, b=2, c=-4, d=0)) == "2x+2y-4z=0"


def test_dump_roundtrip(tmp_path: Path):
    planes = planes_through_123()
    path = tmp_path / "out.txt"
    dump_planes(planes, path)
    loaded = load_planes(path)
    assert [(p.id, p.a, p.b, p.c, p.d) for p in loaded] == [
        (p.id, p.a, p.b, p.c, p.d) for p in planes
    ]


def test_cli_from_file(tmp_path: Path, capsys):
    path = tmp_path / "three.txt"
    dump_planes(planes_through_123(), path)
    assert main([str(path)]) == 0
    out = capsys.readouterr().out
    assert "planes: 3" in out
    assert "seed:" not in out
    assert "P1: x+z-4=0" in out
    assert "cells:" in out
    assert "3D cells:" in out


def test_cli_random_planes_seed_reproducible(capsys):
    assert main(["--random-planes", "3", "--seed", "42"]) == 0
    first = capsys.readouterr().out
    assert main(["--random-planes", "3", "--seed", "42"]) == 0
    second = capsys.readouterr().out
    assert first == second
    assert "planes: 3" in first
    assert "seed: 42" in first


def test_warn_if_not_general_position_simultaneous(capsys):
    result = vertical_decomposition_3d(planes_four_through_123())
    warnings = warn_if_not_general_position(result)
    assert warnings
    err = capsys.readouterr().err
    assert "not in general position" in err
    assert "share z=" in err


def test_cli_seed4_n4_warns_on_stderr(capsys):
    assert main(["--random-planes", "4", "--seed", "4"]) == 0
    err = capsys.readouterr().err
    assert "not in general position" in err


def test_cli_write_planes(tmp_path: Path):
    out = tmp_path / "random.txt"
    assert main(["--random-planes", "2", "--seed", "1", "--write-planes", str(out)]) == 0
    planes = load_planes(out)
    assert len(planes) == 2


def test_cli_requires_input():
    assert main([]) == 1


def test_format_result_matches_sweep():
    result = vertical_decomposition_3d(planes_through_123())
    text = format_result(result)
    assert f"cells: {len(result.cells)}" in text
    assert f"events: {len(result.events)}" in text


def test_resolve_planes_from_namespace(tmp_path: Path):
    path = tmp_path / "p.txt"
    dump_planes(planes_through_123(), path)

    class Args:
        random_planes = None
        seed = 0
        planes_file = path

    planes = resolve_planes(Args())
    assert len(planes) == 3
