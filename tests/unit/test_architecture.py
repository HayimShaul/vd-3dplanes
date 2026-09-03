import ast
from pathlib import Path

import vd3d

_REPO = Path(vd3d.__file__).resolve().parent
_KERNEL_PACKAGES = (
    "geometry",
    "arrangement2d",
    "vertical_decomposition",
    "zone",
    "events",
    "sweep",
    "cells3d",
)
_TWOD_PACKAGES = (
    "arrangement2d",
    "vertical_decomposition",
    "zone",
)
_FORBIDDEN_IN_KERNEL = ("matplotlib", "vd3d.viz", "numpy")
_FORBIDDEN_IN_2D = ("vd3d.sweep", "vd3d.cells3d")


def _iter_python_files(package: str):
    root = _REPO / package
    yield from root.rglob("*.py")


def _imported_modules(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module)
    return names


def _matches_forbidden(imported: str, forbidden: str) -> bool:
    return imported == forbidden or imported.startswith(forbidden + ".")


def test_kernel_packages_do_not_import_viz_or_matplotlib():
    offenders = []
    for package in _KERNEL_PACKAGES:
        for path in _iter_python_files(package):
            for imported in _imported_modules(path):
                for needle in _FORBIDDEN_IN_KERNEL:
                    if _matches_forbidden(imported, needle):
                        rel = path.relative_to(_REPO)
                        offenders.append(f"{rel} imports {imported}")
    assert offenders == []


def test_2d_packages_do_not_import_sweep_or_cells3d():
    offenders = []
    for package in _TWOD_PACKAGES:
        for path in _iter_python_files(package):
            for imported in _imported_modules(path):
                for needle in _FORBIDDEN_IN_2D:
                    if _matches_forbidden(imported, needle):
                        rel = path.relative_to(_REPO)
                        offenders.append(f"{rel} imports {imported}")
    assert offenders == []


def test_viz_to_float_is_a_python_float():
    from fractions import Fraction

    from vd3d.viz.convert import to_float

    value = to_float(Fraction(1, 2))
    assert type(value) is float
    assert value == 0.5
