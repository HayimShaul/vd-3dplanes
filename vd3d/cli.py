"""Compute the 3D vertical decomposition of a list of planes.

Examples::

    python -m vd3d examples/three_planes.txt
    python -m vd3d --random-planes 5 --seed 42
"""

from __future__ import annotations

import argparse
import random
import sys
from collections.abc import Sequence
from pathlib import Path

from vd3d.cells3d.types import Cell3D
from vd3d.geometry.plane import Plane
from vd3d.geometry.sampling import random_general_position_planes
from vd3d.io import dump_planes, format_plane, load_planes
from vd3d.sweep import SweepResult, group_events_by_z, vertical_decomposition_3d


def _plane_ref(plane: Plane | None) -> str:
    if plane is None:
        return "unbounded"
    return f"P{plane.id}:{format_plane(plane)}"


def _z_lower(value) -> str:
    return "-∞" if value is None else str(value)


def _z_upper(value) -> str:
    return "+∞" if value is None else str(value)


def format_cell(cell: Cell3D) -> str:
    walls = ", ".join(_plane_ref(wall) for wall in cell.vertical_walls) or "none"
    return (
        f"cell {cell.id}: "
        f"floor={_plane_ref(cell.floor)} "
        f"ceiling={_plane_ref(cell.ceiling)} "
        f"walls=[{walls}] "
        f"z=({_z_lower(cell.lower_z)}, {_z_upper(cell.upper_z)})"
    )


def format_result(result: SweepResult, *, seed: int | None = None) -> str:
    lines = [
        f"planes: {len(result.planes)}",
        f"events: {len(result.events)}",
        f"intervals: {len(result.intervals)}",
        f"cells: {len(result.cells)}",
    ]
    if seed is not None:
        lines.append(f"seed: {seed}")
    lines.extend(["", "input planes:"])
    for plane in result.planes:
        lines.append(f"  P{plane.id}: {format_plane(plane)}")
    lines.append("")
    lines.append("3D cells:")
    for cell in result.cells:
        lines.append(f"  {format_cell(cell)}")
    return "\n".join(lines) + "\n"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="vd3d",
        description="Compute the vertical decomposition of an arrangement of planes.",
    )
    parser.add_argument(
        "planes_file",
        nargs="?",
        type=Path,
        help="text file of planes (each line: `a b c d` or `id a b c d`)",
    )
    parser.add_argument(
        "--random-planes",
        type=int,
        metavar="N",
        help="ignore the file and use N random general-position planes",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        metavar="S",
        help="RNG seed for --random-planes (default: 0)",
    )
    parser.add_argument(
        "--write-planes",
        type=Path,
        metavar="PATH",
        help="write the input planes (useful with --random-planes) to PATH",
    )
    parser.add_argument(
        "--reference",
        action="store_true",
        help="recompute the 2D VD at every event instead of incremental updates",
    )
    parser.add_argument(
        "--gui",
        action="store_true",
        help="open an interactive 3D view (planes grey, lines dark, one cell red)",
    )
    parser.add_argument(
        "--show-sweep",
        action="store_true",
        help="open before/after sweep-plane views for each event",
    )
    return parser


def resolve_planes(args: argparse.Namespace) -> list[Plane]:
    if args.random_planes is not None:
        if args.random_planes < 1:
            raise ValueError("--random-planes N requires N >= 1")
        rng = random.Random(args.seed)
        return random_general_position_planes(rng, args.random_planes)
    if args.planes_file is None:
        raise ValueError("provide a planes file, or use --random-planes N")
    if not args.planes_file.is_file():
        raise ValueError(f"planes file not found: {args.planes_file}")
    return load_planes(args.planes_file)


def run(planes: Sequence[Plane], *, incremental: bool = True) -> SweepResult:
    return vertical_decomposition_3d(planes, incremental=incremental)


def warn_if_not_general_position(result: SweepResult, *, file=None) -> list[str]:
    """Print a warning for every ``z`` shared by two or more events.

    Returns the warning lines (empty if the instance is in general position).
    """
    if file is None:
        file = sys.stderr
    warnings: list[str] = []
    for group in group_events_by_z(result.events):
        if len(group.events) < 2:
            continue
        kinds = ", ".join(event.type.name for event in group.events)
        msg = (
            f"warning: {len(group.events)} events share z={group.z} ({kinds}); "
            "planes are not in general position"
        )
        warnings.append(msg)
        print(msg, file=file)
    return warnings


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        planes = resolve_planes(args)
        if args.write_planes is not None:
            dump_planes(planes, args.write_planes)
        result = run(planes, incremental=not args.reference)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    seed = args.seed if args.random_planes is not None else None
    warn_if_not_general_position(result)
    sys.stdout.write(format_result(result, seed=seed))
    if args.gui:
        try:
            from vd3d.viz.gui import show_decomposition
        except ImportError as exc:
            print(f"error: matplotlib is required for --gui ({exc})", file=sys.stderr)
            return 1
        try:
            show_decomposition(result, seed=seed)
        except RuntimeError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 1
    if args.show_sweep:
        try:
            from vd3d.viz.sweep_gui import show_sweep
        except ImportError as exc:
            print(
                f"error: matplotlib is required for --show-sweep ({exc})",
                file=sys.stderr,
            )
            return 1
        try:
            show_sweep(result, seed=seed)
        except (RuntimeError, ValueError) as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
