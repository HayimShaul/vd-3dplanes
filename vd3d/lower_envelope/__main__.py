"""CLI for the lower-envelope pipeline.

Examples::

    python -m vd3d.lower_envelope examples/three_planes.txt --gui
    python -m vd3d.lower_envelope examples/three_planes.txt --exact --gui
    python -m vd3d.lower_envelope examples/three_planes.txt --vd --gui
    python -m vd3d.lower_envelope examples/three_planes.txt --exact --vd --gui
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

from vd3d.geometry.plane import Plane
from vd3d.geometry.sampling import random_general_position_planes
from vd3d.io import dump_planes, load_planes
from vd3d.lower_envelope.envelope import (
    compare_envelopes_float_exact,
    format_compare,
    format_envelope,
    lower_envelope,
)
from vd3d.lower_envelope.envelope_float import (
    format_float_envelope,
    lower_envelope_float,
)
from vd3d.lower_envelope.vd import (
    compare_hypograph_float_exact,
    format_hypograph_compare,
    format_hypograph_vd,
    hypograph_vd,
)
from vd3d.lower_envelope.vd_float import (
    format_float_hypograph_vd,
    hypograph_vd_float,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="vd3d.lower_envelope",
        description="Lower envelope of planes and hypograph VD.",
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
        help="write the input planes to PATH",
    )
    parser.add_argument(
        "--exact",
        action="store_true",
        help="use exact Fraction arithmetic (and compare to float)",
    )
    parser.add_argument(
        "--vd",
        action="store_true",
        help="compute the hypograph vertical decomposition",
    )
    parser.add_argument(
        "--gui",
        action="store_true",
        help="open an interactive 3D view",
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


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        planes = resolve_planes(args)
        if args.write_planes is not None:
            dump_planes(planes, args.write_planes)

        float_env = lower_envelope_float(planes)
        exact_env = lower_envelope(planes) if args.exact else None
        env_cmp = (
            compare_envelopes_float_exact(float_env, exact_env)
            if exact_env is not None
            else None
        )

        float_vd = hypograph_vd_float(float_env) if args.vd else None
        exact_vd = hypograph_vd(exact_env) if (args.vd and exact_env is not None) else None
        vd_cmp = (
            compare_hypograph_float_exact(float_vd, exact_vd)
            if float_vd is not None and exact_vd is not None
            else None
        )
    except (OSError, ValueError, RuntimeError, ImportError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    seed = args.seed if args.random_planes is not None else None

    if exact_vd is not None:
        sys.stdout.write(format_hypograph_vd(exact_vd, seed=seed))
        assert vd_cmp is not None
        sys.stdout.write(format_hypograph_compare(vd_cmp))
    elif float_vd is not None:
        sys.stdout.write(format_float_hypograph_vd(float_vd, seed=seed))
    elif exact_env is not None:
        sys.stdout.write(format_envelope(exact_env, seed=seed))
        assert env_cmp is not None
        sys.stdout.write(format_compare(env_cmp))
    else:
        sys.stdout.write(format_float_envelope(float_env, seed=seed))

    if args.gui:
        try:
            if exact_vd is not None:
                from vd3d.viz.hypograph_gui import show_hypograph_vd

                show_hypograph_vd(
                    exact_vd, seed=seed, exact=True, compare=vd_cmp
                )
            elif float_vd is not None:
                from vd3d.viz.hypograph_gui import show_hypograph_vd

                show_hypograph_vd(float_vd, seed=seed, exact=False)
            elif exact_env is not None:
                from vd3d.viz.lower_envelope_gui import show_lower_envelope

                show_lower_envelope(
                    exact_env,
                    seed=seed,
                    float_overlay=float_env,
                    compare=env_cmp,
                    exact=True,
                )
            else:
                from vd3d.viz.lower_envelope_gui import show_lower_envelope

                show_lower_envelope(float_env, seed=seed, exact=False)
        except ImportError as exc:
            print(f"error: matplotlib is required for --gui ({exc})", file=sys.stderr)
            return 1
        except RuntimeError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
