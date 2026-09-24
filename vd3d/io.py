"""Plain-text plane files for the CLI.

Format (one plane per line)::

    # comments and blank lines are ignored
    a b c d          # id assigned 1, 2, 3, …
    id a b c d       # explicit id

Coefficients are exact ``Fraction`` values (``int`` or ratio strings like ``3/2``).
Each line is the oriented plane ``a*x + b*y + c*z + d = 0``.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from vd3d.geometry.plane import Plane
from vd3d.geometry.scalar import Scalar, as_scalar


def load_planes(path: str | Path) -> list[Plane]:
    """Read planes from a text file. Raises ``ValueError`` on bad lines."""
    text = Path(path).read_text(encoding="utf-8")
    planes: list[Plane] = []
    next_id = 1
    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        parts = line.split()
        try:
            if len(parts) == 4:
                plane_id = next_id
                next_id += 1
                coeffs = parts
            elif len(parts) == 5:
                plane_id = int(parts[0])
                coeffs = parts[1:]
            else:
                raise ValueError(
                    f"expected 4 coefficients or `id a b c d`, got {len(parts)} fields"
                )
            a, b, c, d = (as_scalar(token) for token in coeffs)
            planes.append(Plane(id=plane_id, a=a, b=b, c=c, d=d))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{path}:{lineno}: {exc}") from exc
    if not planes:
        raise ValueError(f"{path}: no planes found")
    return planes


def format_plane_record(plane: Plane) -> str:
    """``id a b c d`` with exact Fraction coefficients (file format)."""
    return f"{plane.id} {plane.a} {plane.b} {plane.c} {plane.d}"


def format_plane(plane: Plane) -> str:
    """Human-readable equation, e.g. ``2x+2y-4z=0`` or ``x+y+z-5=0``."""
    return format_plane_equation(plane.a, plane.b, plane.c, plane.d)


def format_plane_equation(
    a: int | Scalar | str,
    b: int | Scalar | str,
    c: int | Scalar | str,
    d: int | Scalar | str,
) -> str:
    """Join nonzero terms with signed ``+`` / ``-`` into ``…=0``."""
    names = ("x", "y", "z", "")
    parts: list[str] = []
    for coeff, name in zip((a, b, c, d), names, strict=True):
        value = as_scalar(coeff)
        if value == 0:
            continue
        if name:
            if value == 1:
                term = name
            elif value == -1:
                term = f"-{name}"
            else:
                term = f"{value}{name}"
        else:
            term = str(value)
        parts.append(term)
    if not parts:
        raise ValueError("zero plane")
    body = parts[0]
    for part in parts[1:]:
        if part.startswith("-"):
            body += part
        else:
            body += f"+{part}"
    return f"{body}=0"


def dump_planes(planes: Sequence[Plane], path: str | Path) -> None:
    """Write planes as ``id a b c d`` lines."""
    body = "\n".join(format_plane_record(plane) for plane in planes) + "\n"
    Path(path).write_text(body, encoding="utf-8")
