"""Write a PNG and a caption that a human can review together."""

from __future__ import annotations

from pathlib import Path

from vd3d.viz.paths import step_dir


def record_figure(fig, *, step: int, name: str, caption: str, dpi: int = 120) -> Path:
    """Save ``fig`` as ``artifacts/visual/step_XX/<name>.png`` plus a ``.md`` caption.

    ``caption`` must state what a human must see if the figure is correct.
    Returns the PNG path.
    """
    if not name or "/" in name or "\\" in name:
        raise ValueError(f"invalid figure name: {name!r}")

    out_dir = step_dir(step)
    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / f"{name}.png"
    md_path = out_dir / f"{name}.md"

    fig.savefig(png_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    md_path.write_text(_format_caption(name, caption), encoding="utf-8")
    return png_path


def _format_caption(name: str, caption: str) -> str:
    body = caption.strip()
    if body.startswith("#"):
        return body + "\n"
    return f"# {name}\n\n{body}\n"
