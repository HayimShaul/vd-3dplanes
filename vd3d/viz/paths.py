"""Filesystem layout for generated review figures."""

from __future__ import annotations

from pathlib import Path


def repo_root() -> Path:
    """Repository root (parent of the ``vd3d`` package)."""
    return Path(__file__).resolve().parents[2]


def visual_root() -> Path:
    return repo_root() / "artifacts" / "visual"


def step_dir(step: int | str) -> Path:
    """Directory for one review step: ``artifacts/visual/step_00``."""
    return visual_root() / f"step_{_step_tag(step)}"


def gallery_index() -> Path:
    return visual_root() / "index.html"


def _step_tag(step: int | str) -> str:
    if isinstance(step, str):
        text = step.strip()
        if text.startswith("step_"):
            text = text[len("step_") :]
        return f"{int(text):02d}"
    return f"{int(step):02d}"
