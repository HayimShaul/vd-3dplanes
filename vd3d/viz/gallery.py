"""Build ``artifacts/visual/index.html`` from recorded figures and captions.

Usage::

    python -m vd3d.viz.gallery
    python -m vd3d.viz.gallery --step 0
"""

from __future__ import annotations

import argparse
import html
import re
from pathlib import Path

from vd3d.viz.paths import gallery_index, visual_root

_STEP_DIR_RE = re.compile(r"^step_(\d+)$")


def build_gallery(*, step: int | None = None) -> Path:
    """Write the HTML index and return its path."""
    root = visual_root()
    root.mkdir(parents=True, exist_ok=True)
    sections = _collect_sections(root, step=step)
    index = gallery_index()
    index.write_text(_render_html(sections, filtered_step=step), encoding="utf-8")
    return index


def _collect_sections(root: Path, *, step: int | None) -> list[tuple[str, list[dict[str, str]]]]:
    dirs = sorted(
        path for path in root.iterdir() if path.is_dir() and _STEP_DIR_RE.match(path.name)
    )
    if step is not None:
        tag = f"{int(step):02d}"
        dirs = [path for path in dirs if path.name == f"step_{tag}"]

    sections: list[tuple[str, list[dict[str, str]]]] = []
    for directory in dirs:
        items: list[dict[str, str]] = []
        stems = sorted({path.stem for path in directory.iterdir() if path.suffix in {".png", ".md"}})
        for stem in stems:
            png = directory / f"{stem}.png"
            md = directory / f"{stem}.md"
            items.append(
                {
                    "name": stem,
                    "png": f"{directory.name}/{stem}.png" if png.exists() else "",
                    "caption": md.read_text(encoding="utf-8") if md.exists() else "",
                    "missing_image": "" if png.exists() else "missing image",
                    "missing_caption": "" if md.exists() else "missing caption",
                }
            )
        sections.append((directory.name, items))
    return sections


def _render_html(sections: list[tuple[str, list[dict[str, str]]]], *, filtered_step: int | None) -> str:
    filter_note = (
        f"<p class='filter'>Showing only <code>step_{int(filtered_step):02d}</code>.</p>"
        if filtered_step is not None
        else ""
    )
    body = "\n".join(_render_section(name, items) for name, items in sections)
    if not body:
        body = (
            "<p class='empty'>No figures yet. Run the visual tests first, "
            "for example <code>pytest tests/visual/test_gallery_smoke.py</code>.</p>"
        )
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>VD3D visual review</title>
  <style>
    body {{ font-family: sans-serif; max-width: 960px; margin: 2rem auto; padding: 0 1rem; line-height: 1.45; }}
    h1 {{ margin-bottom: 0.25rem; }}
    .review {{ background: #fff6d6; border: 1px solid #d4b106; padding: 1rem 1.25rem; }}
    .review ol {{ margin: 0.5rem 0 0; }}
    .filter, .empty {{ color: #444; }}
    section {{ margin: 2.5rem 0; }}
    article {{ border-top: 1px solid #ddd; padding: 1.25rem 0; }}
    img {{ max-width: 100%; height: auto; border: 1px solid #ccc; background: #fff; }}
    .caption {{ white-space: pre-wrap; }}
    .warn {{ color: #a40000; font-weight: bold; }}
    code {{ background: #f4f4f4; padding: 0.1em 0.3em; }}
  </style>
</head>
<body>
  <h1>VD3D visual review</h1>
  <div class="review">
    <strong>How to review a step</strong>
    <ol>
      <li>Run that step's unit and visual tests (see the README).</li>
      <li>Open this page and compare each figure to its caption.</li>
      <li>The caption states what you <em>must</em> see if the step is correct.</li>
      <li>If a figure and its caption disagree, the step is not done. Do not start the next step.</li>
    </ol>
  </div>
  {filter_note}
  {body}
</body>
</html>
"""


def _render_section(name: str, items: list[dict[str, str]]) -> str:
    articles = "\n".join(_render_item(item) for item in items) or "<p class='empty'>No figures in this step.</p>"
    return f"<section>\n  <h2>{html.escape(name)}</h2>\n  {articles}\n</section>"


def _render_item(item: dict[str, str]) -> str:
    warnings = []
    if item["missing_image"]:
        warnings.append(item["missing_image"])
    if item["missing_caption"]:
        warnings.append(item["missing_caption"])
    warn_html = "".join(f"<p class='warn'>{html.escape(w)}</p>" for w in warnings)
    img_html = (
        f"<p><img src='{html.escape(item['png'])}' alt='{html.escape(item['name'])}'></p>"
        if item["png"]
        else ""
    )
    caption = _basic_markdown(item["caption"]) if item["caption"] else ""
    return (
        f"<article>\n"
        f"  <h3>{html.escape(item['name'])}</h3>\n"
        f"  {warn_html}\n"
        f"  {img_html}\n"
        f"  <div class='caption'>{caption}</div>\n"
        f"</article>"
    )


def _basic_markdown(text: str) -> str:
    """Escape HTML, then apply a tiny subset of markdown used in captions."""
    escaped = html.escape(text.strip())
    escaped = re.sub(r"^# (.+)$", r"<strong>\1</strong>", escaped, count=1, flags=re.MULTILINE)
    escaped = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", escaped)
    return escaped


def main(argv: list[str] | None = None) -> Path:
    parser = argparse.ArgumentParser(description="Build the visual review gallery")
    parser.add_argument("--step", type=int, default=None, help="only include this step number")
    args = parser.parse_args(argv)
    path = build_gallery(step=args.step)
    print(path)
    return path


if __name__ == "__main__":
    main()
