#!/usr/bin/env python3
"""Lint a SIMCOM practice report: asset paths, labels, placeholders, coverage."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

LABEL_CMD_RE = re.compile(r"\\label\{([^}]+)\}")
LABEL_OPT_RE = re.compile(r"label\s*=\s*\{?([^,}\]]+)\}?")
REF_RE = re.compile(r"\\(?:eq|page|auto)?ref\{([^}]+)\}")
PH_RE = re.compile(r"\\ph\{([^}]*)\}")
TODO_RE = re.compile(r"%\s*TODO\b")
EX_SECTION_RE = re.compile(r"\\section\{Exercise\s+(\d+)\}", re.IGNORECASE)
NATURAL_EX_RE = re.compile(r"^ex(\d+)$", re.IGNORECASE)
INCLUDE_CMD_RE = re.compile(r"\\(includegraphics|lstinputlisting)(?![A-Za-z])")
CAPTION_RE = re.compile(r"\\caption\{")
GENERIC_CAPTION_RE = re.compile(r"^\s*(Output|Screenshot) of exercise\s+\d+\s*\.?\s*$", re.IGNORECASE)
FIGURE_ENV_RE = re.compile(r"\\begin\{figure\}(.*?)\\end\{figure\}", re.DOTALL)
SKIP_DIRS = {"report", "data", "material"}
NUM_RE = re.compile(r"(\d+)")


def strip_comments(text: str) -> str:
    lines = []
    for line in text.splitlines():
        out = []
        i = 0
        while i < len(line):
            ch = line[i]
            if ch == "%" and (i == 0 or line[i - 1] != "\\"):
                break
            out.append(ch)
            i += 1
        lines.append("".join(out))
    return "\n".join(lines)


def read_delimited(text: str, start: int, open_ch: str, close_ch: str) -> tuple[str, int]:
    """Read a possibly nested delimited group. `start` points at `open_ch`."""
    depth = 0
    i = start
    while i < len(text):
        ch = text[i]
        if ch == open_ch:
            depth += 1
        elif ch == close_ch:
            depth -= 1
            if depth == 0:
                return text[start + 1 : i], i + 1
        i += 1
    return text[start + 1 :], len(text)


def iter_includes(text: str) -> list[tuple[str, str]]:
    """Return (optional_args, path) for each includegraphics/lstinputlisting."""
    found: list[tuple[str, str]] = []
    for m in INCLUDE_CMD_RE.finditer(text):
        i = m.end()
        opt = ""
        if i < len(text) and text[i] == "[":
            opt, i = read_delimited(text, i, "[", "]")
        if i < len(text) and text[i] == "{":
            path, _ = read_delimited(text, i, "{", "}")
            found.append((opt, path.strip()))
    return found


def natural_key(name: str) -> list:
    return [int(tok) if tok.isdigit() else tok.lower() for tok in NUM_RE.split(name)]


def discover_exercises(practice_root: Path) -> list[str]:
    """Exercise directories: `ex*` if any, else any subdir holding a .f90."""
    if not practice_root.is_dir():
        return []
    subdirs = [
        p
        for p in practice_root.iterdir()
        if p.is_dir() and p.name.lower() not in SKIP_DIRS and not p.name.startswith(".")
    ]
    ex_dirs = [p.name for p in subdirs if NATURAL_EX_RE.match(p.name)]
    if ex_dirs:
        return sorted(ex_dirs, key=natural_key)
    code_dirs = [p.name for p in subdirs if any(p.glob("*.f90"))]
    if code_dirs:
        return sorted(code_dirs, key=natural_key)
    f90s = [p.stem for p in practice_root.glob("*.f90") if p.is_file()]
    return sorted(f90s, key=natural_key)


def mentioned_exercises(tex: str) -> set[str]:
    found: set[str] = set()
    for m in EX_SECTION_RE.finditer(tex):
        found.add(f"ex{m.group(1)}")
    for _, path in iter_includes(tex):
        parts = Path(path.replace("\\", "/")).parts
        if len(parts) >= 2 and parts[0] == "..":
            found.add(parts[1])
        found.add(Path(path).stem)
    return found


def iter_captions(text: str) -> list[str]:
    caps = []
    for m in CAPTION_RE.finditer(text):
        body, _ = read_delimited(text, m.end() - 1, "{", "}")
        caps.append(body)
    return caps


def unused_pngs(practice_root: Path, exercises: list[str], report_dir: Path, includes) -> list[str]:
    used = {(report_dir / rel.replace("\\", "/")).resolve() for _o, rel in includes}
    out = []
    for ex in exercises:
        ex_dir = practice_root / ex
        if not ex_dir.is_dir():
            continue
        for png in sorted(ex_dir.glob("*.png"), key=lambda p: natural_key(p.name)):
            if png.resolve() not in used:
                out.append(f"{ex}/{png.name}")
    return out


def collect_labels(stripped: str) -> set[str]:
    labels = set(LABEL_CMD_RE.findall(stripped))
    for opt, _path in iter_includes(stripped):
        labels.update(LABEL_OPT_RE.findall(opt))
    return labels


def lint(tex_path: Path) -> int:
    if not tex_path.is_file():
        print(f"ERROR: report not found: {tex_path}", file=sys.stderr)
        return 2

    raw = tex_path.read_text(encoding="utf-8")
    stripped = strip_comments(raw)
    report_dir = tex_path.resolve().parent
    practice_root = report_dir.parent
    includes = iter_includes(stripped)

    missing_assets: list[str] = []
    for _opt, rel in includes:
        rel_norm = rel.replace("\\", "/")
        target = (report_dir / rel_norm).resolve()
        if not target.is_file():
            missing_assets.append(rel_norm)

    labels = collect_labels(stripped)
    refs = REF_RE.findall(stripped)
    unresolved = sorted({r for r in refs if r not in labels})

    placeholders = PH_RE.findall(stripped)
    todos = []
    for i, line in enumerate(raw.splitlines(), 1):
        if TODO_RE.search(line):
            todos.append(f"L{i}: {line.strip()}")

    expected = discover_exercises(practice_root)
    mentioned = mentioned_exercises(stripped)
    missing_ex = [e for e in expected if e not in mentioned]

    generic_caps = [c for c in iter_captions(stripped) if GENERIC_CAPTION_RE.match(c)]
    ref_set = set(refs)
    unreferenced_figs = []
    for fig in FIGURE_ENV_RE.findall(stripped):
        for lab in LABEL_CMD_RE.findall(fig):
            if lab not in ref_set:
                unreferenced_figs.append(lab)
    unused = unused_pngs(practice_root, expected, report_dir, includes)

    print(f"Report: {tex_path}")
    print(f"Practice root: {practice_root}")
    print()

    if missing_assets:
        print(f"Missing assets ({len(missing_assets)}):")
        for rel in missing_assets:
            print(f"  - {rel}")
    else:
        print("Assets: OK")

    print()
    if unresolved:
        print(f"Unresolved refs ({len(unresolved)}):")
        for r in unresolved:
            print(f"  - {r}")
    else:
        print("Refs: OK")

    print()
    if placeholders:
        print(f"Placeholders ({len(placeholders)}):")
        for p in placeholders:
            print(f"  - [{p}]")
    else:
        print("Placeholders: none")

    print()
    if todos:
        print(f"TODO markers ({len(todos)}):")
        for t in todos:
            print(f"  - {t}")
    else:
        print("TODO markers: none")

    print()
    if missing_ex:
        print(f"Exercises missing from report ({len(missing_ex)}):")
        for e in missing_ex:
            print(f"  - {e}")
    else:
        print("Coverage: all discovered exercises are mentioned")

    print()
    print("Style warnings (non-fatal):")
    if not (generic_caps or unreferenced_figs or unused):
        print("  none")
    for c in generic_caps:
        print(f"  - generic caption, describe the panels instead: \"{c.strip()}\"")
    for lab in unreferenced_figs:
        print(f"  - figure never cited with \\ref: {lab}")
    for png in unused:
        print(f"  - plot not included in report: {png}")

    return 1 if missing_assets else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tex", type=Path, help="Path to P*_report.tex")
    args = parser.parse_args()
    return lint(args.tex)


if __name__ == "__main__":
    sys.exit(main())
