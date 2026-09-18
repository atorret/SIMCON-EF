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
EX_FILE_RE = re.compile(r"ex(\d+)", re.IGNORECASE)
NATURAL_EX_RE = re.compile(r"^ex(\d+)$", re.IGNORECASE)
INCLUDE_CMD_RE = re.compile(r"\\(includegraphics|lstinputlisting)(?![A-Za-z])")


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


def natural_ex_key(name: str) -> tuple[int, str]:
    m = NATURAL_EX_RE.match(name)
    return (int(m.group(1)), name) if m else (10**9, name)


def discover_exercises(practice_root: Path) -> list[str]:
    if not practice_root.is_dir():
        return []
    dirs = [
        p.name
        for p in practice_root.iterdir()
        if p.is_dir() and NATURAL_EX_RE.match(p.name)
    ]
    if dirs:
        return sorted(dirs, key=natural_ex_key)
    f90s = [p.stem for p in practice_root.glob("*.f90") if p.is_file()]
    return sorted(f90s, key=natural_ex_key)


def mentioned_exercises(tex: str) -> set[str]:
    found: set[str] = set()
    for m in EX_SECTION_RE.finditer(tex):
        found.add(f"ex{m.group(1)}")
    for _, path in iter_includes(tex):
        stem = Path(path.replace("\\", "/")).stem
        em = EX_FILE_RE.match(stem)
        if em:
            found.add(f"ex{em.group(1)}")
    return found


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

    return 1 if missing_assets else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tex", type=Path, help="Path to P*_report.tex")
    args = parser.parse_args()
    return lint(args.tex)


if __name__ == "__main__":
    sys.exit(main())
