#!/usr/bin/env python3
"""
make_plots.py - turn every data file under ``Plots/`` into tidy figures and
one quick-read PDF report.

Layout it understands (discovered automatically, nothing is hard-coded):

    Plots/
      Step3/                 <- every first-level folder is a *section*
        Energy.dat           <- loose files: one panel each (overview page)
        Temperature.dat
        G(r)/                <- every sub-folder is a *group* (grid page)
          Ow_Hw.dat
          ...
      Step4/  Question1/ ... <- new folders are picked up on the next run
      _output/               <- generated (ignored by the discovery)

Data files: *.dat, *.txt, *.csv, *.xvg. Lines starting with # @ ; ! % are
comments; an optional first non-numeric line is used as column names.
Column 0 is x and every further column is a series (a single column is
plotted against its index).

How a file is drawn is decided by RULES below (first match wins). Add a line
there to teach the script a new kind of file; unmatched files still get a
clean generic plot.

Usage:
    python make_plots.py                      # everything
    python make_plots.py --only Step4         # just some sections
    python make_plots.py --no-individual      # only the PDF report
"""
from __future__ import annotations

import argparse
import datetime as dt
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.ticker import FuncFormatter, MaxNLocator  # noqa: E402

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
DATA_EXT = {".dat", ".txt", ".csv", ".xvg"}
COMMENT_PREFIXES = ("#", "@", ";", "!", "%")
OUTPUT_DIRNAME = "_output"
REPORT_NAME = "Plots_report.pdf"
A4_LANDSCAPE = (11.69, 8.27)
MAX_PANELS_PER_PAGE = 12

# Palette (light, print-friendly). Categorical slots in fixed order.
SURFACE = "#ffffff"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8984"
AXIS = "#b9b8b3"
GRID = "#ebeae6"
ALERT = "#c4352f"
BLUE, BLUE_LIGHT, BLUE_DARK = "#2a78d6", "#9ec5f4", "#1c5cab"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
          "#e87ba4", "#008300", "#6250d6", "#e34948"]


@dataclass(frozen=True)
class Rule:
    """How to draw files whose name/path matches ``pattern`` (case-insensitive).

    on="name": pattern is searched in the file stem (e.g. "Energy").
    on="path": pattern is searched in the path relative to the section
               (e.g. "G(r)/Ow_Hw.dat"), useful to match by folder.
    kind: "timeseries" | "rdf" | "cumulative" | "line"
    """
    pattern: str
    kind: str = "line"
    xlabel: str | None = None
    ylabel: str | None = None
    title: str | None = None
    series: tuple[str, ...] = ()
    on: str = "name"


RULES: list[Rule] = [
    Rule(r"(^|/)g\(r\)/[^/]+$", "rdf", "r (Å)", "g(r)", on="path"),
    Rule(r"(^|/)n\(r\)/[^/]+$", "cumulative", "r (Å)", "N(r)", on="path"),
    Rule(r"energy", "timeseries", "Timestep", "Energy (kcal/mol)", "Energy"),
    Rule(r"temp", "timeseries", "Timestep", "Temperature (K)", "Temperature"),
    Rule(r"pres", "timeseries", "Timestep", "Pressure (bar)", "Pressure",
         series=("instantaneous", "PRESSAVG")),
    Rule(r"hydrogen|h_?bond", "timeseries", "Frame", "Number of H-bonds",
         "Hydrogen bonds"),
]
DEFAULT_RULE = Rule(r".*", "line")

# Readable titles for group folders (matched on the folder path in a section).
GROUP_TITLES: list[tuple[str, str]] = [
    (r"(^|/)g\(r\)$", "Radial distribution functions g(r)"),
    (r"(^|/)n\(r\)$", "Running coordination numbers N(r)"),
]

# Atom-name tokens rendered with subscripts in pair titles ("Ow_Hw").
TOKEN_TEX = {
    "ow": r"O$_\mathrm{w}$", "hw": r"H$_\mathrm{w}$",
    "oh2": r"O$_\mathrm{w}$", "h1": r"H$_1$", "h2": r"H$_2$",
}


def apply_style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Segoe UI", "Helvetica Neue", "Helvetica", "Arial",
                            "DejaVu Sans"],
        "mathtext.fontset": "custom",
        "mathtext.rm": "sans",
        "mathtext.it": "sans:italic",
        "mathtext.bf": "sans:bold",
        "font.size": 9,
        "text.color": INK,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "axes.edgecolor": AXIS,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "axes.labelcolor": INK2,
        "axes.labelsize": 8.5,
        "axes.titlesize": 10.5,
        "axes.titleweight": "semibold",
        "axes.titlelocation": "left",
        "axes.titlecolor": INK,
        "axes.formatter.useoffset": False,
        "axes.formatter.use_mathtext": True,
        "axes.unicode_minus": True,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "legend.fontsize": 7.5,
        "legend.frameon": False,
        "lines.solid_capstyle": "round",
        "figure.dpi": 110,
        "savefig.dpi": 200,
        "pdf.fonttype": 42,
    })


# --------------------------------------------------------------------------
# Data model
# --------------------------------------------------------------------------
@dataclass
class DataFile:
    path: Path
    section: str
    rel: str                     # path relative to the section, posix
    rule: Rule
    x: np.ndarray | None = None
    ys: list[np.ndarray] = field(default_factory=list)
    names: list[str] = field(default_factory=list)
    xname: str | None = None     # from a header line, if any
    status: str = "ok"          # ok | zero | constant | empty | error
    message: str = ""
    stats: dict = field(default_factory=dict)

    @property
    def stem(self) -> str:
        return self.path.stem

    @property
    def kind(self) -> str:
        return self.rule.kind

    @property
    def title(self) -> str:
        return self.rule.title or pretty_name(self.stem)

    @property
    def flagged(self) -> bool:
        return self.status != "ok"


@dataclass
class Group:
    rel: str                     # "." for loose files of the section
    files: list[DataFile]

    @property
    def title(self) -> str:
        if self.rel == ".":
            return "Overview"
        for pat, title in GROUP_TITLES:
            if re.search(pat, self.rel, re.I):
                return title
        return self.rel


@dataclass
class Section:
    name: str
    groups: list[Group]

    @property
    def files(self) -> list[DataFile]:
        return [f for g in self.groups for f in g.files]


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def natural_key(s: str):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


def pretty_name(stem: str) -> str:
    tokens = stem.split("_")
    if len(tokens) == 2 and all(0 < len(t) <= 4 and t.isalnum() for t in tokens):
        return " – ".join(TOKEN_TEX.get(t.lower(), t) for t in tokens)
    return stem.replace("_", " ")


def fmt_num(v: float, sig: int = 4) -> str:
    if v is None or not np.isfinite(v):
        return "–"
    if v == 0:
        return "0"
    mag = math.floor(math.log10(abs(v)))
    decimals = max(0, sig - 1 - mag)
    return f"{v:,.{decimals}f}".replace(",", " ").replace("-", "−")


def match_rule(rel: str, stem: str) -> Rule:
    for rule in RULES:
        target = rel if rule.on == "path" else stem
        if re.search(rule.pattern, target, re.I):
            return rule
    return DEFAULT_RULE


def read_table(path: Path) -> tuple[np.ndarray, list[str] | None]:
    rows: list[list[float]] = []
    header: list[str] | None = None
    for raw in path.read_text(errors="replace").splitlines():
        s = raw.strip()
        if not s or s.startswith(COMMENT_PREFIXES):
            continue
        tokens = [t for t in re.split(r"[,\s]+", s) if t]
        try:
            rows.append([float(t) for t in tokens])
        except ValueError:
            if not rows and header is None:
                header = tokens
                continue
            raise ValueError(f"non-numeric line: {s[:40]!r}")
    if not rows:
        return np.empty((0, 0)), header
    width = len(rows[0])
    if any(len(r) != width for r in rows):
        raise ValueError("rows have different numbers of columns")
    return np.asarray(rows, dtype=float), header


def load_file(path: Path, section_dir: Path, section: str) -> DataFile:
    rel = path.relative_to(section_dir).as_posix()
    f = DataFile(path=path, section=section, rel=rel, rule=match_rule(rel, path.stem))
    try:
        arr, header = read_table(path)
    except Exception as exc:  # one bad file must not stop the run
        f.status, f.message = "error", str(exc)
        return f
    if arr.size == 0:
        f.status, f.message = "empty", "no numeric rows"
        return f
    if arr.shape[1] == 1:
        f.x, f.ys = np.arange(arr.shape[0], dtype=float), [arr[:, 0]]
    else:
        f.x, f.ys = arr[:, 0], [arr[:, i] for i in range(1, arr.shape[1])]

    n = len(f.ys)
    if len(f.rule.series) == n:
        f.names = list(f.rule.series)
    elif header and len(header) == arr.shape[1]:
        f.names = header[1:] if arr.shape[1] > 1 else header
        if arr.shape[1] > 1:
            f.xname = header[0]
    elif n == 1:
        f.names = [f.rule.ylabel or "value"]
    else:
        f.names = [f"column {i + 2}" for i in range(n)]

    if all(np.all(y == 0) for y in f.ys):
        f.status, f.message = "zero", f"all {len(f.x)} values are 0"
    elif all(np.ptp(y) == 0 for y in f.ys):
        f.status, f.message = "constant", "values are constant"
    return f


def discover(root: Path, out_dir: Path) -> list[Section]:
    def skipped(name: str) -> bool:
        return name.startswith((".", "_"))

    sections: list[Section] = []
    candidates: list[tuple[str, Path, bool]] = []
    loose_root = [p for p in root.iterdir()
                  if p.is_file() and p.suffix.lower() in DATA_EXT]
    if loose_root:
        candidates.append((root.name, root, False))
    for d in sorted((p for p in root.iterdir() if p.is_dir()),
                    key=lambda p: natural_key(p.name)):
        if skipped(d.name) or d.resolve() == out_dir.resolve():
            continue
        candidates.append((d.name, d, True))

    for name, d, recursive in candidates:
        files = d.rglob("*") if recursive else d.iterdir()
        by_dir: dict[str, list[Path]] = {}
        for p in files:
            if not p.is_file() or p.suffix.lower() not in DATA_EXT:
                continue
            rel_parts = p.relative_to(d).parts[:-1]
            if any(skipped(part) for part in rel_parts):
                continue
            by_dir.setdefault(p.parent.relative_to(d).as_posix(), []).append(p)
        if not by_dir:
            continue
        groups = []
        for rel in sorted(by_dir, key=lambda r: (r != ".", natural_key(r))):
            paths = sorted(by_dir[rel], key=lambda p: natural_key(p.name))
            groups.append(Group(rel, [load_file(p, d, name) for p in paths]))
        sections.append(Section(name, groups))
    return sections


# --------------------------------------------------------------------------
# Analysis
# --------------------------------------------------------------------------
def smooth3(y: np.ndarray) -> np.ndarray:
    if len(y) < 3:
        return y.copy()
    s = np.convolve(y, np.ones(3) / 3, mode="same")
    s[0], s[-1] = y[0], y[-1]
    return s


def running_mean(x: np.ndarray, y: np.ndarray, w: int):
    if w < 2 or len(y) < w:
        return x, y
    ym = np.convolve(y, np.ones(w) / w, mode="valid")
    off = w // 2
    return x[off:off + len(ym)], ym


def rdf_features(r: np.ndarray, g: np.ndarray) -> dict:
    s = smooth3(g)
    gmax = s.max()
    peak = None
    for i in range(1, len(s) - 1):
        if s[i] >= s[i - 1] and s[i] > s[i + 1] and s[i] >= 0.5 * gmax:
            peak = i
            break
    if peak is None:
        peak = int(np.argmax(s))
    # smoothing locates the feature; the raw data gives its exact position
    lo = max(0, peak - 2)
    peak_raw = lo + int(np.argmax(g[lo:peak + 3]))
    out = {"r_peak": float(r[peak_raw]), "g_peak": float(g[peak_raw])}
    for j in range(peak + 1, len(s) - 1):
        if s[j] <= s[j - 1] and s[j] < s[j + 1] and s[j] < 0.9 * s[peak]:
            lo = max(peak_raw + 1, j - 2)
            j_raw = lo + int(np.argmin(g[lo:j + 3]))
            out["r_min"], out["g_min"] = float(r[j_raw]), float(g[j_raw])
            break
    return out


def analyse(sections: list[Section], equil_frac: float) -> None:
    for sec in sections:
        for f in sec.files:
            if f.status in ("error", "empty", "zero"):
                continue
            if f.kind == "timeseries":
                x0, x1 = f.x.min(), f.x.max()
                start = x0 + equil_frac * (x1 - x0)
                mask = f.x >= start
                if mask.sum() < 2:
                    mask = np.ones_like(f.x, dtype=bool)
                f.stats["equil_end"] = float(start)
                f.stats["mean"] = [float(np.mean(y[mask])) for y in f.ys]
                f.stats["std"] = [float(np.std(y[mask], ddof=1)) for y in f.ys]
            elif f.kind == "rdf":
                f.stats.update(rdf_features(f.x, f.ys[0]))

        # Coordination numbers: N(r) evaluated at the first minimum of the
        # g(r) file with the same name (any folder of the same section).
        rdfs = [f for f in sec.files if f.kind == "rdf" and "r_min" in f.stats]
        for f in sec.files:
            if f.kind != "cumulative" or f.status in ("error", "empty", "zero"):
                continue
            partners = [g for g in rdfs if g.stem.lower() == f.stem.lower()]
            if not partners:
                continue
            same_parent = [g for g in partners
                           if g.path.parent.parent == f.path.parent.parent]
            partner = (same_parent or partners)[0]
            r_min = partner.stats["r_min"]
            f.stats["r_min"] = r_min
            f.stats["n_min"] = float(np.interp(r_min, f.x, f.ys[0]))


def subtitle(f: DataFile) -> str:
    if f.status in ("error", "empty", "zero"):
        return ""
    st = f.stats
    if f.kind == "timeseries" and "mean" in st:
        unit = ""
        if f.rule.ylabel and "(" in f.rule.ylabel:
            unit = " " + f.rule.ylabel[f.rule.ylabel.rfind("(") + 1:-1]
        parts = [f"{fmt_num(m)} ± {fmt_num(s, 2)}{unit}"
                 for m, s in zip(st["mean"], st["std"])]
        if len(parts) == 1:
            return "mean " + parts[0]
        return "   ".join(f"{n}: {p}" for n, p in zip(f.names, parts))
    if f.kind == "rdf" and "r_peak" in st:
        s = f"1st max {fmt_num(st['r_peak'], 3)} Å (g = {fmt_num(st['g_peak'], 3)})"
        if "r_min" in st:
            s += f"   1st min {fmt_num(st['r_min'], 3)} Å"
        return s
    if f.kind == "cumulative":
        if "n_min" in st:
            return f"N at 1st min of g(r) ({fmt_num(st['r_min'], 3)} Å) = {fmt_num(st['n_min'], 3)}"
        return f"N({fmt_num(f.x[-1], 3)} Å) = {fmt_num(f.ys[0][-1], 3)}"
    return ""


# --------------------------------------------------------------------------
# Drawing
# --------------------------------------------------------------------------
def scaled_axis(ax, x: np.ndarray, label: str | None) -> str:
    """Show large integer-like x values as x / 10^k to keep ticks short."""
    label = label or ""
    span = float(np.nanmax(np.abs(x))) if len(x) else 0.0
    if span >= 1e4:
        k = 3 * int(math.log10(span) // 3)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v / 10**k:g}"))
        label = f"{label} (×10$^{{{k}}}$)".strip()
    return label


def put_title(ax, f: DataFile, compact: bool) -> None:
    sub = subtitle(f)
    ax.set_title(f.title, pad=16 if sub else 6,
                 fontsize=9.5 if compact else 11)
    if sub:
        ax.text(0, 1.015, sub, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=7 if compact else 8, color=INK2)


def draw_placeholder(ax, f: DataFile, compact: bool) -> None:
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    for side, sp in ax.spines.items():
        sp.set_visible(True)
        sp.set_color(AXIS)
        sp.set_linestyle((0, (3, 3)))
    head = {"zero": "No signal", "empty": "Empty file",
            "error": "Could not read"}.get(f.status, f.status)
    ax.text(0.5, 0.56, head, transform=ax.transAxes, ha="center", va="center",
            color=ALERT, fontsize=11 if compact else 14, fontweight="semibold")
    ax.text(0.5, 0.40, f.message, transform=ax.transAxes, ha="center",
            va="center", color=INK2, fontsize=7.5 if compact else 9, wrap=True)
    ax.set_title(f.title, pad=6, fontsize=9.5 if compact else 11)
    if f.rule.xlabel:
        ax.set_xlabel(f.rule.xlabel)
    if f.rule.ylabel:
        ax.set_ylabel(f.rule.ylabel)


def draw_timeseries(ax, f: DataFile, compact: bool) -> None:
    x = f.x
    if "equil_end" in f.stats and len(x) >= 10:
        ax.axvspan(x.min(), f.stats["equil_end"], color=INK2, alpha=0.06, lw=0)
        ax.text(x.min(), 1.0, " equil.", transform=ax.get_xaxis_transform(),
                ha="left", va="top", fontsize=6.5, color=MUTED)
    if len(f.ys) == 1:
        y = f.ys[0]
        w = max(3, int(round(len(y) / 15)) | 1)
        ax.plot(x, y, color=BLUE_LIGHT, lw=1.0, label="data")
        xm, ym = running_mean(x, y, w)
        ax.plot(xm, ym, color=BLUE_DARK, lw=2.0, label=f"running mean ({w} pts)")
        if "mean" in f.stats:
            ax.axhline(f.stats["mean"][0], color=INK2, lw=0.9, ls=(0, (4, 3)),
                       label="mean (after equil.)")
    else:
        for i, (y, name) in enumerate(zip(f.ys, f.names)):
            ax.plot(x, y, color=SERIES[i % len(SERIES)],
                    lw=1.1 if i == 0 else 1.8, alpha=0.85 if i == 0 else 1.0,
                    label=name, zorder=2 + i)
    ax.set_xlim(x.min(), x.max())
    # headroom at the top so the legend never sits on the data
    lo, hi = min(np.min(y) for y in f.ys), max(np.max(y) for y in f.ys)
    pad = (hi - lo) or abs(hi) or 1.0
    ax.set_ylim(lo - 0.05 * pad, hi + (0.38 if compact else 0.22) * pad)
    ax.set_xlabel(scaled_axis(ax, x, f.rule.xlabel))
    ax.set_ylabel(f.rule.ylabel or "")
    ax.legend(loc="upper right", ncol=1 if compact else 3, handlelength=1.6,
              borderaxespad=0.3)


def draw_rdf(ax, f: DataFile, compact: bool) -> None:
    r, g = f.x, f.ys[0]
    ax.fill_between(r, g, color=BLUE, alpha=0.12, lw=0)
    ax.plot(r, g, color=BLUE, lw=1.8)
    ax.axhline(1.0, color=INK2, lw=0.8, ls=":")
    st = f.stats
    if "r_min" in st:
        ax.axvline(st["r_min"], color=INK2, lw=0.8, ls=(0, (2, 2)), alpha=0.7)
    if "r_peak" in st:
        ax.plot(st["r_peak"], st["g_peak"], "o", ms=6, mfc=BLUE, mec=SURFACE,
                mew=1.5, zorder=5)
        ax.annotate(f"{fmt_num(st['r_peak'], 3)} Å", (st["r_peak"], st["g_peak"]),
                    xytext=(7, 0), textcoords="offset points", va="center",
                    fontsize=7 if compact else 8, color=INK)
    ax.set_xlim(0, r.max())
    ax.set_ylim(0, max(g.max(), 1.0) * 1.12)
    ax.set_xlabel(f.rule.xlabel or "r")
    ax.set_ylabel(f.rule.ylabel or "g(r)")


def draw_cumulative(ax, f: DataFile, compact: bool) -> None:
    r, n = f.x, f.ys[0]
    ax.plot(r, n, color=BLUE, lw=1.8)
    st = f.stats
    if "n_min" in st:
        ax.axvline(st["r_min"], color=INK2, lw=0.8, ls=(0, (2, 2)), alpha=0.7)
        ax.plot(st["r_min"], st["n_min"], "o", ms=6, mfc=BLUE, mec=SURFACE,
                mew=1.5, zorder=5)
        ax.annotate(f"N = {fmt_num(st['n_min'], 3)}", (st["r_min"], st["n_min"]),
                    xytext=(-7, 6), textcoords="offset points", ha="right",
                    fontsize=7 if compact else 8, color=INK)
    ax.set_xlim(0, r.max())
    ax.set_ylim(bottom=0)
    ax.set_xlabel(f.rule.xlabel or "r")
    ax.set_ylabel(f.rule.ylabel or "N(r)")


def draw_line(ax, f: DataFile, compact: bool) -> None:
    for i, (y, name) in enumerate(zip(f.ys, f.names)):
        ax.plot(f.x, y, color=SERIES[i % len(SERIES)], lw=1.6, label=name)
    ax.set_xlim(f.x.min(), f.x.max())
    ax.set_xlabel(scaled_axis(ax, f.x, f.rule.xlabel or f.xname or "column 1"))
    ax.set_ylabel(f.rule.ylabel or (f.names[0] if len(f.ys) == 1 else ""))
    if len(f.ys) > 1:
        ax.legend(loc="best")


DRAWERS = {"timeseries": draw_timeseries, "rdf": draw_rdf,
           "cumulative": draw_cumulative, "line": draw_line}


def draw(ax, f: DataFile, compact: bool = False) -> None:
    if f.status in ("error", "empty", "zero"):
        draw_placeholder(ax, f, compact)
        return
    DRAWERS.get(f.kind, draw_line)(ax, f, compact)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5 if compact else 6))
    put_title(ax, f, compact)


# --------------------------------------------------------------------------
# Report pages
# --------------------------------------------------------------------------
def grid_shape(n: int) -> tuple[int, int]:
    cols = {1: 1, 2: 2, 3: 3, 4: 2, 5: 3, 6: 3, 7: 4, 8: 4, 9: 3}.get(n, 4)
    return math.ceil(n / cols), cols


@dataclass
class PageSpec:
    section: str
    group: Group
    files: list[DataFile]
    part: int
    parts: int
    number: int = 0


def plan_pages(sections: list[Section]) -> list[PageSpec]:
    pages = []
    for sec in sections:
        for grp in sec.groups:
            chunks = [grp.files[i:i + MAX_PANELS_PER_PAGE]
                      for i in range(0, len(grp.files), MAX_PANELS_PER_PAGE)]
            for k, chunk in enumerate(chunks, 1):
                pages.append(PageSpec(sec.name, grp, chunk, k, len(chunks)))
    return pages


def page_frame(fig, title: str, subtitle_text: str, page_no: int, stamp: str):
    fig.text(0.045, 0.945, title, fontsize=18, fontweight="bold", color=INK)
    if subtitle_text:
        fig.text(0.045, 0.912, subtitle_text, fontsize=10, color=INK2)
    fig.add_artist(plt.Line2D([0.045, 0.955], [0.895, 0.895], color=GRID,
                              lw=1, transform=fig.transFigure))
    fig.text(0.045, 0.025, f"Projecte Turri · {stamp}", fontsize=7, color=MUTED)
    fig.text(0.955, 0.025, str(page_no), fontsize=8, color=MUTED, ha="right")


def render_grid_page(pdf: PdfPages, spec: PageSpec, stamp: str) -> None:
    fig = plt.figure(figsize=A4_LANDSCAPE)
    n = len(spec.files)
    rows, cols = grid_shape(n)
    gtitle = spec.group.title
    if spec.parts > 1:
        gtitle += f" ({spec.part}/{spec.parts})"
    flagged = sum(f.flagged for f in spec.files)
    sub = f"{gtitle}  ·  {n} file{'s' if n != 1 else ''}"
    if spec.group.rel != ".":
        sub += f"  ·  folder {spec.group.rel}/"
    if flagged:
        sub += f"  ·  {flagged} flagged"
    page_frame(fig, spec.section, sub, spec.number, stamp)
    # cap the row height so a single row of panels is not stretched
    bottom = max(0.085, 0.83 - 0.40 * rows)
    gs = fig.add_gridspec(rows, cols, left=0.065, right=0.965, top=0.83,
                          bottom=bottom, wspace=0.30, hspace=0.75 if rows > 1 else 0.3)
    for i, f in enumerate(spec.files):
        ax = fig.add_subplot(gs[i // cols, i % cols])
        draw(ax, f, compact=cols > 2 or rows > 2)
    pdf.savefig(fig)
    plt.close(fig)


def index_lines(sections: list[Section], pages: list[PageSpec]) -> list[tuple]:
    first_page = {}
    for p in pages:
        first_page.setdefault((p.section, p.group.rel), p.number)
    lines: list[tuple] = []
    for sec in sections:
        nflag = sum(f.flagged for f in sec.files)
        lines.append(("section", sec.name,
                      f"{len(sec.files)} files" + (f", {nflag} flagged" if nflag else "")))
        for grp in sec.groups:
            lines.append(("group", grp.title, f"p. {first_page[(sec.name, grp.rel)]}"))
            for f in grp.files:
                npts = "" if f.x is None else str(len(f.x))
                info = f.message if f.flagged else subtitle(f)
                lines.append(("file", f.path.name, f.kind, npts, info, f.status))
        lines.append(("gap",))
    return lines


LINES_PER_INDEX_PAGE = 34


def render_index(pdf: PdfPages, lines: list[tuple], n_pages: int, stamp: str,
                 root: Path) -> None:
    chunks = [lines[i:i + LINES_PER_INDEX_PAGE]
              for i in range(0, len(lines), LINES_PER_INDEX_PAGE)] or [[]]
    for k, chunk in enumerate(chunks, 1):
        fig = plt.figure(figsize=A4_LANDSCAPE)
        page_frame(fig, "Simulation plots — index",
                   f"{'/'.join(root.parts[-2:])}  ·  generated {stamp}  ·  {n_pages} pages",
                   k, stamp)
        cols = [0.06, 0.30, 0.40, 0.47]
        y = 0.855
        fig.text(cols[0], y, "FILE", fontsize=7, color=MUTED, fontweight="bold")
        fig.text(cols[1], y, "KIND", fontsize=7, color=MUTED, fontweight="bold")
        fig.text(cols[2], y, "POINTS", fontsize=7, color=MUTED, fontweight="bold")
        fig.text(cols[3], y, "KEY NUMBERS", fontsize=7, color=MUTED, fontweight="bold")
        y -= 0.03
        step = 0.0215
        for line in chunk:
            tag = line[0]
            if tag == "section":
                y -= 0.004
                fig.text(0.045, y, line[1], fontsize=12, fontweight="bold")
                fig.text(0.955, y, line[2], fontsize=8, color=INK2, ha="right")
            elif tag == "group":
                fig.text(0.05, y, line[1], fontsize=9, fontweight="semibold",
                         color=INK2)
                fig.text(0.955, y, line[2], fontsize=8, color=BLUE_DARK, ha="right")
            elif tag == "file":
                _, name, kind, npts, info, status = line
                color = ALERT if status in ("zero", "empty", "error") else INK
                fig.text(cols[0], y, name, fontsize=8.5, color=color)
                fig.text(cols[1], y, kind, fontsize=8, color=INK2)
                fig.text(cols[2], y, npts, fontsize=8, color=INK2)
                prefix = {"zero": "NO SIGNAL — ", "empty": "EMPTY — ",
                          "error": "ERROR — ", "constant": "CONSTANT — "}.get(status, "")
                fig.text(cols[3], y, prefix + info, fontsize=8,
                         color=ALERT if prefix else INK)
            y -= step
        pdf.savefig(fig)
        plt.close(fig)


def save_individual(f: DataFile, out_dir: Path) -> list[Path]:
    fig, ax = plt.subplots(figsize=(6.4, 4.0), layout="constrained")
    draw(ax, f, compact=False)
    target = out_dir / f.section / Path(f.rel).with_suffix("")
    target.parent.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext in (".png", ".pdf"):
        p = target.with_suffix(ext)
        fig.savefig(p)
        paths.append(p)
    plt.close(fig)
    return paths


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--root", type=Path, default=here,
                    help="folder that contains the section folders (default: script folder)")
    ap.add_argument("--out", type=Path, default=None,
                    help=f"output folder (default: <root>/{OUTPUT_DIRNAME})")
    ap.add_argument("--only", nargs="+", metavar="SECTION",
                    help="only these sections (folder names)")
    ap.add_argument("--no-individual", action="store_true",
                    help="skip the per-figure PNG/PDF files")
    ap.add_argument("--equil-frac", type=float, default=0.10,
                    help="fraction of a time series discarded before mean/std (default 0.10)")
    args = ap.parse_args(argv)

    root = args.root.resolve()
    out_dir = (args.out or root / OUTPUT_DIRNAME).resolve()
    apply_style()

    sections = discover(root, out_dir)
    if args.only:
        wanted = {s.lower() for s in args.only}
        sections = [s for s in sections if s.name.lower() in wanted]
    if not sections:
        print(f"No data files found under {root}")
        return 1
    analyse(sections, args.equil_frac)

    stamp = dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    pages = plan_pages(sections)
    lines = index_lines(sections, pages)  # page numbers filled below
    n_index = max(1, math.ceil(len(lines) / LINES_PER_INDEX_PAGE))
    for i, p in enumerate(pages):
        p.number = n_index + 1 + i
    lines = index_lines(sections, pages)

    out_dir.mkdir(parents=True, exist_ok=True)
    report = out_dir / REPORT_NAME
    with PdfPages(report) as pdf:
        render_index(pdf, lines, n_index + len(pages), stamp, root)
        for spec in pages:
            render_grid_page(pdf, spec, stamp)
        meta = pdf.infodict()
        meta["Title"] = "Projecte Turri - simulation plots"

    n_fig = 0
    if not args.no_individual:
        for sec in sections:
            for f in sec.files:
                save_individual(f, out_dir)
                n_fig += 1

    print(f"Report: {report}  ({n_index + len(pages)} pages)")
    for sec in sections:
        print(f"  {sec.name}: " + ", ".join(
            f"{'top level' if g.rel == '.' else g.rel} ({len(g.files)})" for g in sec.groups))
        for f in sec.files:
            if f.flagged:
                print(f"    ! {f.rel}: {f.status} - {f.message}")
    if n_fig:
        print(f"Individual figures: {n_fig} x (png + pdf) in {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
