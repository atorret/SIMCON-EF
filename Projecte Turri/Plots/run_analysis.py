#!/usr/bin/env python3
"""
run_analysis.py - analyse every finished NAMD run and rebuild the PDF report.

For each run in RUNS (or given with --run) it writes into Plots/<Section>/:

    Energy.dat                    TS  TOTAL                (from output.log)
    Temperature.dat               TS  TEMP                 (from output.log)
    PresVSPressAvg.dat            TS  PRESSURE  PRESSAVG   (from output.log)
    Area_per_lipid.dat            TS  XY box area / lipids per leaflet (Å²)
                                  (box of every DCD frame = every 2500 steps, fallback
                                  step7_production.xst; leaflets from step5_input.pdb)
    G(r)/<pair>.dat, N(r)/<pair>.dat                       (VMD, measure gofr)
    Hydrogen_bonds.dat, Hydrogen_bonds_per_water.dat       (VMD, hbonds plugin)

then runs make_plots.py, which makes Plots/_output/Plots_report.pdf.
The VMD part is in vmd_analysis.tcl (selections, pairs and cutoffs are there).

Runs whose production has not finished yet are skipped with a message.

Usage (from any folder):
    python Plots/run_analysis.py                     # every run in RUNS, then the report
    python Plots/run_analysis.py --only Step3 Q2_N64
    python Plots/run_analysis.py --run Test=v1/charmm-gui-DMPC-NaCl-membrane/namd-ORG --only Test
    python Plots/run_analysis.py --skip-vmd          # only energy/temperature/pressure
    python Plots/run_analysis.py --step 10           # every 10th frame (quick test)
"""
from __future__ import annotations

import argparse
import math
import os
import struct
import subprocess
import sys
import tempfile
import time
from pathlib import Path

PLOTS_DIR = Path(__file__).resolve().parent
PROJECT_DIR = PLOTS_DIR.parent
VMD_SCRIPT = PLOTS_DIR / "vmd_analysis.tcl"
DEFAULT_VMD = r"C:\Program Files\University of Illinois\VMD2\vmd.exe"

# Section name in Plots/  ->  NAMD run folder (glob, relative to "Projecte Turri").
# The section name is what appears in the PDF; "Step*" sections come first.
RUNS: dict[str, str] = {
    "Step3":    "charmm-ORG/charmm-gui-*/namd",             # 128 DMPC, 303.15 K
    "Q1_T283K": "charmm-ORG/charmm-gui-*/namd-TEMP-283",    # 128 DMPC, 283 K
    "Q1_T323K": "charmm-ORG/charmm-gui-*/namd-TEMP-323",    # 128 DMPC, 323 K
    "Q2_N64":   "charmm-N-64/charmm-gui-*/namd",            #  64 DMPC, 303.15 K
    "Q2_N192":  "charmm-N-192/charmm-gui-*/namd",           # 192 DMPC, 303.15 K
}

PSF_NAME = "step5_input.psf"
PDB_NAME = "step5_input.pdb"
DCD_NAME = "step7_production.dcd"
XST_NAME = "step7_production.xst"
LIPID_RESNAME = "DMPC"
LIPID_HEAD_ATOM = "P"      # its z sign tells the leaflet (CHARMM-GUI centres the bilayer at z = 0)
LOG_NAMES = ("output.log", "step7_production.out")

# file name -> columns of the NAMD ENERGY table (x is always TS)
LOG_OUTPUTS = {
    "Energy.dat": ("TOTAL",),
    "Temperature.dat": ("TEMP",),
    "PresVSPressAvg.dat": ("PRESSURE", "PRESSAVG"),
}


def resolve_run(pattern: str) -> Path | None:
    p = Path(pattern)
    if not p.is_absolute():
        matches = sorted(PROJECT_DIR.glob(pattern))
        return matches[0] if matches else None
    return p if p.exists() else None


def find_log(run: Path) -> Path | None:
    for name in LOG_NAMES:
        if (run / name).is_file():
            return run / name
    return None


def check_run(run: Path | None) -> tuple[bool, str]:
    if run is None or not run.is_dir():
        return False, "folder not found"
    missing = [n for n in (PSF_NAME, DCD_NAME) if not (run / n).is_file()]
    log = find_log(run)
    if log is None:
        missing.append(" or ".join(LOG_NAMES))
    if missing:
        return False, "not simulated yet (missing " + ", ".join(missing) + ")"
    text = log.read_text(errors="replace")
    if "FATAL ERROR" in text:
        return False, f"{log.name} contains FATAL ERROR"
    if "WallClock:" not in text:
        return False, f"production still running or interrupted ({log.name} has no WallClock line)"
    return True, ""


# --------------------------------------------------------------------------
# NAMD log -> Energy / Temperature / Pressure
# --------------------------------------------------------------------------
def parse_namd_log(log: Path) -> tuple[list[str], list[list[str]]]:
    titles: list[str] = []
    rows: list[list[str]] = []
    for line in log.read_text(errors="replace").splitlines():
        if line.startswith("ETITLE:") and not titles:
            titles = line.split()[1:]
        elif line.startswith("ENERGY:"):
            rows.append(line.split()[1:])
    if not titles or not rows:
        raise ValueError(f"no ETITLE/ENERGY lines in {log}")
    return titles, rows


def write_log_data(log: Path, out: Path) -> None:
    titles, rows = parse_namd_log(log)
    for fname, cols in LOG_OUTPUTS.items():
        idx = [titles.index("TS")] + [titles.index(c) for c in cols]
        with open(out / fname, "w") as fh:
            for r in rows:
                fh.write(" ".join(r[i] for i in idx) + "\n")
    print(f"    log: {len(rows)} energy rows from {log.name} -> "
          + ", ".join(LOG_OUTPUTS))


# --------------------------------------------------------------------------
# Area per lipid: XY box area / lipids per leaflet
# --------------------------------------------------------------------------
def lipids_per_leaflet(pdb: Path) -> tuple[int, int]:
    upper = lower = 0
    for line in pdb.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:21].strip() != LIPID_RESNAME or line[12:16].strip() != LIPID_HEAD_ATOM:
            continue
        if float(line[46:54]) > 0:
            upper += 1
        else:
            lower += 1
    return upper, lower


def xy_area_from_xst(xst: Path) -> tuple[list[tuple[int, float]], int]:
    """(step, XY area) for each readable line of the .xst, and the number of bad lines."""
    rows, bad = [], 0
    for line in xst.read_text(errors="replace").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        try:
            v = line.split()
            step, ax, ay, bx, by = int(float(v[0])), *map(float, (v[1], v[2], v[4], v[5]))
        except (ValueError, IndexError):
            bad += 1                      # e.g. NUL bytes left by a crash
            continue
        rows.append((step, abs(ax * by - ay * bx)))   # |a x b|, any XY cell
    return rows, bad


def xy_area_from_dcd(dcd: Path) -> list[tuple[int, float]]:
    """(step, XY area) for each frame, from the unit cell stored in a CHARMM/NAMD DCD.

    Stops at the first frame with an empty (zero) box: after a crash the rest of
    the file can be NUL bytes, and those frames are not real data.
    """
    rows = []
    with open(dcd, "rb") as fh:
        def record() -> bytes:
            (n,) = struct.unpack("<i", fh.read(4))
            data = fh.read(n)
            fh.read(4)
            return data
        head = record()
        if head[:4] != b"CORD":
            raise ValueError(f"{dcd.name} is not a DCD file")
        icntrl = struct.unpack("<20i", head[4:84])
        nset, istart, nsavc, has_cell = icntrl[0], icntrl[1], icntrl[2], icntrl[10]
        if not has_cell:
            raise ValueError(f"{dcd.name} has no unit cell data")
        record()                                         # title
        (natom,) = struct.unpack("<i", record())
        coord_block = 3 * (4 + 4 * natom + 4)
        for i in range(nset):
            raw = fh.read(4)
            if len(raw) < 4:
                break
            cell = struct.unpack("<6d", fh.read(48))
            fh.read(4)
            a, gamma, b = cell[0], cell[1], cell[2]
            if a <= 0 or b <= 0:
                break
            # NAMD stores the angle as a cosine; older files store degrees
            sin_g = math.sqrt(max(0.0, 1 - gamma**2)) if abs(gamma) <= 1 \
                else math.sin(math.radians(gamma))
            rows.append((istart + i * nsavc, a * b * sin_g))
            fh.seek(coord_block, 1)
    return rows


def write_area_per_lipid(run: Path, out: Path) -> None:
    target = out / "Area_per_lipid.dat"
    target.unlink(missing_ok=True)
    pdb = run / PDB_NAME
    if not pdb.is_file():
        print(f"    area: skipped ({PDB_NAME} missing)")
        return
    upper, lower = lipids_per_leaflet(pdb)
    if upper == 0 or lower == 0:
        print(f"    area: skipped (found {upper} upper / {lower} lower {LIPID_RESNAME})")
        return
    n_leaflet = (upper + lower) / 2   # equal leaflets in CHARMM-GUI builds

    # The DCD has a box for every saved frame (dcdfreq = every 2500 steps), the
    # same frames as g(r) and H-bonds; the .xst (every 5000 steps) is the fallback.
    rows, source = [], DCD_NAME
    try:
        rows = xy_area_from_dcd(run / DCD_NAME)
    except (OSError, ValueError, struct.error) as exc:
        print(f"    area: could not read the box from {DCD_NAME} ({str(exc)[:100]})")
    if not rows and (run / XST_NAME).is_file():
        rows, bad = xy_area_from_xst(run / XST_NAME)
        source = XST_NAME
        if bad:
            print(f"    area: {bad} unreadable line(s) in {XST_NAME}")
    if not rows:
        print("    area: skipped (no readable box in the DCD or the XST)")
        return

    target.write_text("".join(f"{s} {area / n_leaflet:.4f}\n" for s, area in rows))
    note = "" if upper == lower else f" (unequal leaflets: {upper}/{lower}, using the mean)"
    print(f"    area: {len(rows)} boxes from {source}, {upper}+{lower} lipids -> "
          f"Area_per_lipid.dat{note}")


# --------------------------------------------------------------------------
# VMD
# --------------------------------------------------------------------------
def tcl_path(p: Path) -> str:
    return "{" + p.as_posix() + "}"


def dcd_frame_counts(dcd: Path) -> tuple[int, int, int]:
    """(frames in header, valid frames, last valid step) of a DCD."""
    with open(dcd, "rb") as fh:
        fh.read(4)
        head = fh.read(84)
    nset = struct.unpack("<i", head[4:8])[0]
    valid = xy_area_from_dcd(dcd)
    return nset, len(valid), (valid[-1][0] if valid else -1)


def run_vmd(vmd: str, run: Path, out: Path, step: int, last: int = -1) -> bool:
    # Remove old structure files so a skipped pair cannot leave stale data behind.
    for sub in ("G(r)", "N(r)"):
        for f in (out / sub).glob("*.dat"):
            f.unlink()

    header = (
        f"set PSF {tcl_path(run / PSF_NAME)}\n"
        f"set DCD {tcl_path(run / DCD_NAME)}\n"
        f"set OUT {tcl_path(out)}\n"
        f"set STEP {step}\n"
        f"set LAST {last}\n"
        f"if {{[catch {{source {tcl_path(VMD_SCRIPT)}}} err]}} {{\n"
        f"    puts \"ANALYSIS) ERROR $err\"\n"
        f"}}\n"
        f"quit\n"
    )
    fd, header_path = tempfile.mkstemp(suffix=".tcl", prefix="vmd_run_")
    with os.fdopen(fd, "w") as fh:
        fh.write(header)

    ok = True
    t0 = time.time()
    try:
        proc = subprocess.Popen(
            [vmd, "-dispdev", "text", "-e", header_path],
            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True, errors="replace",
            cwd=str(out))
        for line in proc.stdout:
            line = line.rstrip()
            if line.startswith("ANALYSIS)"):
                print("    vmd: " + line[len("ANALYSIS)"):].strip())
                if "ERROR" in line:
                    ok = False
            elif line.startswith("ERROR)") and "Duplicate resname" not in line:
                print("    vmd: " + line)
                ok = False
        proc.wait()
    finally:
        os.unlink(header_path)

    if not (out / "Hydrogen_bonds.dat").is_file():
        print("    vmd: ERROR no Hydrogen_bonds.dat written")
        ok = False
    print(f"    vmd: finished in {time.time() - t0:.0f} s")
    return ok


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", action="append", default=[], metavar="NAME=PATH",
                    help="extra run: section name and NAMD folder (repeatable)")
    ap.add_argument("--only", nargs="+", metavar="NAME",
                    help="only these sections")
    ap.add_argument("--skip-vmd", action="store_true",
                    help="only the NAMD log data (no g(r), N(r), H-bonds)")
    ap.add_argument("--no-report", action="store_true",
                    help="do not run make_plots.py at the end")
    ap.add_argument("--step", type=int, default=1,
                    help="use every N-th trajectory frame (default 1 = all)")
    ap.add_argument("--vmd", default=os.environ.get("VMD", DEFAULT_VMD),
                    help="path to vmd.exe (default: $VMD or the VMD2 install)")
    args = ap.parse_args(argv)

    runs = dict(RUNS)
    for item in args.run:
        name, sep, path = item.partition("=")
        if not sep or not name or not path:
            ap.error(f"--run expects NAME=PATH, got {item!r}")
        runs[name] = path
    if args.only:
        unknown = [n for n in args.only if n not in runs]
        if unknown:
            ap.error("unknown section(s): " + ", ".join(unknown)
                     + ". Known: " + ", ".join(runs))
        runs = {n: runs[n] for n in args.only}

    if not args.skip_vmd and not Path(args.vmd).is_file():
        print(f"VMD not found at {args.vmd}. Use --vmd or set VMD, or --skip-vmd.")
        return 1

    done, failed, skipped = [], [], []
    for name, pattern in runs.items():
        run = resolve_run(pattern)
        print(f"[{name}] {run.relative_to(PROJECT_DIR) if run and run.is_relative_to(PROJECT_DIR) else (run or pattern)}")
        ok, why = check_run(run)
        if not ok:
            print(f"    skipped: {why}")
            skipped.append(name)
            continue
        out = PLOTS_DIR / name
        out.mkdir(exist_ok=True)
        try:
            write_log_data(find_log(run), out)
            write_area_per_lipid(run, out)
        except (ValueError, OSError, IndexError, struct.error) as exc:
            print(f"    ERROR reading log/box: {str(exc)[:200]}")
            failed.append(name)
            continue
        if not args.skip_vmd:
            last = -1
            try:
                nset, nvalid, last_step = dcd_frame_counts(run / DCD_NAME)
            except (OSError, ValueError, struct.error) as exc:
                print(f"    ERROR reading {DCD_NAME}: {str(exc)[:200]}")
                failed.append(name)
                continue
            if nvalid == 0:
                print(f"    ERROR {DCD_NAME} has no valid frames")
                failed.append(name)
                continue
            if nvalid < nset:
                last = nvalid - 1
                print(f"    WARNING {DCD_NAME}: only {nvalid} of {nset} frames are valid "
                      f"(empty after step {last_step}, probably a crash). "
                      f"Analysing frames 0-{last} only.")
            if not run_vmd(args.vmd, run, out, args.step, last):
                failed.append(name)
                continue
        done.append(name)

    print()
    print(f"Analysed: {', '.join(done) or '-'}")
    if skipped:
        print(f"Skipped:  {', '.join(skipped)}")
    if failed:
        print(f"FAILED:   {', '.join(failed)}")

    if not args.no_report and done:
        print("\nBuilding report ...")
        sys.path.insert(0, str(PLOTS_DIR))
        import make_plots
        make_plots.main([])
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
