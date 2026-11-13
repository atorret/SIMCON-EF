"""Serie de simulaciones MC (NVT) a T* fija para varias densidades.

Para cada rho*: genera una red FCC de 500 atomos (mc-conf.data), ajusta delta
hasta ~50 % de aceptacion con tandas cortas, lanza la produccion con mcnvt.exe
y guarda todo en runs/n<ncycl>/rho<rho>/. Resume en runs/summary_n<ncycl>.csv
(ncycl, delta, aceptacion, cutoff, U, sigma_U, P, sigma_P) leyendo la salida
por pantalla del programa, sin tocar el Fortran.

Uso:  python run_series.py [--ncycl 1000000] [--rhos 0.1 0.2 ...] [--tstar 2.0]
"""
import argparse
import itertools
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
GFORTRAN_BIN = r"C:\Program Files (x86)\Simply Fortran 3\mingw-w64\bin"
SIGMA, EPSIL = 3.405, 119.8
NATOMS = 500
KEYS = {
    "T*": r"T\*\s*=\s*(\S+)",
    "rho*": r"rho\*\s*=\s*(\S+)",
    "rc*": r"rc\*\s*=\s*(\S+)",
    "delta*": r"delta\*\s*=\s*(\S+)",
    "acc": r"acceptance\s*=\s*(\S+)",
    "U": r"<U\*>/N\s*=\s*(\S+)\s+sigma_U\s*=\s*(\S+)",
    "P": r"<P\*>\s*=\s*(\S+)\s+sigma_P\s*=\s*(\S+)",
}

try:
    sys.path.insert(0, str(HERE.parents[2] / "scripts"))
    from progreso import marcar, terminar
except Exception:  # la barra es opcional
    def marcar(*a, **k): pass
    def terminar(*a, **k): pass


def compile_exe() -> Path:
    RUNS.mkdir(exist_ok=True)
    exe = RUNS / "mcnvt.exe"
    env = dict(os.environ, PATH=os.environ["PATH"] + os.pathsep + GFORTRAN_BIN)
    # el include '../chdir_to_code.inc' se resuelve desde P3/
    subprocess.run(["gfortran", "-static", "-O2", str(HERE / "mcnvt.f90"), "-o", str(exe)],
                   check=True, cwd=HERE, env=env)
    return exe


def write_fcc(path: Path, rho: float, n: int = 5) -> None:
    box = (4 * n**3 / rho) ** (1 / 3)  # en sigma
    a = box / n
    basis = [(0, 0, 0), (.5, .5, 0), (.5, 0, .5), (0, .5, .5)]
    with open(path, "w") as f:
        for i, j, k in itertools.product(range(n), repeat=3):
            for x, y, z in basis:
                f.write(" ".join(f"{(c + b + .25) * a * SIGMA:.10f}"
                                 for c, b in zip((i, j, k), (x, y, z), strict=True)) + "\n")
        f.write(f"{box * SIGMA:.10f}\n")


def run_mc(exe: Path, d: Path, ncycl: int, nsamp: int, tstar: float, delta: float) -> dict:
    d.mkdir(parents=True, exist_ok=True)
    shutil.copy(exe, d / "mcnvt.exe")
    (d / "mc-lj.data").write_text(
        f"{ncycl},{nsamp}\n{NATOMS}\n{tstar * EPSIL}d0\n{SIGMA}d0,{EPSIL}d0\n{delta}d0\n")
    out = subprocess.run([str(d / "mcnvt.exe")], capture_output=True, text=True, cwd=d, check=True).stdout
    (d / "stdout.txt").write_text(out)
    res = {}
    for k, pat in KEYS.items():
        m = re.search(pat, out)
        vals = [float(v.replace("D", "E")) for v in m.groups()]
        res[k] = vals if len(vals) > 1 else vals[0]
    return res


def one_density(args) -> list:
    exe, rho, ncycl, tstar = args
    base = RUNS / f"n{ncycl}" / f"rho{rho:.2f}"
    base.mkdir(parents=True, exist_ok=True)
    write_fcc(base / "mc-conf.data", rho)
    # ajuste de delta (aceptacion ~ 50 %) con tandas cortas
    delta = 0.75
    for _ in range(6):
        t = base / "tune"
        t.mkdir(exist_ok=True)
        shutil.copy(base / "mc-conf.data", t / "mc-conf.data")
        acc = run_mc(exe, t, 30000, 10000, tstar, delta)["acc"]
        if abs(acc - 0.5) < 0.03:
            break
        delta = min(max(delta * (max(acc, 0.05) / 0.5) ** 1.5, 0.05), 3.0)
    shutil.rmtree(base / "tune")
    res = run_mc(exe, base, ncycl, max(ncycl // 2000, 100), tstar, delta)
    return [ncycl, delta, res["acc"], res["rc*"], res["rho*"], res["T*"],
            *res["U"], *res["P"]]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncycl", type=int, default=1_000_000)
    ap.add_argument("--rhos", type=float, nargs="+", default=[.1, .2, .3, .4, .5, .6, .7, .8])
    ap.add_argument("--tstar", type=float, default=2.0)
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    exe = compile_exe()
    tarea = f"MC n={a.ncycl:g}"
    marcar(tarea, 0, len(a.rhos), "ajustando delta y simulando")
    rows = []
    with ProcessPoolExecutor(a.workers) as ex:
        for i, r in enumerate(ex.map(one_density, [(exe, rho, a.ncycl, a.tstar) for rho in a.rhos]), 1):
            rows.append(r)
            marcar(tarea, i, len(a.rhos), f"rho*={r[4]:.2f} listo")
    terminar(tarea)
    out = RUNS / f"summary_n{a.ncycl}.csv"
    with open(out, "w") as f:
        f.write("ncycl,delta_A,acceptance,rc,rho,T,U_per_atom,sigma_U,P,sigma_P\n")
        for r in rows:
            f.write(",".join(str(x) for x in r) + "\n")
    print("Resumen en", out)


if __name__ == "__main__":
    main()
