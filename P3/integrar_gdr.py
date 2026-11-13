"""Apartado 7: U y P por integracion numerica de g(r), comparadas con el calculo directo.

U/N = 2 pi rho  Int u(r) g(r) r^2 dr
P   = rho T - (2 pi/3) rho^2  Int u'(r) g(r) r^3 dr      (unidades reducidas)
Se integra hasta rc = 2.5 (el mismo potencial truncado del MC) con la regla del
rectangulo sobre los bins de g(r), y se suman las correcciones de cola.

Uso:  python integrar_gdr.py [ncycl ...]      (por defecto 1000000 y 10000000)
"""
import sys
from pathlib import Path

import numpy as np

RUNS = Path(__file__).resolve().parent / "runs"
N, T, RC, NHIS = 500, 2.0, 2.5, 500


def integrar(path, rho):
    r, g = np.loadtxt(path, unpack=True)
    delg = (N / rho) ** (1 / 3) / (2 * NHIS)
    w = np.clip(RC - (r - delg / 2), 0, delg)  # ancho del bin dentro de [0, rc]
    u = 4 * (r**-12 - r**-6)
    du = -24 * (2 * r**-13 - r**-7)
    uN = 2 * np.pi * rho * np.sum(u * g * r**2 * w)
    p = rho * T - 2 * np.pi / 3 * rho**2 * np.sum(du * g * r**3 * w)
    rc3 = RC**-3
    uN += 8 / 3 * np.pi * rho * (rc3**3 / 3 - rc3)
    p += 16 / 3 * np.pi * rho**2 * (2 / 3 * rc3**3 - rc3)
    return uN, p


for n in sys.argv[1:] or ["1000000", "10000000"]:
    filas = np.genfromtxt(RUNS / f"summary_n{n}.csv", delimiter=",", names=True)
    out = RUNS / f"integral_n{n}.csv"
    with open(out, "w") as f:
        f.write("rho,U_direct,U_integral,P_direct,P_integral\n")
        for row in filas:
            rho = round(row["rho"], 2)
            ui, pi_ = integrar(RUNS / f"n{n}" / f"rho{rho:.2f}" / "gdr-mc.dat", rho)
            f.write(f"{rho},{row['U_per_atom']},{ui},{row['P']},{pi_}\n")
            print(f"n={n} rho={rho:.1f}  U: {row['U_per_atom']:8.4f} {ui:8.4f}   P: {row['P']:7.4f} {pi_:7.4f}")
