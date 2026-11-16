"""Genera report/tables.tex a partir de runs/*.csv y de la tabla de Johnson et al. (1993)."""
from pathlib import Path
import numpy as np
R = Path(__file__).resolve().parents[1] / "runs"
# rho: (P, sigma, U, sigma) de Johnson et al., T*=2.0
J = {0.1:("0.1777",2,-0.667,1),0.2:("0.3290",6,-1.306,1),0.4:("0.705",1,-2.538,2),0.5:("1.069",3,-3.149,2),
     0.6:("1.756",7,-3.746,2),0.7:("3.024",7,-4.300,2),0.8:("5.28",1,-4.753,3)}
def load(n): return np.genfromtxt(R/f"summary_n{n}.csv",delimiter=",",names=True)
def integ(n): return np.genfromtxt(R/f"integral_n{n}.csv",delimiter=",",names=True)
o=[]
for n,lab in ((1000000,"10^6"),(10000000,"10^7")):
    d=load(n)
    o.append(r"\begin{table}[H]\centering\small")
    o.append(rf"\caption{{Technical parameters and results for $n_{{\mathrm{{cycl}}}}={lab}$ ($N=500$, $T^*=2.0$, $r_c^*=2.5$). $\sigma_U$ and $\sigma_P$ are the standard deviations of the sampled values (2000 samples per run, one every $n_{{\mathrm{{cycl}}}}/2000$ moves). $U$ is per atom and includes the tail correction.}}")
    o.append(rf"\label{{tab:n{n}}}")
    o.append(r"\begin{tabular}{cccccccc}\hline")
    o.append(r"$\rho^*$ & $\Delta^*$ & acc.\ (\%) & $U^*/N$ & $\sigma_U$ & $P^*$ & $\sigma_P$ \\ \hline")
    for r in d:
        o.append(f"{r['rho']:.1f} & {r['delta_A']/3.405:.3f} & {100*r['acceptance']:.1f} & {r['U_per_atom']:.4f} & {r['sigma_U']:.3f} & {r['P']:.4f} & {r['sigma_P']:.3f} \\\\".replace(" & & "," & "))
    o.append(r"\hline\end{tabular}\end{table}")
# fix: header has 7 columns
txt="\n".join(o).replace("{cccccccc}","{ccccccc}")
# delta_A column holds delta in sigma units already? see below
d=load(10000000); i=integ(10000000)
o=[r"\begin{table}[H]\centering\small",
 r"\caption{Direct averages (Step~4, $n_{\mathrm{cycl}}=10^7$) against numerical integration of $g(r)$ (Step~7), both with tail corrections.}",
 r"\label{tab:integral}\begin{tabular}{ccccccc}\hline",
 r"$\rho^*$ & $U_{\mathrm{dir}}$ & $U_{\mathrm{int}}$ & $\Delta U$ & $P_{\mathrm{dir}}$ & $P_{\mathrm{int}}$ & $\Delta P$ \ \hline"]
for r in i:
    o.append(f"{r['rho']:.1f} & {r['U_direct']:.4f} & {r['U_integral']:.4f} & {r['U_integral']-r['U_direct']:+.4f} & {r['P_direct']:.4f} & {r['P_integral']:.4f} & {r['P_integral']-r['P_direct']:+.4f} \\\\")
o.append(r"\hline\end{tabular}\end{table}")
txt+="\n"+"\n".join(o)
o=[r"\begin{table}[H]\centering\small",
 r"\caption{Comparison with Johnson, Zollweg and Gubbins (1993) at $T^*=2.0$ ($10^7$ cycles; their $r_c$ is $5\sigma$ to $4.27\sigma$, ours $2.5\sigma$ plus tail). The number in parentheses is the error in the last digit(s).}",
 r"\label{tab:johnson}\begin{tabular}{cccccc}\hline",
 r"$\rho^*$ & $P^*$ (this work) & $P^*$ (Johnson) & $U^*/N$ (this work) & $U^*/N$ (Johnson) & $\Delta P/P$ (\%) \ \hline"]
for r in d:
    k=round(r['rho'],1)
    if k in J:
        p,sp,u,su=J[k]
        o.append(f"{k:.1f} & {r['P']:.4f} & {p}({sp}) & {r['U_per_atom']:.4f} & {u:.3f}({su}) & {100*(r['P']-float(p))/float(p):+.1f} \\\\")
    else:
        o.append(f"{k:.1f} & {r['P']:.4f} & --- & {r['U_per_atom']:.4f} & --- & --- \\\\")
o.append(r"\hline\end{tabular}\end{table}")
txt+="\n"+"\n".join(o)
o=[r"\begin{table}[H]\centering\small",
 r"\caption{Ideal-gas test at $T^*=2.0$ ($10^7$ cycles): compressibility factor $Z=P^*/(\rho^*T^*)$; $Z=1$ for an ideal gas.}",
 r"\label{tab:ideal}\begin{tabular}{cccc}\hline $\rho^*$ & $\rho^*T^*$ & $P^*$ & $Z$ \ \hline"]
for r in d:
    o.append(f"{r['rho']:.1f} & {2*r['rho']:.2f} & {r['P']:.4f} & {r['P']/(2*r['rho']):.3f} \\\\")
o.append(r"\hline\end{tabular}\end{table}")
txt+="\n"+"\n".join(o)
(Path(__file__).parent/"tables.tex").write_text(txt,encoding="utf-8")
print(txt[:900])
