# Independent check of the Fortran S(k): full k sphere (not half space), plain numpy
import sys, numpy as np
run = sys.argv[1]
sigma = 3.405
raw = np.loadtxt(f"{run}/newconf.data", max_rows=500)          # final configuration, Angstrom
L = float(open(f"{run}/newconf.data").read().split()[-1]) / sigma
r = raw / sigma; N = len(r); ncell = 5
dk = 2*np.pi/L; nmax = int(20.0/dk)
g = np.arange(-nmax, nmax+1)
n = np.array(np.meshgrid(g, g, g, indexing='ij')).reshape(3, -1).T
n2 = (n**2).sum(1); keep = (n2 > 0) & (n2 <= nmax**2); n, n2 = n[keep], n2[keep]
rho = np.exp(1j*dk*(n @ r.T)).sum(1)                              # rho(k) for every k
S = np.abs(rho)**2 / N
# shell average over n^2
ref = np.loadtxt(f"{run}/sk-mc_T1.00.dat")
shells = np.unique(n2)
Spy = np.array([S[n2 == s].mean() for s in shells]); nv = np.array([(n2 == s).sum() for s in shells])
assert np.allclose(ref[:, 0], dk*np.sqrt(shells)) and np.all(ref[:, 2]*2 == nv)
print("S(k) shells: max |diff| =", np.abs(ref[:, 1]-Spy).max(), "  max S =", Spy.max())
# Bragg families
hkl_ok = np.all(n % ncell == 0, 1)
h = n // ncell; par = np.all(h % 2 == (h[:, :1] % 2), 1); bragg = hkl_ok & par
m2 = (h**2).sum(1)
br = np.loadtxt(f"{run}/bragg-mc_T1.00.dat")
for row in br:
    sel = bragg & (m2 == int(row[0]))
    print(f"  hkl^2={int(row[0]):3d}  S(G)/N  fortran {row[2]:.10f}  python {S[sel].mean()/N:.10f}")
sel = bragg & (m2 == 3)
print("Phi python =", np.abs(rho[sel]).mean()/N, "  (fortran stdout 0.94886)")
