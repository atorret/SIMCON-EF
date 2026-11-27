# Independent check of the Fortran fcc energy: direct lattice sum (no PBC code involved)
import numpy as np
rho, rc = 1.1, 3.8444                      # rc printed by the Fortran run (= L/2)
a = (4/rho)**(1/3); rnn = a/np.sqrt(2)
n = 8                                      # enough cells to cover far beyond rc
g = np.arange(-n, n+1)
cells = np.array(np.meshgrid(g, g, g, indexing='ij')).reshape(3, -1).T
basis = np.array([[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5]])
R = (cells[:, None, :] + basis[None]).reshape(-1, 3) * a
d = np.linalg.norm(R, axis=1); d = d[d > 1e-9]   # distances from atom at origin
# 1) structure: nearest neighbours
print(f"a = {a:.5f}  r_nn = {rnn:.5f}  min dist = {d.min():.5f}  #NN = {(abs(d-rnn)<1e-6).sum()}")
u = lambda r: 4*(r**-12 - r**-6)
# 2) truncated lattice sum, half because each pair is shared by 2 atoms
U_trunc = 0.5*u(d[d < rc]).sum()
# 3) full lattice sum (converged) and Lennard-Jones-Ingham constants A12, A6
dd = d/rnn
A12, A6 = (dd**-12).sum(), (dd**-6)[d < n*a*0.9].sum()
U_inf = 2*(A12*rnn**-12 - A6*rnn**-6)
utail = 8/3*np.pi*rho*(1/3/rc**9 - 1/rc**3)
print(f"A12 = {A12:.4f} (lit. 12.1319)   A6 = {A6:.4f} (lit. 14.4539)")
print(f"U/N truncated lattice sum  = {U_trunc:.6f}   (Fortran: -8.436868)")
print(f"U/N truncated + tail       = {U_trunc+utail:.6f}   (Fortran: -8.599041)")
print(f"U/N infinite lattice       = {U_inf:.6f}   -> tail-correction error = {U_trunc+utail-U_inf:+.4f}")
