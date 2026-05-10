import numpy as np
import subprocess
import os
import matplotlib.pyplot as plt

# ══════════════════════════════════════════════════════════════════════════════
# engine2:  eps_r(r) = 2  (r<b)  or  2 + 4*((r-b)/(R-b))^2  (r>=b)
#           rho/eps0(r) = a0*(1 - r/b)^2  (r<b)  or  0  (r>=b)
# ══════════════════════════════════════════════════════════════════════════════
params = dict(b=0.05, R=0.1, V0=0, a0=1e4, N2=5)

N_phi_Er = 2**np.arange(1, 9)   # N = 2, 4, ..., 256
N_conv   = N_phi_Er              # same range for convergence (N=256 as reference)

b = params['b']
R = params['R']
N_all = np.unique(np.concatenate([N_phi_Er, N_conv])).astype(int)

executable = './engine2'
outstr = f"engine2_b_{b:.2g}_R_{R:.2g}"
outdir = f"Scan_N1_{outstr}"
os.makedirs(outdir, exist_ok=True)
print(f"engine2  →  {outdir}")

# ── Run scan ──────────────────────────────────────────────────────────────────
for N in N_all:
    p = params.copy()
    p['N1'] = int(N)
    p['N2'] = int(N)
    out_prefix   = os.path.join(outdir, f"{outstr}_N1_{int(N)}")
    param_string = " ".join(f"{k}={v}" for k, v in p.items())
    cmd = f"{executable} trivial.in {param_string} output={out_prefix}"
    print(cmd)
    subprocess.run(cmd, shell=True)

# ── Couleurs ──────────────────────────────────────────────────────────────────
cmap       = plt.cm.viridis
n_phi      = len(N_phi_Er)
colors_phi = [cmap(i / max(n_phi - 1, 1)) for i in range(n_phi)]

# ══════════════════════════════════════════════════════════════════════════════
# Figure 1 : phi(r)
# ══════════════════════════════════════════════════════════════════════════════
fig1, ax1 = plt.subplots(figsize=(7, 5))
for N, col in zip(N_phi_Er, colors_phi):
    fname = os.path.join(outdir, f"{outstr}_N1_{int(N)}_phi.out")
    data  = np.loadtxt(fname)
    ax1.plot(data[:, 0], data[:, 1], color=col, lw=1.5, label=f"N={int(N)}")
ax1.axvline(b, color="gray", ls=":", lw=0.8, label=f"b = {b} m")
ax1.set_xlabel("r [m]")
ax1.set_ylabel(r"$\phi(r)$ [V]")
ax1.legend(fontsize=13)
ax1.grid(True, alpha=0.3)
fig1.tight_layout()
fig1.savefig(f"scan_N1_phi_{outstr}.png", dpi=150)

# ══════════════════════════════════════════════════════════════════════════════
# Figure 2 : Er(r) and Dr/eps0(r) side by side
# ══════════════════════════════════════════════════════════════════════════════
fig2, (ax2, ax2b) = plt.subplots(1, 2, figsize=(12, 5))
for N, col in zip(N_phi_Er, colors_phi):
    fname = os.path.join(outdir, f"{outstr}_N1_{int(N)}_ErDr.out")
    data  = np.loadtxt(fname)
    ax2.plot(data[:, 0], data[:, 1], color=col, lw=1.5, label=f"N={int(N)}")
    ax2b.plot(data[:, 0], data[:, 2], color=col, lw=1.5, label=f"N={int(N)}")

for ax in (ax2, ax2b):
    ax.axvline(b, color="gray", ls=":", lw=0.8, label=f"b = {b} m")
    ax.set_xlabel("r [m]")
    ax.legend(fontsize=13)
    ax.grid(True, alpha=0.3)
ax2.set_ylabel(r"$E_r$ [V/m]")
ax2b.set_ylabel(r"$D_r/\varepsilon_0$ [V/m]")
fig2.tight_layout()
fig2.savefig(f"scan_N1_ErDr_{outstr}.png", dpi=150)

plt.show()
