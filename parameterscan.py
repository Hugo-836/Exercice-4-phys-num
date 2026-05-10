import numpy as np
import subprocess
import os
import matplotlib.pyplot as plt

# ══════════════════════════════════════════════════════════════════════════════
# MODE : choisir 'trivial' ou 'nontrivial'
# ══════════════════════════════════════════════════════════════════════════════
MODE = 'trivial'   # <-- changer ici

# ── Paramètres selon le mode ──────────────────────────────────────────────────
if MODE == 'trivial':
    params = dict(b=0.05, R=0.1, V0=0, a0=1, trivial='true', N2=5)
    N_phi_Er    = np.array([5])             # N used for phi and Er plots
    N_conv      = 2**np.arange(1, 9)       # N = 2, 4, ..., 256 for convergence
    phi0_exact  = params['R']**2 / 4       # analytical phi(0)
    input_file  = 'trivial.in'
else:
    params = dict(b=0.02, R=0.1, V0=0, a0=1e4, trivial='false', N2=5)
    N_phi_Er    = 2**np.arange(1, 9)       # N = 2, 4, ..., 256
    N_conv      = N_phi_Er
    phi0_exact  = None                     # no analytical solution, use N=256
    input_file  = 'nontrivial.in'

b = params['b']
R = params['R']
N_all = np.unique(np.concatenate([N_phi_Er, N_conv])).astype(int)

executable = './engine'
outstr = (f"electrostatics_b_{b:.2g}_R_{R:.2g}_trivial_{params['trivial']}")
outdir = f"Scan_N1_{outstr}"
os.makedirs(outdir, exist_ok=True)
print(f"Mode: {MODE}  →  {outdir}")

# ── Run scan ──────────────────────────────────────────────────────────────────
for N in N_all:
    p = params.copy()
    p['N1'] = int(N)
    if p['N2'] != 0:
        p['N2'] = int(N)
    out_prefix   = os.path.join(outdir, f"{outstr}_N1_{int(N)}")
    param_string = " ".join(f"{k}={v}" for k, v in p.items())
    cmd = f"{executable} {input_file} {param_string} output={out_prefix}"
    print(cmd)
    subprocess.run(cmd, shell=True)

# ── Couleurs ──────────────────────────────────────────────────────────────────
cmap        = plt.cm.viridis
n_phi       = len(N_phi_Er)
colors_phi  = [cmap(i / max(n_phi - 1, 1)) for i in range(n_phi)]

# ══════════════════════════════════════════════════════════════════════════════
# Figure 1 : phi(r)
# ══════════════════════════════════════════════════════════════════════════════
fig1, ax1 = plt.subplots(figsize=(7, 5))
for N, col in zip(N_phi_Er, colors_phi):
    fname = os.path.join(outdir, f"{outstr}_N1_{int(N)}_phi.out")
    data  = np.loadtxt(fname)
    ax1.plot(data[:, 0], data[:, 1], color=col, lw=1.5, label=f"N={int(N)} (numerical)")

if MODE == 'trivial':
    r_ana = np.linspace(0, R, 500)
    ax1.plot(r_ana, (R**2 - r_ana**2) / 4, "k--", lw=1.5, label="analytical")

ax1.axvline(b, color="gray", ls=":", lw=0.8, label=f"b = {b} m")
ax1.set_xlabel("r [m]")
ax1.set_ylabel(r"$\phi(r)$ [V]")
ax1.legend(fontsize=13)
ax1.grid(True, alpha=0.3)
fig1.tight_layout()
fig1.savefig(f"scan_N1_phi_{MODE}.png", dpi=150)

# ══════════════════════════════════════════════════════════════════════════════
# Figure 2 : Er(r)  [+ Dr/eps0 pour non-trivial]
# ══════════════════════════════════════════════════════════════════════════════
if MODE == 'nontrivial':
    fig2, (ax2, ax2b) = plt.subplots(1, 2, figsize=(12, 5))
else:
    fig2, ax2 = plt.subplots(figsize=(7, 5))

for N, col in zip(N_phi_Er, colors_phi):
    fname = os.path.join(outdir, f"{outstr}_N1_{int(N)}_ErDr.out")
    data  = np.loadtxt(fname)
    ax2.plot(data[:, 0], data[:, 1], color=col, lw=1.5, label=f"N={int(N)} (numerical)")
    if MODE == 'nontrivial':
        ax2b.plot(data[:, 0], data[:, 2], color=col, lw=1.5, label=f"N={int(N)} (numerical)")

if MODE == 'trivial':
    r_ana = np.linspace(0, R, 500)
    ax2.plot(r_ana, r_ana / 2, "k--", lw=1.5, label="analytical")

ax2.axvline(b, color="gray", ls=":", lw=0.8, label=f"b = {b} m")
ax2.set_xlabel("r [m]")
ax2.set_ylabel(r"$E_r$ [V/m]")
ax2.legend(fontsize=13)
ax2.grid(True, alpha=0.3)

if MODE == 'nontrivial':
    ax2b.axvline(b, color="gray", ls=":", lw=0.8, label=f"b = {b} m")
    ax2b.set_xlabel("r [m]")
    ax2b.set_ylabel(r"$D_r/\varepsilon_0$ [V/m]")
    ax2b.legend(fontsize=13)
    ax2b.grid(True, alpha=0.3)

fig2.tight_layout()
if MODE == 'trivial':
    fig2.savefig(f"scan_N1_Er_{MODE}.png", dpi=150)
else:
    fig2.savefig(f"scan_N1_ErDr_{MODE}.png", dpi=150)

# ══════════════════════════════════════════════════════════════════════════════
# Figure 3 : divD/eps0 et rho/eps0  (non-trivial uniquement)
# ══════════════════════════════════════════════════════════════════════════════
if MODE == 'nontrivial':
    n_all   = len(N_phi_Er)
    colors_all = [cmap(i / max(n_all - 1, 1)) for i in range(n_all)]
    fig3, ax3 = plt.subplots(figsize=(7, 5))
    for N, col in zip(N_phi_Er, colors_all):
        fname = os.path.join(outdir, f"{outstr}_N1_{int(N)}_divD_rho.out")
        data  = np.loadtxt(fname)
        ax3.plot(data[:, 0], data[:, 1], color=col, lw=1.5, label=f"N={int(N)}")
    ax3.plot(data[:, 0], data[:, 2], "k--", lw=1.5, label=r"$\rho/\varepsilon_0$")
    ax3.axvline(b, color="gray", ls=":", lw=0.8, label=f"b = {b} m")
    ax3.set_xlabel("r [m]")
    ax3.set_ylabel(r"[V/m²]")
    ax3.legend(fontsize=13)
    ax3.grid(True, alpha=0.3)
    fig3.tight_layout()
    fig3.savefig(f"scan_N1_divD_rho_{MODE}.png", dpi=150)

# ══════════════════════════════════════════════════════════════════════════════
# Figure 4 : convergence phi(0)
# ══════════════════════════════════════════════════════════════════════════════
phi0_vals = []
for N in N_conv:
    fname = os.path.join(outdir, f"{outstr}_N1_{int(N)}_phi.out")
    phi0_vals.append(np.loadtxt(fname)[0, 1])

N_arr = np.array(N_conv, dtype=float)

if phi0_exact is not None:
    ref       = phi0_exact
    ref_label = r"$|\phi(0) - \phi_{\rm exact}|$"
    err       = np.abs(np.array(phi0_vals) - ref)
    N_plot    = N_arr
else:
    ref       = phi0_vals[-1]   # N=256 as reference
    ref_label = r"$|\phi(0) - \phi_{N=256}|$"
    err       = np.abs(np.array(phi0_vals[:-1]) - ref)
    N_plot    = N_arr[:-1]

fig4, ax4 = plt.subplots(figsize=(7, 5))
ax4.loglog(N_plot, err + 1e-16, "o-", color="steelblue", ms=8, label=ref_label)
x_g = np.array([N_plot[0], N_plot[-1]])
ax4.loglog(x_g, (err[0] + 1e-16) * (x_g / N_plot[0])**(-2),
           "k--", lw=0.8, label=r"$\propto N^{-2}$")
ax4.set_xlabel("N  (= N1 = N2)")
ax4.set_ylabel(r"Error on $\phi(0)$")
ax4.legend(fontsize=13)
ax4.grid(True, which="both", alpha=0.3)
fig4.tight_layout()
fig4.savefig(f"scan_N1_convergence_{MODE}.png", dpi=150)

plt.show()
