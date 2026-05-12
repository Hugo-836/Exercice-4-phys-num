import numpy as np
import subprocess
import os
import matplotlib.pyplot as plt

# ══════════════════════════════════════════════════════════════════════════════
# MODE : choisir 'trivial' ou 'nontrivial'
# ══════════════════════════════════════════════════════════════════════════════
MODE = 'nontrivial'   # <-- changer ici

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
    ax1.plot(data[:, 0], data[:, 1], color=col, lw=1.5, label=f"N={int(N)}")

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
    fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(12, 5))
else:
    fig2, ax2a = plt.subplots(figsize=(7, 5))

for N, col in zip(N_phi_Er, colors_phi):
    fname = os.path.join(outdir, f"{outstr}_N1_{int(N)}_ErDr.out")
    data  = np.loadtxt(fname)
    ax2a.plot(data[:, 0], data[:, 1], color=col, lw=1.5, label=f"N={int(N)}")
    if MODE == 'nontrivial':
        ax2b.plot(data[:, 0], data[:, 2], color=col, lw=1.5, label=f"N={int(N)}")

if MODE == 'trivial':
    r_ana = np.linspace(0, R, 500)
    ax2a.plot(r_ana, r_ana / 2, "k--", lw=1.5, label="analytical")

ax2a.axvline(b, color="gray", ls=":", lw=0.8, label=f"b = {b} m")
ax2a.set_xlabel("r [m]")
ax2a.set_ylabel(r"$E_r$ [V/m]")
ax2a.legend(fontsize=13)
ax2a.grid(True, alpha=0.3)

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
    fig3, ax3a = plt.subplots(figsize=(7, 5))
    for N, col in zip(N_phi_Er, colors_all):
        fname = os.path.join(outdir, f"{outstr}_N1_{int(N)}_divD_rho.out")
        data  = np.loadtxt(fname)
        ax3a.plot(data[:, 0], data[:, 1], color=col, lw=1.5, label=f"N={int(N)}")
    ax3a.plot(data[:, 0], data[:, 2], "k--", lw=1.5, label=r"$\rho/\varepsilon_0$")
    ax3a.axvline(b, color="gray", ls=":", lw=0.8, label=f"b = {b} m")
    ax3a.set_xlabel("r [m]")
    ax3a.set_ylabel(r"[V/m²]")
    ax3a.legend(fontsize=13)
    ax3a.grid(True, alpha=0.3)
    fig3.tight_layout()
    fig3.savefig(f"scan_N1_divD_rho_{MODE}.png", dpi=150)

if MODE == 'nontrivial':
    n_all   = len(N_phi_Er)
    colors_all = [cmap(i / max(n_all - 1, 1)) for i in range(n_all)]
    fig3, ax3b = plt.subplots(figsize=(7, 5))
    for N, col in zip(N_phi_Er, colors_all):
        fname = os.path.join(outdir, f"{outstr}_N1_{int(N)}_divD_rho.out")
        data  = np.loadtxt(fname)
        ax3b.plot(data[:, 0], np.abs(data[:, 1]-data[:, 2]), color=col, lw=1.5, label=f"N={int(N)}")
    ax3b.axvline(b, color="gray", ls=":", lw=0.8, label=f"b = {b} m")
    ax3b.set_xlabel("r [m]")
    ax3b.set_ylabel(r"[V/m²]")
    ax3b.legend(fontsize=13)
    ax3b.grid(True, alpha=0.3)
    ax3b.loglog()
    fig3.tight_layout()
    fig3.savefig(f"scan_N1_error_divD_rho_{MODE}.png", dpi=150)

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
    ref_label = r"$|\phi(0) - \phi_{\rm ana}|$"
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

# ══════════════════════════════════════════════════════════════════════════════
# Figure 5 : diagramme en bâton — Q_pol, Q_lib, Q_tot  (non-trivial)
#   Chacune calculée indépendamment :
#   Q_lib : ∫ ρ_lib/ε₀ · r dr   (intégrale directe de la densité)
#   Q_tot : r · E_r(R)           (loi de Gauss dans le vide à r=R)
#   Q_pol : −r · P_r(R)         (flux de polarisation à r=R, P_r = D_r − E_r)
# ══════════════════════════════════════════════════════════════════════════════
if MODE == 'nontrivial':
    N_hist  = int(N_conv[-1])
    fname   = os.path.join(outdir, f"{outstr}_N1_{N_hist}_ErDr.out")
    data    = np.loadtxt(fname)
    r_mid_h = data[:, 0]
    Er_h    = data[:, 1]
    Dr_h    = data[:, 2]   # D_r = ε_r · E_r

    # Largeurs des intervalles
    h1_val = 2.0 * r_mid_h[0]
    h2_val = 2.0 * (R - r_mid_h[-1])
    N1_h   = int(np.round(b / h1_val))
    h_arr  = np.where(np.arange(len(r_mid_h)) < N1_h, h1_val, h2_val)

    eps0 = 8.854187817e-12   # F/m

    # 1. Q_lib/L_z : intégrale de ρ_lib · r sur chaque cellule × 2π
    rho_h = np.where(r_mid_h <= b,
                     params['a0'] * np.sin(np.pi * r_mid_h / b),
                     0.0)
    Q_lib = 2 * np.pi * eps0 * np.sum(rho_h * r_mid_h * h_arr)

    # 2. Q_tot/L_z : loi de Gauss avec E_r au bord r ≈ R
    Q_tot = 2 * np.pi * eps0 * r_mid_h[-1] * Er_h[-1]

    # 3. Q_pol/L_z : flux de polarisation à r = b  (premier midpoint extérieur, P_r = D_r − E_r)
    Q_pol = -2 * np.pi * eps0 * r_mid_h[N1_h] * (Dr_h[N1_h] - Er_h[N1_h])

    # ── Valeurs analytiques ───────────────────────────────────────────────────
    a0  = params['a0']
    # ε_r(b⁺) = 3,  ε_r(R) = 9  (formule engine : 3 + 6*(r-b)/(R-b))
    eps_b = 3.0 + 6.0 * 0.0          # = 3
    eps_R = 3.0 + 6.0 * 1.0          # = 9
    # Q_lib_ana/L_z = 2πε₀ · a0·b²/π = 2ε₀·a0·b²
    Q_lib_ana = 2 * eps0 * a0 * b**2
    # D_r(b) = a0·b/π  →  E_r(b⁺) = a0·b/(ε_b·π)  →  P_r(b⁺) = (ε_b−1)·E_r(b⁺)
    # Q_pol_ana/L_z = −2πε₀·b·P_r(b⁺) = −4ε₀·a0·b²/3
    Q_pol_ana = -2 * eps0 * a0 * b**2 * (eps_b - 1) / eps_b
    # Q_tot_ana/L_z = 2πε₀·R·E_r(R) = 2ε₀·a0·b²/ε_R
    Q_tot_ana = 2 * eps0 * a0 * b**2 / eps_R

    fig5, ax5 = plt.subplots(figsize=(6, 5))
    labels   = [r"$Q_{\rm pol}$", r"$Q_{\rm lib}$", r"$Q_{\rm tot}$"]
    num_vals = [Q_pol, Q_lib, Q_tot]
    ana_vals = [Q_pol_ana, Q_lib_ana, Q_tot_ana]
    colors   = ["orange", "tomato", "steelblue"]

    ax5.bar(labels, num_vals, color=colors, width=0.4, edgecolor="k", linewidth=0.8,
            label="numérique")
    for i, val in enumerate(ana_vals):
        ax5.hlines(val, i - 0.3, i + 0.3, colors="k", linestyles="--", lw=1.8,
                   label="analytique" if i == 0 else "")
    ax5.axhline(0, color="k", lw=0.8)
    ax5.set_ylabel(r"$Q / L_z$  [C/m]")
    ax5.legend(fontsize=11)
    ax5.grid(True, alpha=0.3, axis='y')
    fig5.tight_layout()
    fig5.savefig(f"scan_N1_charges_{MODE}.png", dpi=150)

plt.show()
