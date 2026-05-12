import numpy as np
import subprocess
import os
import matplotlib.pyplot as plt
from scipy.special import sici

executable = './engine_RK4'
N_values   = np.array([64])                    # N = 64  (phi/Er plots)
N_conv     = 2**np.arange(1, 9)               # N = 2, 4, ..., 256  (convergence)

# ── Définition des 3 cas ────────────────────────────────────────────────────
def _Er_cylindre(r, a0, R):
    """Analytical E_r for full cylinder (4.1.c): eps_r=1, rho=eps0*a0*sin(pi*r/R)."""
    r = np.asarray(r, dtype=float)
    result = np.where(r > 1e-30,
                      a0*(R/np.pi)**2 * (np.sin(np.pi*r/R)/r - (np.pi/R)*np.cos(np.pi*r/R)),
                      0.0)
    return result

def _phi_cylindre(r, a0, R):
    """Analytical phi for full cylinder (4.1.c)."""
    r = np.asarray(r, dtype=float)
    Si_pi   = sici(np.pi)[0]
    Si_pirR = sici(np.pi * r / R)[0]
    return a0*(R/np.pi)**2 * (np.sin(np.pi*r/R) + Si_pi - Si_pirR)

cases = [
    {
        "label"    : "trivial",
        "params"   : dict(b=0.05, R=0.1, V0=0, a0=1,   trivial='true',  N2='N'),
        "phi_exact": lambda N, R=0.1: R**2 / 4,
        "phi_ana"  : lambda r, R=0.1: (R**2 - np.asarray(r)**2) / 4,
        "Er_ana"   : lambda r: np.asarray(r) / 2,
    },
    {
        "label"    : "cylindre_plein",
        "params"   : dict(b=0.1,  R=0.1, V0=0, a0=1e4, trivial='false', N2=0),
        "phi_exact": None,
        "phi_ana"  : lambda r, a0=1e4, R=0.1: _phi_cylindre(r, a0, R),
        "Er_ana"   : lambda r, a0=1e4, R=0.1: _Er_cylindre(r, a0, R),
    },
]

# ── Couleurs communes ────────────────────────────────────────────────────────
cmap   = plt.cm.viridis
colors = [cmap(i / max(len(N_values) - 1, 1)) for i in range(len(N_values))]

for case in cases:
    label  = case["label"]
    params = case["params"]
    b = params["b"];  R = params["R"];  trivial_str = params["trivial"]

    outstr = f"RK4_b_{b:.2g}_R_{R:.2g}_trivial_{trivial_str}"
    outdir = f"Scan_N1_{outstr}"
    os.makedirs(outdir, exist_ok=True)
    print(f"\n{'='*60}")
    print(f"Cas : {label}  →  {outdir}")
    print(f"{'='*60}")

    N_all = np.unique(np.concatenate([N_values, N_conv])).astype(int)

    # ── Scan ─────────────────────────────────────────────────────────────────
    for N in N_all:
        p = {k: v for k, v in params.items()}
        p["N1"] = int(N)
        if p["N2"] == 'N':
            p["N2"] = int(N)

        out_prefix = os.path.join(outdir, f"{outstr}_N1_{int(N)}")
        param_str  = " ".join(f"{k}={v}" for k, v in p.items())
        cmd = f"{executable} trivial.in {param_str} output={out_prefix}"
        print(cmd)
        subprocess.run(cmd, shell=True)

    # ── Figure 1 : phi(r) ────────────────────────────────────────────────────
    fig1, ax1 = plt.subplots(figsize=(7, 5))
    for N, col in zip(N_values, colors):
        fname = os.path.join(outdir, f"{outstr}_N1_{int(N)}_phi.out")
        data  = np.loadtxt(fname)
        ax1.plot(data[:, 0], data[:, 1], color=col, lw=1.2, label=f"N={int(N)}")
    if case["phi_ana"] is not None:
        r_ana = np.linspace(0, R, 500)
        ax1.plot(r_ana, case["phi_ana"](r_ana), "k--", lw=1.5, label="analytical")
    ax1.axvline(b, color="gray", ls=":", lw=0.8, label=f"b = {b}")
    ax1.set_xlabel("r [m]")
    ax1.set_ylabel(r"$\phi(r)$ [V]")
    ax1.legend(fontsize=13);  ax1.grid(True, alpha=0.3)
    fig1.tight_layout()
    fig1.savefig(f"scan_N1_phi_{label}_RK4.png", dpi=150)
    print(f"  Saved: scan_N1_phi_{label}_RK4.png")
    plt.close(fig1)

    # ── Figure 2 : Er ────────────────────────────────────────────────────────
    fig2, ax2 = plt.subplots(figsize=(7, 5))
    for N, col in zip(N_values, colors):
        fname = os.path.join(outdir, f"{outstr}_N1_{int(N)}_ErDr.out")
        data  = np.loadtxt(fname)
        ax2.plot(data[:, 0], data[:, 1], color=col, lw=1.2, label=f"N={int(N)}")
    if case["Er_ana"] is not None:
        r_ana = np.linspace(0, R, 500)
        ax2.plot(r_ana, case["Er_ana"](r_ana), "k--", lw=1.5, label="analytical")
    ax2.axvline(b, color="gray", ls=":", lw=0.8, label=f"b={b}")
    ax2.set_xlabel("r [m]");  ax2.set_ylabel(r"$E_r$ [V/m]")
    ax2.legend(fontsize=13);  ax2.grid(True, alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(f"scan_N1_Er_{label}_RK4.png", dpi=150)
    plt.close(fig2)

    # ── Figure 3 : convergence  |phi(r) - phi_ana(r)|  vs N ─────────────────
    phi_ana_0  = float(np.atleast_1d(case["phi_ana"](np.array([0.0])))[0])
    phi_ana_R2 = float(np.atleast_1d(case["phi_ana"](np.array([R / 2])))[0])

    err_r0  = []
    err_rR2 = []
    for N in N_conv:
        fname = os.path.join(outdir, f"{outstr}_N1_{int(N)}_phi.out")
        data  = np.loadtxt(fname)
        phi_r0  = data[0, 1]
        phi_rR2 = np.interp(R / 2, data[:, 0], data[:, 1])
        err_r0 .append(abs(phi_r0  - phi_ana_0 ))
        err_rR2.append(abs(phi_rR2 - phi_ana_R2))

    N_arr = np.array(N_conv, dtype=float)

    # régression log-log sur phi(0) (sans le premier point N=2)
    log_N   = np.log(N_arr[1:])
    log_err = np.log(np.array(err_r0[1:]) + 1e-16)
    slope0, intercept0 = np.polyfit(log_N, log_err, 1)
    err_fit0 = np.exp(intercept0) * N_arr**slope0

    fig3, ax3 = plt.subplots(figsize=(7, 5))
    ax3.loglog(N_arr, np.array(err_r0)  + 1e-16, "o-", color="steelblue", ms=7,
               label=r"$|\phi(0) - \phi_{\rm ana}(0)|$")
    ax3.loglog(N_arr, np.array(err_rR2) + 1e-16, "s-", color="tomato",    ms=7,
               label=r"$|\phi(R/2) - \phi_{\rm ana}(R/2)|$")
    ax3.loglog(N_arr, err_fit0, "b--", lw=1.2,
               label=fr"régression $\phi(0)$: $N^{{{slope0:.2f}}}$")
    ref0 = err_r0[0] + 1e-16
    ax3.loglog(N_arr, ref0 * (N_arr / N_arr[0])**(-4),
               "k--", lw=0.8, label=r"$\propto N^{-4}$")
    ax3.set_xlabel("N")
    ax3.set_ylabel(r"Error on $\phi$")
    ax3.legend(fontsize=13)
    ax3.grid(True, which="both", alpha=0.3)
    fig3.tight_layout()
    fig3.savefig(f"scan_N1_convergence_{label}_RK4.png", dpi=150)
    print(f"  Saved: scan_N1_convergence_{label}_RK4.png")
    plt.close(fig3)

