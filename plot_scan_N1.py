import numpy as np
import matplotlib.pyplot as plt
import os

scan_dir = "Scan_N1_electrostatics_b_0.02_R_0.1_trivial_false"
prefix   = "electrostatics_b_0.02_R_0.1_trivial_false_N1_"
N_values = [2, 4, 8, 16, 32, 64, 128, 256]

b  = 0.02
R  = 0.1

# suffixe ajouté aux noms des figures pour ne pas écraser le cas trivial
fig_suffix = "_nontrivial"

# ── couleurs dégradées ─────────────────────────────────────────────────────
cmap   = plt.cm.viridis
colors = [cmap(i / (len(N_values) - 1)) for i in range(len(N_values))]

# ══════════════════════════════════════════════════════════════════════════
# Figure 1 : phi(r) pour chaque N1
# ══════════════════════════════════════════════════════════════════════════
fig1, ax1 = plt.subplots(figsize=(7, 5))

for N, col in zip(N_values, colors):
    fname = os.path.join(scan_dir, f"{prefix}{N}_phi.out")
    data  = np.loadtxt(fname)
    r, phi = data[:, 0], data[:, 1]
    ax1.plot(r, phi, color=col, lw=1.2, label=f"N={N}")

ax1.axvline(b, color="gray", ls="--", lw=0.8, label=f"b = {b}")
ax1.set_xlabel("r [m]")
ax1.set_ylabel(r"$\phi(r)$ [V]")
ax1.set_title("Potentiel électrostatique — scan N1")
ax1.legend(fontsize=8, ncol=2)
ax1.grid(True, alpha=0.3)
fig1.tight_layout()
fig1.savefig(f"scan_N1_phi{fig_suffix}.png", dpi=150)
print("Saved: scan_N1_phi.png")

# ══════════════════════════════════════════════════════════════════════════
# Figure 2 : Er et Dr pour le plus grand N
# ══════════════════════════════════════════════════════════════════════════
fig2, axes = plt.subplots(1, 2, figsize=(11, 4))

for N, col in zip(N_values, colors):
    fname = os.path.join(scan_dir, f"{prefix}{N}_ErDr.out")
    data  = np.loadtxt(fname)
    r_mid, Er, Dr = data[:, 0], data[:, 1], data[:, 2]
    axes[0].plot(r_mid, Er, color=col, lw=1.2, label=f"N={N}")
    axes[1].plot(r_mid, Dr, color=col, lw=1.2, label=f"N={N}")

for ax, ylabel, title in zip(axes,
        [r"$E_r$ [V/m]", r"$D_r$ [C/m²]"],
        ["Champ électrique Er", "Déplacement électrique Dr"]):
    ax.axvline(b, color="gray", ls="--", lw=0.8, label=f"b = {b}")
    ax.set_xlabel("r [m]")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3)

fig2.tight_layout()
fig2.savefig(f"scan_N1_ErDr{fig_suffix}.png", dpi=150)
print("Saved: scan_N1_ErDr.png")

# ══════════════════════════════════════════════════════════════════════════
# Figure 3 : div D / rho (doit être ≈ 1 partout)
# ══════════════════════════════════════════════════════════════════════════
fig3, ax3 = plt.subplots(figsize=(7, 5))

for N, col in zip(N_values, colors):
    fname = os.path.join(scan_dir, f"{prefix}{N}_divD_rho.out")
    data  = np.loadtxt(fname)
    r_mid, ratio = data[:, 0], data[:, 1]
    ax3.plot(r_mid, ratio, color=col, lw=1.2, label=f"N={N}")

ax3.axhline(1.0, color="k", ls="--", lw=0.8, label="valeur exacte = 1")
ax3.axvline(b, color="gray", ls=":", lw=0.8, label=f"b = {b}")
ax3.set_xlabel("r [m]")
ax3.set_ylabel(r"$(\nabla \cdot \mathbf{D})\,/\,\rho$")
ax3.set_title(r"Vérification $\nabla \cdot D = \rho$ — scan N1")
ax3.legend(fontsize=8, ncol=2)
ax3.grid(True, alpha=0.3)
fig3.tight_layout()
fig3.savefig(f"scan_N1_divD_rho{fig_suffix}.png", dpi=150)
print("Saved: scan_N1_divD_rho.png")

# ══════════════════════════════════════════════════════════════════════════
# Figure 4 : convergence de phi(0) et erreur max sur divD/rho
# ══════════════════════════════════════════════════════════════════════════
phi0_vals      = []
err_divD_max   = []
err_phi_rms    = []

# solution exacte de référence (N=256 comme proxy, ou valeur analytique)
fname_ref = os.path.join(scan_dir, f"{prefix}256_phi.out")
ref_phi   = np.loadtxt(fname_ref)
phi0_exact = ref_phi[0, 1]      # phi(r=0) exact ≈ 0.005

for N in N_values:
    fname = os.path.join(scan_dir, f"{prefix}{N}_phi.out")
    data  = np.loadtxt(fname)
    phi0_vals.append(data[0, 1])

    fname_d = os.path.join(scan_dir, f"{prefix}{N}_divD_rho.out")
    data_d  = np.loadtxt(fname_d)
    err_divD_max.append(np.max(np.abs(data_d[:, 1] - 1.0)))

N_arr  = np.array(N_values)
phi0   = np.array(phi0_vals)
err_phi = np.abs(phi0 - phi0_exact)

fig4, (ax4a, ax4b) = plt.subplots(1, 2, figsize=(11, 4))

# --- phi(0) vs N ---
ax4a.loglog(N_arr, err_phi + 1e-16, "o-", color="steelblue", label=r"|$\phi(0) - \phi_{exact}$|")
# guide pente -2
x_g = np.array([N_arr[0], N_arr[-2]])
ax4a.loglog(x_g, err_phi[0] * (x_g / N_arr[0])**(-2),
            "k--", lw=0.8, label=r"$\propto N^{-2}$")
ax4a.set_xlabel("N  (= N1 = N2)")
ax4a.set_ylabel(r"Erreur sur $\phi(0)$")
ax4a.set_title(r"Convergence de $\phi(0)$")
ax4a.legend(fontsize=9)
ax4a.grid(True, which="both", alpha=0.3)

# --- erreur max divD/rho vs N ---
ax4b.loglog(N_arr, err_divD_max, "s-", color="tomato",
            label=r"$\max|\nabla\cdot D/\rho - 1|$")
x_g2 = np.array([N_arr[0], N_arr[-2]])
ax4b.loglog(x_g2, err_divD_max[0] * (x_g2 / N_arr[0])**(-2),
            "k--", lw=0.8, label=r"$\propto N^{-2}$")
ax4b.set_xlabel("N  (= N1 = N2)")
ax4b.set_ylabel(r"Erreur max $(\nabla\cdot D)/\rho$")
ax4b.set_title(r"Convergence $\nabla\cdot D = \rho$")
ax4b.legend(fontsize=9)
ax4b.grid(True, which="both", alpha=0.3)

fig4.tight_layout()
fig4.savefig(f"scan_N1_convergence{fig_suffix}.png", dpi=150)
print("Saved: scan_N1_convergence.png")

plt.show()
