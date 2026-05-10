#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "common/ConfigFile.h"

using namespace std;

const double PI = 3.1415926535897932384626433832795028841971e0;

// ── Profils physiques (identiques à engine.cpp) ────────────────────────────

double epsilon_r(bool trivial, double r, double b, double R)
{
    if (trivial) return 1.0;
    if (r <= b)  return 1.0;
    return 3.0 + 6.0 * ((r - b) / (R - b));
}

double rho_lib(bool trivial, double r, double b, double R, double a0)
{
    if (trivial)  return 1.0;
    if (r <= b)   return a0 * sin(PI * r / b);
    return 0.0;
}

// ── Membre de droite de l'ODE ──────────────────────────────────────────────
//  y = { phi, F }  avec  F = r · ε_r · dφ/dr
//  Retourne dy/dr = { dφ/dr, dF/dr }

array<double,2> computeRhs(double r, const array<double,2>& y,
                            bool trivial, double b, double R, double a0)
{
    const double eps = epsilon_r(trivial, r, b, R);
    const double rho = rho_lib  (trivial, r, b, R, a0);

    // En r=0 : F=0 et r=0 → limite 0/0 = 0 par continuité (symétrie cylindrique)
    const double dphi_dr = (r < 1e-30) ? 0.0 : y[1] / (r * eps);
    const double dF_dr   = -r * rho;

    return {dphi_dr, dF_dr};
}

// ── Un pas RK4 classique (longueur h, départ en r) ────────────────────────

array<double,2> rk4Step(double r, double h, const array<double,2>& y,
                         bool trivial, double b, double R, double a0)
{
    auto k1 = computeRhs(r,         y,                        trivial, b, R, a0);

    array<double,2> y2 = {y[0] + 0.5*h*k1[0], y[1] + 0.5*h*k1[1]};
    auto k2 = computeRhs(r + 0.5*h, y2,                      trivial, b, R, a0);

    array<double,2> y3 = {y[0] + 0.5*h*k2[0], y[1] + 0.5*h*k2[1]};
    auto k3 = computeRhs(r + 0.5*h, y3,                      trivial, b, R, a0);

    array<double,2> y4 = {y[0] + h*k3[0], y[1] + h*k3[1]};
    auto k4 = computeRhs(r + h,      y4,                     trivial, b, R, a0);

    return {
        y[0] + (h / 6.0) * (k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0]),
        y[1] + (h / 6.0) * (k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1])
    };
}

// =============================================================================
int main(int argc, char* argv[])
{
    string inputPath = "trivial.in";
    if (argc > 1) inputPath = argv[1];

    ConfigFile configFile(inputPath);
    for (int i = 2; i < argc; ++i)
        configFile.process(argv[i]);

    // ── Paramètres physiques ─────────────────────────────────────────────────
    const double b       = configFile.get<double>("b");
    const double R       = configFile.get<double>("R");
    const double V0      = configFile.get<double>("V0");
    const double a0      = configFile.get<double>("a0");
    const bool   trivial = configFile.get<bool>  ("trivial");

    // ── Paramètres numériques ────────────────────────────────────────────────
    const int    N1      = configFile.get<int>   ("N1");
    const int    N2      = configFile.get<int>   ("N2");
    const string output  = configFile.get<string>("output");

    // ── Construction de la grille (identique à engine.cpp) ──────────────────
    const int    ninters = N1 + N2;
    const int    npoints = ninters + 1;
    const double h1      = b / N1;
    const double h2      = (N2 > 0) ? (R - b) / N2 : 0.0;

    vector<double> r(npoints);
    r[0] = 0.0;
    for (int i = 1; i <= N1; ++i) r[i]      = i * h1;
    for (int i = 1; i <= N2; ++i) r[N1 + i] = b + i * h2;

    vector<double> h(ninters), midPoint(ninters);
    for (int i = 0; i < N1; ++i) {
        h[i]        = h1;
        midPoint[i] = 0.5 * (r[i] + r[i+1]);
    }
    for (int i = N1; i < ninters; ++i) {
        h[i]        = h2;
        midPoint[i] = 0.5 * (r[i] + r[i+1]);
    }

    // ── Intégration RK4 de r=0 à r=R ────────────────────────────────────────
    //  Conditions initiales : φ(0)=0 (arbitraire), F(0)=0 (symétrie)
    vector<array<double,2>> y(npoints);
    y[0] = {0.0, 0.0};

    for (int i = 0; i < ninters; ++i)
        y[i+1] = rk4Step(r[i], h[i], y[i], trivial, b, R, a0);

    // ── Correction affine pour satisfaire φ(R) = V0 ─────────────────────────
    //  L'équation est linéaire : ajouter (V0 - φ_calc(R)) à tous les φ est exact.
    const double correction = V0 - y[npoints - 1][0];

    vector<double> phi(npoints);
    for (int i = 0; i < npoints; ++i)
        phi[i] = y[i][0] + correction;

    // ── Champ électrique Er et déplacement Dr aux points milieux ─────────────
    vector<double> rmid(ninters), Er(ninters), Dr(ninters);
    for (int k = 0; k < ninters; ++k) {
        rmid[k] = midPoint[k];
        Er[k]   = -(phi[k+1] - phi[k]) / h[k];
        Dr[k]   = epsilon_r(trivial, midPoint[k], b, R) * Er[k];
    }

    // ── Divergence de D et densité de charge aux points milieu-milieux ───────
    vector<double> rmidmid(ninters - 1), div_Dr(ninters - 1), rho_at_midmid(ninters - 1);
    for (int k = 0; k < ninters - 1; ++k) {
        rmidmid[k] = 0.5 * (rmid[k] + rmid[k+1]);
        // (1/r) d(r·Dr)/dr discrétisé au point rmidmid[k]
        div_Dr[k] = (rmid[k+1] * Dr[k+1] - rmid[k] * Dr[k])
                    / (rmidmid[k] * (rmid[k+1] - rmid[k]));
        rho_at_midmid[k] = rho_lib(trivial, rmidmid[k], b, R, a0);
    }

    // ── Écriture des fichiers de sortie ──────────────────────────────────────
    {
        ofstream ofs(output + "_phi.out");
        ofs.precision(15);
        for (int i = 0; i < npoints; ++i)
            ofs << r[i] << " " << phi[i] << "\n";
    }
    {
        ofstream ofs(output + "_ErDr.out");
        ofs.precision(15);
        for (int k = 0; k < ninters; ++k)
            ofs << rmid[k] << " " << Er[k] << " " << Dr[k] << "\n";
    }
    {
        ofstream ofs(output + "_divD_rho.out");
        ofs.precision(15);
        for (int k = 0; k < ninters - 1; ++k)
            ofs << rmidmid[k] << " " << div_Dr[k]
                << " " << rho_at_midmid[k] << "\n";
    }

    return 0;
}
