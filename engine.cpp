#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../../common/ConfigFile.h"

using namespace std;

const double PI = 3.1415926535897932384626433832795028841971e0;

// Resolution d'un systeme d'equations lineaires par elimination de Gauss-Jordan
// (tridiagonal system: diag, lower, upper, rhs all of consistent sizes)
template<class T>
vector<T> solve(const vector<T>& diag,
                const vector<T>& lower,
                const vector<T>& upper,
                const vector<T>& rhs)
{
    vector<T> solution(diag.size());
    vector<T> new_diag(diag);
    vector<T> new_rhs(rhs);

    for (int i = 1; i < (int)diag.size(); ++i) {
        double pivot = lower[i - 1] / new_diag[i - 1];
        new_diag[i] -= pivot * upper[i - 1];
        new_rhs[i]  -= pivot * new_rhs[i - 1];
    }

    solution[diag.size() - 1] = new_rhs[diag.size() - 1] / new_diag[diag.size() - 1];

    for (int i = (int)diag.size() - 2; i >= 0; --i)
        solution[i] = (new_rhs[i] - upper[i] * solution[i + 1]) / new_diag[i];

    return solution;
}

// TODO: Implement the relative permittivity epsilon_r(r).
//       Should allow for a trivial test case (trivial=true) 
double epsilon_r(/* TODO: add arguments */)
{
    return 0.0; // TODO: replace
}

// TODO: Implement the normalised free charge density rho_lib(r) / epsilon_0.
//       Should allow for a trivial test case (trivial=true) 
double rho_lib(/* TODO: add arguments */)
{
    return 0.0; // TODO: replace
}

int main(int argc, char* argv[])
{
    // USAGE: ./engine [configuration-file] [<key=value> ...]

    string inputPath = "trivial.in";
    if (argc > 1)
        inputPath = argv[1];

    ConfigFile configFile(inputPath);
    for (int i = 2; i < argc; ++i)
        configFile.process(argv[i]);

    // Physical parameters
    const double b   = configFile.get<double>("b");   // Inner radius [m]
    const double R   = configFile.get<double>("R");   // Outer radius [m]
    const double V0  = configFile.get<double>("V0");  // Boundary potential at r=R [V]
    const double a0  = configFile.get<double>("a0");  // Free charge density scale [V/m^2]
    const bool trivial = configFile.get<bool>("trivial"); // true: uniform test case

    // Discretisation
    const int N1 = configFile.get<int>("N1"); // Intervals in [0, b]
    const int N2 = configFile.get<int>("N2"); // Intervals in [b, R]

    // Output file prefix
    const string output = configFile.get<string>("output");

    // ---------------------------------------------------------------
    // Build grid
    // ---------------------------------------------------------------
    const int ninters = N1 + N2;         // Total number of intervals
    const int npoints = ninters + 1;     // Total number of grid points
    const double h1 = b / N1;            // Step size in inner region
    const double h2 = (R - b) / N2;      // Step size in outer region

    vector<double> r(npoints);
    // TODO: fill r[0..N1] and r[N1..npoints-1]

    vector<double> h(ninters);           // Interval widths
    vector<double> midPoint(ninters);    // Midpoints of each interval
    // TODO: fill h[i]  and  midPoint[i]

    // ---------------------------------------------------------------
    // Assemble the tridiagonal system  A * phi = rhs
    // ---------------------------------------------------------------
    vector<double> diag(npoints, 0.0);   // Main diagonal
    vector<double> lower(ninters, 0.0);  // Sub-diagonal  (lower[i] links row i+1 to col i)
    vector<double> upper(ninters, 0.0);  // Super-diagonal (upper[i] links row i to col i+1)
    vector<double> rhs(npoints, 0.0);    // Right-hand side

    for (int k = 0; k < ninters; ++k) {
        // TODO: compute alpha_k and beta_k
        //       then add their contributions to diag, lower, upper, and rhs
    }

    // TODO: enforce the Dirichlet BC at r = R

    // ---------------------------------------------------------------
    // Solve the linear system
    // ---------------------------------------------------------------
    vector<double> phi = solve(diag, lower, upper, rhs);

    // ---------------------------------------------------------------
    // Compute electric field E_r and displacement D_r (normalised by eps0)
    // ---------------------------------------------------------------
    vector<double> rmid(ninters);
    vector<double> Er(ninters, 0.0);
    vector<double> Dr(ninters, 0.0);
    for (int k = 0; k < ninters; ++k) {
        rmid[k] = midPoint[k];
        // TODO: compute Er[k] and Dr[k] 
    }

    // ---------------------------------------------------------------
    // Compute div(D_r)/eps0 and compare to rho_lib/eps0
    // using finite differences on the midpoint values
    // ---------------------------------------------------------------
    vector<double> rmidmid(ninters - 1);
    vector<double> div_Dr(ninters - 1, 0.0);
    vector<double> rho_at_midmid(ninters - 1, 0.0);
    for (int k = 0; k < ninters - 1; ++k) {
        rmidmid[k] = 0.5 * (rmid[k] + rmid[k + 1]);
        // TODO: compute div_Dr[k] and rho_at_midmid[k]
    }

    // ---------------------------------------------------------------
    // Write output files
    // ---------------------------------------------------------------
    {
        // 1. Electric potential: columns  r  phi
        ofstream ofs(output + "_phi.out");
        ofs.precision(15);
        for (int i = 0; i < npoints; ++i)
            ofs << r[i] << " " << phi[i] << "\n";
    }
    {
        // 2. Electric field and displacement: columns  r_mid  E_r  D_r/eps0
        ofstream ofs(output + "_ErDr.out");
        ofs.precision(15);
        for (int k = 0; k < ninters; ++k)
            ofs << rmid[k] << " " << Er[k] << " " << Dr[k] << "\n";
    }
    {
        // 3. Divergence check: columns  r_midmid  div(D_r)/eps0  rho_lib/eps0
        ofstream ofs(output + "_divD_rho.out");
        ofs.precision(15);
        for (int k = 0; k < ninters - 1; ++k)
            ofs << rmidmid[k] << " " << div_Dr[k]
                << " " << rho_at_midmid[k] << "\n";
    }

    return 0;
}
