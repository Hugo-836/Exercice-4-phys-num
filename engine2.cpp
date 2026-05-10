#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "common/ConfigFile.h"

using namespace std;

// epsilon_r(r) = 2            if r < b
//              = 6 + 4*((r-b)/(R-b))^2  if b <= r <= R  (jump 2->6 at b)
double epsilon_r(double r, double b, double R)
{
    if (r < b)
        return 2.0;
    else
        return 6.0 + 4.0 * pow((r - b) / (R - b), 2);
}

// rho_lib(r)/eps0 = a0*(1 - r/b)^2  if r < b
//                 = 0                if b <= r <= R
double rho_lib(double r, double b, double a0)
{
    if (r < b)
        return a0 * pow(1.0 - r / b, 2);
    else
        return 0.0;
}

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

int main(int argc, char* argv[])
{
    string inputPath = "trivial.in";
    if (argc > 1)
        inputPath = argv[1];

    ConfigFile configFile(inputPath);
    for (int i = 2; i < argc; ++i)
        configFile.process(argv[i]);

    const double b   = configFile.get<double>("b");
    const double R   = configFile.get<double>("R");
    const double V0  = configFile.get<double>("V0");
    const double a0  = configFile.get<double>("a0");
    const int    N1  = configFile.get<int>("N1");
    const int    N2  = configFile.get<int>("N2");
    const string output = configFile.get<string>("output");

    // ---------------------------------------------------------------
    // Build grid
    // ---------------------------------------------------------------
    const int ninters = N1 + N2;
    const int npoints = ninters + 1;
    const double h1 = b / N1;
    const double h2 = (R - b) / N2;

    vector<double> r(npoints);
    r[0] = 0.0;
    for (int i = 1; i <= N1; ++i)
        r[i] = i * h1;
    for (int i = 1; i <= N2; ++i)
        r[N1 + i] = b + i * h2;

    vector<double> h(ninters);
    vector<double> midPoint(ninters);
    for (int i = 0; i < N1; ++i) {
        h[i]        = h1;
        midPoint[i] = (r[i] + r[i+1]) / 2.0;
    }
    for (int i = N1; i < ninters; ++i) {
        h[i]        = h2;
        midPoint[i] = (r[i] + r[i+1]) / 2.0;
    }

    // ---------------------------------------------------------------
    // Assemble the tridiagonal system  A * phi = rhs
    // ---------------------------------------------------------------
    vector<double> diag(npoints, 0.0);
    vector<double> lower(ninters, 0.0);
    vector<double> upper(ninters, 0.0);
    vector<double> rhs(npoints, 0.0);

    for (int k = 0; k < ninters; ++k) {
        double eps_mid = epsilon_r(midPoint[k], b, R);
        double rho_mid = rho_lib(midPoint[k], b, a0);
        double alpha_k = midPoint[k] * eps_mid / h[k];
        double beta_k  = midPoint[k] * rho_mid * h[k] / 2.0;
        diag[k]   += alpha_k;
        upper[k]  += -alpha_k;
        diag[k+1] += alpha_k;
        lower[k]  += -alpha_k;
        rhs[k]    += beta_k;
        rhs[k+1]  += beta_k;
    }

    // Dirichlet BC at r = R
    int last = npoints - 1;
    diag[last]      = 1.0;
    rhs[last]       = V0;
    upper[last - 1] = 0.0;
    lower[last - 1] = 0.0;

    // ---------------------------------------------------------------
    // Solve
    // ---------------------------------------------------------------
    vector<double> phi = solve(diag, lower, upper, rhs);

    // ---------------------------------------------------------------
    // Compute E_r and D_r/eps0
    // ---------------------------------------------------------------
    vector<double> rmid(ninters);
    vector<double> Er(ninters, 0.0);
    vector<double> Dr(ninters, 0.0);
    for (int k = 0; k < ninters; ++k) {
        rmid[k] = midPoint[k];
        Er[k]   = -(phi[k+1] - phi[k]) / h[k];
        Dr[k]   = epsilon_r(midPoint[k], b, R) * Er[k];
    }

    // ---------------------------------------------------------------
    // Compute div(D_r)/eps0 and rho_lib/eps0
    // ---------------------------------------------------------------
    vector<double> rmidmid(ninters - 1);
    vector<double> div_Dr(ninters - 1, 0.0);
    vector<double> rho_at_midmid(ninters - 1, 0.0);
    for (int k = 0; k < ninters - 1; ++k) {
        rmidmid[k]       = 0.5 * (rmid[k] + rmid[k + 1]);
        div_Dr[k]        = (rmid[k+1] * Dr[k+1] - rmid[k] * Dr[k])
                           / (rmidmid[k] * (rmid[k+1] - rmid[k]));
        rho_at_midmid[k] = rho_lib(rmidmid[k], b, a0);
    }

    // ---------------------------------------------------------------
    // Write output files
    // ---------------------------------------------------------------
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
