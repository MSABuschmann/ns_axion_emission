#include <cmath>
#include <numeric>

#include "pbf_process.h"
#include "utils.h"

#define keV2K 11605000.

double PbfProcess::Integrand(double E, double r) {
    double T = nscool->GetT(r);
    double Tc = nscool->GetTcn(r);
    double DeltaT = GetDeltaT(T, Tc);
    if (DeltaT <= 0) {
        return 0.;
    }
    double ephi = nscool->GetEphi(r);
    double dvdr = nscool->GetDvdr(r);
    double omega = E * keV2K / ephi;
    return J(omega, T, DeltaT) * dvdr * ephi;
}

std::vector<double> PbfProcess::GetSpectrum(std::vector<double>& E_bins,
                                            int N) {
    std::cout << "Compute spectrum with N_r = " << N << " ..." << std::endl;

    double rmin=0, rmax=0;
    if (state == SF_3p2) {
        rmin = 0;
        rmax = nscool->Get1s03p2Boundary();
    } else if (state == SF_1s0) {
        rmin = nscool->Get1s03p2Boundary();
        rmax = nscool->GetRMax();
    }
    std::vector<double> r = CreateVector(rmin, rmax, N);
    std::vector<double> spectrum(E_bins.size(), 0.);
    for (unsigned int i = 0; i < E_bins.size(); ++i) {
        for (int j = 0; j < N; ++j) {
            spectrum[i] += Integrand(E_bins[i], r[j]);
        }
    }
    Normalize(E_bins, spectrum);
    return spectrum;
}

void PbfProcess::Normalize(std::vector<double>& x, std::vector<double>& y) {
    double sum = std::accumulate(y.begin(), y.end(), 0.);
    double norm = 1./(sum*(x.back()-x[0]));
    for (unsigned int i = 0; i < y.size(); ++i) {
        y[i] *= norm;
    }
}

std::vector<double> PbfProcess::GetDeltaT(std::vector<double> T,
                                          std::vector<double> Tc) {
    std::vector<double> DeltaT(T.size(), 0.);
    for(unsigned int i = 0; i < T.size(); ++i) {
        DeltaT[i] = GetDeltaT(T[i], Tc[i]);
    }
    return DeltaT;
}

double PbfProcess_1s0::J(double omega, double T, double DeltaT) {
    double argwt = omega / (2.*DeltaT);
    if (argwt <= 1) {
        return 0;
    }

    double fermi = std::exp(omega/(2.*T))+1.;
    fermi = 1./(fermi*fermi);

    return 1./(DeltaT) *
            argwt*argwt*argwt / std::sqrt(argwt*argwt-1) * fermi;
}

double PbfProcess_1s0::GetDeltaT(double T, double Tc) {
    double tau = T/Tc;
    return tau >=1 ? 0. : T * std::sqrt(1.-tau)*(
                            1.456 - 0.157/std::sqrt(tau) + 1.764/tau);
}
