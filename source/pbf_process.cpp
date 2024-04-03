#include <cmath>

#include "pbf_process.h"

#define keV2K 11605000.

std::vector<double> PbfProcess::GetSpectrum(std::vector<double>& E_bins) {
    const std::vector<double> T ;//= nscool->GetT();
    const std::vector<double> Tc ;//= nscool->GetTc();
    const std::vector<double> ephi ;//= nscool->GetEphi();
    const std::vector<double> DeltaT = GetDeltaT(T, Tc);

    std::vector<double> spectrum(E_bins.size(), 0.);
    for (unsigned int i = 0; i < E_bins.size(); ++i) {
        for (unsigned int j = 0; j < T.size(); ++j) {
            double omega = E_bins[i] * keV2K / ephi[j];
            spectrum[i] += J(omega, T[j], DeltaT[j]);
        }
    }
    return spectrum;
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

    return keV2K/(2.*DeltaT) *
            argwt*argwt*argwt / std::sqrt(argwt*argwt-1) * fermi;
}

double PbfProcess_1s0::GetDeltaT(double T, double Tc) {
    double tau = T/Tc;
    return tau >=1 ? 0. : T * std::sqrt(1.-tau)*(
                            1.456 - 0.157/std::sqrt(tau) + 1.764/tau);
}
