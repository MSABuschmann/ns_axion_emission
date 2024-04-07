#include <cmath>

#include "1s0.h"

double PbfProcess_1s0::Integrand(double r, double E) {
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

double PbfProcess_1s0::J(double omega, double T, double DeltaT) {
    double argwt = omega / (2.*DeltaT);
    if (argwt <= 1) {
        return 0;
    }

    double fermi = std::exp(omega/(2.*T))+1.;
    fermi = 1./(fermi*fermi);
    return 1./(DeltaT)/keV2K *
            argwt*argwt*argwt / std::sqrt(argwt*argwt-1) * fermi;
}

double PbfProcess_1s0::GetDeltaT(double T, double Tc) {
    double tau = T/Tc;
    return tau >=1 ? 0. : T * std::sqrt(1.-tau)*(
                            1.456 - 0.157/std::sqrt(tau) + 1.764/tau);
}
