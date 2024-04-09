#include <cmath>

#include "pbf_1s0.h"
#include "utils.h"

double Pbf_1s0::Integrand(double r, double E) {
    double T = nscool->GetT(r);
    double Tc;

    if (source == "n") {
        Tc = nscool->GetTcn(r);
    } else if (source == "p") {
        Tc = nscool->GetTcp(r);
    } else {
        return 0;
    }

    double DeltaT = GetDeltaT(T, Tc);
    if (DeltaT <= 0) {
        return 0.;
    }

    double ephi = nscool->GetEphi(r);
    double dvdr = nscool->GetDvdr(r);
    double omega = E * keV2K / ephi;
    return J(omega, T, DeltaT) * dvdr * ephi;
}

double Pbf_1s0::J(double omega, double T, double DeltaT) {
    double argwt = omega / (2.*DeltaT);
    if (argwt <= 1) {
        return 0;
    }

    double fermi = std::exp(omega/(2.*T))+1.;
    fermi = 1./(fermi*fermi);
    return 1./(DeltaT)/keV2K *
            argwt*argwt*argwt / std::sqrt(argwt*argwt-1) * fermi;
}

double Pbf_1s0::GetDeltaT(double T, double Tc) {
    double tau = T/Tc;
    return tau >=1 ? 0. : T * std::sqrt(1.-tau)*(
                            1.456 - 0.157/std::sqrt(tau) + 1.764/tau);
}

void Pbf_1s0::GetBoundaries(double* rmin, double* rmax) {
    if (source == "p") {
        (*rmin) = 0;
        (*rmax) = nscool->Get1s03p2Boundary();
    } else if (source == "n") {
        (*rmin) = nscool->Get1s03p2Boundary();
        (*rmax) = nscool->GetRMax();
    }
}
