#include <cmath>

#include "pbf_3p2.h"
#include "utils.h"

double Pbf_3p2::Integrand(double r, double E) {
    double T = nscool->GetT(r);
    double Tc = nscool->GetTcn(r);
    double DeltaT = GetDeltaT(T, Tc);
    double DeltaT0 = GetDeltaT(T, Tc);
    if (DeltaT <= 0) {
        return 0.;
    }

    double ephi = nscool->GetEphi(r);
    double dvdr = nscool->GetDvdr(r);
    double omega = E * keV2K / ephi;
    const double I = Ian(T, DeltaT);
    const double I0 = Ian(T, DeltaT0);
    return J(omega, T, DeltaT, I0) * Epsilon(r, T, I) * dvdr * ephi;
}

double Pbf_3p2::Epsilon(double r, double T, double I) {
    const double kfn = nscool->GetKfn(r);
    const double mstn = nscool->GetMstn(r);
    const double T8 = T / 1e8;
    const double gann10 = gann * 1e10;
    const double kfn168 = kfn / 1.68;
    const double gamma = 1 / (1. + mstn * kfn168 / 3.);
    const double epsilon = 3.769e13 * gann10 * gann10 * kfn168 * T8 * T8 * T8 *
                           T8 * T8 / mstn * gamma * gamma * (I / 0.022);
    return epsilon;
}

inline double Pbf_3p2::Ian(double T, double DeltaT) {
    double z = DeltaT / T;
    if (source == "A") {
        return (0.316302 * z * z + 0.0255728 * z * z * z * z) *
               (1. + 2.22858 * z * z) /
               (1. + 0.856577 * z * z + 0.000449543 * z * z * z * z) *
               std::exp(2.22569 - std::sqrt(4. * z * z + 4.9536959));
    } else if (source == "B") {
        return (0.105433 * z * z - 0.000271463 * z * z * z * z) *
               std::sqrt(1. + 0.006347022 * z * z) /
               (1 - 0.043745 * z * z + 0.0216661 * z * z * z * z);
    }
    return 0;
}

inline double Pbf_3p2::J(double omega, double T, double DeltaT, double I) {
    double argwt = omega / (2. * DeltaT);
    if (argwt <= 1) {
        return 0;
    }

    double z = T / DeltaT;
    double N = I / (T * T * T * T * T);

    double fermi = std::exp(omega / (2. * T)) + 1.;
    fermi = 1. / (fermi * fermi);
    return N / (DeltaT * 2.) / keV2K * argwt * argwt * argwt /
           std::sqrt(argwt * argwt - 1) * fermi;
}

double Pbf_3p2::GetDeltaT(double T, double Tc) { return 0; }

double Pbf_3p2::GetDeltaT(double T, double Tc, double theta) {
    double tau = T / Tc;
    return tau >= 1 ? 0.
                    : T * std::sqrt(1. - tau) *
                          (1.456 - 0.157 / std::sqrt(tau) + 1.764 / tau);
}

void Pbf_3p2::GetBoundaries(double *rmin, double *rmax) {
    (*rmin) = 0;
    (*rmax) = nscool->Get1s03p2Boundary();
}
