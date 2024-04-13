#include <cmath>

#include "pbf_1s0.h"
#include "utils.h"

double Pbf_1s0::dI_dE_dr(double r, double E) {
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
    const double I = Ias(T, DeltaT);
    // Internally, everything is consistently in metres and Kelvin, only the
    // input 'E' is in keV. The 1e-6 * keV2K factor is to convert the final
    // result from erg/(m^3 K sec) to erg/(cm^3 keV sec).
    return J(omega, T, DeltaT, I) * Epsilon(r, T, I) * dvdr * ephi * ephi *
           1e-6 * keV2K;
}

double Pbf_1s0::Epsilon(double r, double T, double I) {
    const double kfn = nscool->GetKfn(r);
    const double kfp = nscool->GetKfp(r);
    const double mstn = nscool->GetMstn(r);
    const double mstp = nscool->GetMstp(r);
    const double T8 = T / 1e8;
    const double gaNN10 = gaNN * 1e10;
    const double kfn168 = kfn / 1.68;
    const double kfp168 = kfp / 1.68;
    // Eq.(S6) of M. Buschmann et al, Phys. Rev. Lett. 128 (2022) 091102
    // [2111.09892]
    const double gamma = 1 / (1. + mstn * kfn168 / 3.);

    double epsilon = 0;
    if (source == "n") {
        // Eq.(S7) of M. Buschmann et al, Phys. Rev. Lett. 128 (2022) 091102
        // [2111.09892]
        // Internally everything is consistently in metres not centimetres.
        epsilon = (4.692e12 * 1e6) * gaNN10 * gaNN10 * kfn168 * kfn168 *
                  kfn168 * T8 * T8 * T8 * T8 * T8 / mstn * gamma * gamma *
                  (I / 0.022);
    } else if (source == "p") {
        // Eq.(S8) of M. Buschmann et al, Phys. Rev. Lett. 128 (2022) 091102
        // [2111.09892]
        // Internally everything is consistently in metres not centimetres.
        epsilon = (4.711e12 * 1e6) * gaNN10 * gaNN10 * kfp168 * kfp168 *
                  kfp168 * T8 * T8 * T8 * T8 * T8 / mstp * gamma * gamma *
                  (I / 0.022);
    }
    return epsilon;
}

inline double Pbf_1s0::Ias(double T, double DeltaT) {
    double z = DeltaT / T;
    if (z >= 15 || z <= 0) {
        return 0;
    }

    // Eq.(22) of A. Sedrakian, Phys. Rev. D93, 065044 (2016), arXiv:1512.07828
    const double a = 0.158151;
    const double c = 0.543166;
    const double h = 0.0535359;
    const double f = M_PI / (4. * c * c);
    return (a * z * z + c * z * z * z * z) * std::sqrt(1 + f * z) *
           std::exp(h - std::sqrt(4. * z * z + h * h));
}

inline double Pbf_1s0::J(double omega, double T, double DeltaT, double I) {
    double argwt = omega / (2. * DeltaT);
    if (argwt <= 1 || I <= 0) {
        return 0;
    }
    double z = DeltaT / T;

    // Eq.(S4) of M. Buschmann et al, Phys. Rev. Lett. 126, 021102 (2021),
    // arXiv:1910.04164
    double N = z * z * z * z * z / I;
    double fermi = Fermi(omega / (2. * T));
    return N / (2. * DeltaT) * argwt * argwt * argwt /
           std::sqrt(argwt * argwt - 1) * fermi * fermi;
}

double Pbf_1s0::GetDeltaT(double T, double Tc) {
    const double tau = T / Tc;
    if (tau >= 1) {
        return 0;
    }
    // Eq.(23) of D. G. Yakovlev et al, Astronomy and Astrophysics 297,
    // 717 (1995).
    return T * std::sqrt(1. - tau) *
           (1.456 - 0.157 / std::sqrt(tau) + 1.764 / tau);
}

void Pbf_1s0::GetBoundaries(double *rmin, double *rmax) {
    if (source == "p") {
        (*rmin) = 0;
        (*rmax) = nscool->Get1s03p2Boundary();
    } else if (source == "n") {
        (*rmin) = nscool->Get1s03p2Boundary();
        (*rmax) = nscool->GetRMax();
    }
}
