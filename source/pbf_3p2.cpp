#include <cmath>
#include <gsl/gsl_integration.h>

#include "pbf_3p2.h"
#include "utils.h"

double Pbf_3p2::Integrand(double r, double E) {
    const double T = nscool->GetT(r);
    const double Tc = nscool->GetTcn(r);
    const double DeltaT0 = GetDeltaT(T, Tc);
    const double I0 = Ian(T, DeltaT0);
    const double ephi = nscool->GetEphi(r);
    const double dvdr = nscool->GetDvdr(r);
    const double omega = E * keV2K / ephi;

    double JTimesEpsilon = 0;
    GslThetaIntegrationParams gsl_params;
    gsl_params.pthis = this;
    gsl_params.r = r;
    gsl_params.omega = omega;
    gsl_params.T = T;
    gsl_params.Tc = Tc;
    gsl_params.I0 = I0;
    gsl_function F;
    F.function = &Pbf_3p2::JTimesEpsilonIntegrand;
    F.params = &gsl_params;
    size_t neval;
    gsl_integration_romberg_workspace *w = gsl_integration_romberg_alloc(15);
    gsl_integration_romberg(&F, -1, 1, 0, 1e-3, &JTimesEpsilon, &neval, w);

    // Internally, everything is consistently in metres and Kelvin, only the
    // input 'E' is in keV. The 1e-6 * keV2K factor is to convert the final
    // result from erg/(m^3 K sec) to erg/(cm^3 keV sec).
    return JTimesEpsilon * dvdr * ephi * ephi * 1e-6 * keV2K;
}

double Pbf_3p2::JTimesEpsilonIntegrand(double cos_theta, void *params) {
    GslThetaIntegrationParams *gsl_params =
        static_cast<GslThetaIntegrationParams *>(params);
    const double r = gsl_params->r;
    const double omega = gsl_params->omega;
    const double T = gsl_params->T;
    const double Tc = gsl_params->Tc;
    const double DeltaT = gsl_params->pthis->GetDeltaT(T, Tc, cos_theta);
    if (DeltaT <= 0) {
        return 0;
    }

    const double I = gsl_params->pthis->Ian(T, DeltaT);
    const double I0 = gsl_params->I0;

    if (gsl_params->pthis->source == "B" && false) {
#pragma omp critical
        {
            std::cout << "B: " << DeltaT << " " << I << " " << I0 << " "
                      << gsl_params->pthis->J(omega, T, DeltaT, I0) << " "
                      << gsl_params->pthis->Epsilon(r, T, I) << std::endl;
        }
    }

    return gsl_params->pthis->J(omega, T, DeltaT, I0) *
           gsl_params->pthis->Epsilon(r, T, I);
}

double Pbf_3p2::Epsilon(double r, double T, double I) {
    const double kfn = nscool->GetKfn(r);
    const double mstn = nscool->GetMstn(r);
    const double T8 = T / 1e8;
    const double gann10 = gann * 1e10;
    const double kfn168 = kfn / 1.68;
    // Eq.(S6) of M. Buschmann et al, Phys. Rev. Lett. 128 (2022) 091102
    // [2111.09892]
    const double gamma = 1 / (1. + mstn * kfn168 / 3.);

    // Eq.(S7) of M. Buschmann et al, Phys. Rev. Lett. 128 (2022) 091102
    // [2111.09892]
    // Internally everything is consistently in metres not centimetres.
    const double epsilon = (3.769e13 * 1e6) * gann10 * gann10 * kfn168 * T8 *
                           T8 * T8 * T8 * T8 * mstn * gamma * gamma *
                           (I / 0.022);
    return epsilon;
}

inline double Pbf_3p2::J(double omega, double T, double DeltaT, double I0) {
    double argwt = omega / (2. * DeltaT);
    if (argwt <= 1 || DeltaT <= 0 || I0 <= 0) {
        return 0;
    }

    // Eq.(S7) of M. Buschmann et al, Phys. Rev. Lett. 126, 021102 (2021),
    // arXiv:1910.04164
    double N = 1. / (T * T * T * T * T * I0);
    double fermi = Fermi(omega / (2. * T));
    return 1. / 4. * DeltaT * DeltaT * DeltaT * DeltaT * N * argwt * argwt *
           argwt / std::sqrt(argwt * argwt - 1) * fermi * fermi;
}

inline double Pbf_3p2::Ian(double T, double DeltaT) {
    const double z = DeltaT / T;
    if (z >= 15 || z <= 0) {
        return 0;
    }

    if (source == "A") {
        // Eq.(23) of A. Sedrakian, Phys. Rev. D93, 065044 (2016),
        // arXiv:1512.07828
        const double a = 2. * 0.158151;
        const double b = 0.856577;
        const double c = 0.0255728;
        const double f = 2.22858;
        const double g = 0.000449543;
        const double h = 2.22569;
        return (a * z * z + c * z * z * z * z) * (1 + f * z * z) /
               (1 + b * z * z + g * z * z * z * z) *
               std::exp(h - std::sqrt(4 * z * z + h * h));
    } else if (source == "B") {
        // Eq.(24) of A. Sedrakian, Phys. Rev. D93, 065044 (2016),
        // arXiv:1512.07828
        const double a = (2. / 3.) * 0.158151;
        const double b = -0.043745;
        const double c = -0.000271463;
        const double f = 0.0063470221;
        const double g = 0.0216661;
        return (a * z * z + c * z * z * z * z) * std::sqrt(1 + f * z * z) /
               (1 + b * z * z + g * z * z * z * z);
    }
    return 0;
}

double Pbf_3p2::GetDeltaT(double T, double Tc) { return GetDeltaT(T, Tc, 0); }

double Pbf_3p2::GetDeltaT(double T, double Tc, double cos_theta) {
    double tau = T / Tc;
    if (tau >= 1) {
        return 0;
    }
    // Eq.(10) and (11) of D. G. Yakolev et al, Phys.Usp. 42 (1999) 737-778,
    // astro-ph/9906456
    if (source == "A") {
        return T * std::sqrt(1. - tau) * (0.7893 + 1.188 / tau) *
               std::sqrt(1. + 3. * cos_theta * cos_theta);
    } else if (source == "B") {
        const double tau4 = tau * tau * tau * tau;
        return T * std::sqrt(1. - tau4) / tau *
               (2.030 - 0.4903 * tau4 + 0.1727 * tau4 * tau4) *
               std::sqrt(1 - cos_theta * cos_theta);
    }
    return 0;
}

void Pbf_3p2::GetBoundaries(double *rmin, double *rmax) {
    (*rmin) = 0;
    (*rmax) = nscool->Get1s03p2Boundary();
}
