#include <cmath>
#include <numeric>
#include <gsl/gsl_integration.h>
#include <chrono>

#include "pbf_process.h"
#include "utils.h"

#define keV2K 11605000.
#define QAGS 1
#define QAGP 2
#define CQUAD 3
#define ROMBERG 4
#define SCHEME ROMBERG

double PbfProcess::Integrand(double r, double E, PbfProcess* pthis) {
    double T = pthis->nscool->GetT(r);
    double Tc = pthis->nscool->GetTcn(r);
    double DeltaT = pthis->GetDeltaT(T, Tc);
    if (DeltaT <= 0) {
        return 0.;
    }
    double ephi = pthis->nscool->GetEphi(r);
    double dvdr = pthis->nscool->GetDvdr(r);
    double omega = E * keV2K / ephi;
    return pthis->J(omega, T, DeltaT) * dvdr * ephi;
}

double PbfProcess::GslIntegrand(double r, void* params) {
    GslIntegrationParams* gsl_params =
            static_cast<GslIntegrationParams*>(params);
    return Integrand(r, gsl_params->E, gsl_params->pthis);
}

std::vector<double> PbfProcess::GetSpectrum(std::vector<double>& E_bins) {
    std::cout << "Compute spectrum with GSL ..." << std::endl;
    double rmin=0, rmax=0;
    GetBoundaries(&rmin, &rmax);
    //nscool->DetermineDeltaTInfty(this);

#if SCHEME<CQUAD
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
#elif SCHEME==CQUAD
    gsl_integration_cquad_workspace* w = gsl_integration_cquad_workspace_alloc(10000);
#elif SCHEME==ROMBERG
    gsl_integration_romberg_workspace* w = gsl_integration_romberg_alloc(19);
#else
    std::cout << "No apropriate integrator selected during compilation!"
              << std::endl;
    return NULL;
#endif
    gsl_function F;
    F.function = &PbfProcess::GslIntegrand;

    std::vector<double> spectrum(E_bins.size(), 0.);
    for (unsigned int i = 0; i < E_bins.size(); ++i) {
        GslIntegrationParams gsl_params;
        gsl_params.pthis = this;
        gsl_params.E = E_bins[i];
        F.params = &gsl_params;

#if SCHEME==QAGS
        double error;
        gsl_integration_qags(&F, rmin, rmax, 0.1, 1e-2, 1000, w, &spectrum[i],
                             &error);
#elif SCHEME==QAPG
        std::vector<double> resonances =
                nscool->GetResonanceLayer(E_bins[i] * keV2K);
        std::vector<double> pts = {rmin};
        for (double res : resonances) {
            if (rmin < res && res < rmax) {
                pts.push_back(res);
            }
        }
        pts.push_back(rmax);
        double error;
        gsl_integration_qagp(&F, &pts[0], pts.size(), 1e-2, 1e-3, 1000, w,
                             &spectrum[i], &error);
#elif SCHEME==CQUAD
        double error;
        size_t neval;
        gsl_integration_cquad(&F, rmin, rmax, 0, 1e-6, w, &spectrum[i],
                              &error,  &neval);
#elif SCHEME==ROMBERG
        size_t neval;
        gsl_integration_romberg(&F, rmin, rmax, 0, 1e-3, &spectrum[i], &neval,
                                w);
#endif
    }
    Normalize(E_bins, spectrum);
    return spectrum;
}

std::vector<double> PbfProcess::GetSpectrum(std::vector<double>& E_bins,
                                            int N) {
    std::cout << "Compute spectrum with N_r = " << N << " ..." << std::endl;
    double rmin=0, rmax=0;
    GetBoundaries(&rmin, &rmax);
    std::vector<double> r = CreateVector(rmin, rmax, N);
    //nscool->DetermineDeltaTInfty(this);

    std::vector<double> spectrum(E_bins.size(), 0.);
    for (unsigned int i = 0; i < E_bins.size(); ++i) {
        for (int j = 0; j < N; ++j) {
            spectrum[i] += Integrand(r[j], E_bins[i], this);
        }
    }
    Normalize(E_bins, spectrum);
    return spectrum;
}

void PbfProcess::GetBoundaries(double* rmin, double* rmax) {
    if (state == SF_3p2) {
        (*rmin) = 0;
        (*rmax) = nscool->Get1s03p2Boundary();
    } else if (state == SF_1s0) {
        (*rmin) = nscool->Get1s03p2Boundary();
        (*rmax) = nscool->GetRMax();
    }
}

void PbfProcess::Normalize(std::vector<double>& x, std::vector<double>& y) {
    double sum = std::accumulate(y.begin(), y.end(), 0.);
    double norm = 1./(sum*(x.back()-x[0]));
    for (unsigned int i = 0; i < y.size(); ++i) {
        y[i] *= norm;
    }
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
