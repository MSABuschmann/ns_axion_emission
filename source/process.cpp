#include <chrono>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <numeric>

#include "process.h"
#include "utils.h"

#define QAGS 1
#define QAGP 2
#define CQUAD 3
#define ROMBERG 4
#define SCHEME ROMBERG

// Rescaling can help numerical stability
#define FUDGE 1e26

double Process::GslIntegrand(double r, void *params) {
    GslIntegrationParams *gsl_params =
        static_cast<GslIntegrationParams *>(params);
    return gsl_params->pthis->Integrand(r, gsl_params->E) / FUDGE;
}

std::vector<double> Process::GetSpectrum(std::vector<double> &E_bins) {
    std::cout << "Compute spectrum with GSL ..." << std::endl;
    std::chrono::steady_clock::time_point start_time =
        std::chrono::steady_clock::now();
    double rmin = 0, rmax = 0;
    GetBoundaries(&rmin, &rmax);
    std::vector<double> spectrum(E_bins.size(), 0.);
#if SCHEME == QAPG
    nscool->DetermineDeltaTInfty(this);
#endif

#pragma omp parallel for
    for (size_t i = 0; i < E_bins.size(); ++i) {
        GslIntegrationParams gsl_params;
        gsl_params.pthis = this;
        gsl_params.E = E_bins[i];
        gsl_function F;
        F.function = &Process::GslIntegrand;
        F.params = &gsl_params;

#if SCHEME == QAGS
        double error;
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
        gsl_integration_qags(&F, rmin, rmax, 0.1, 1e-2, 1000, w, &spectrum[i],
                             &error);
#elif SCHEME == QAPG
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
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
        gsl_integration_qagp(&F, &pts[0], pts.size(), 1e-2, 1e-3, 1000, w,
                             &spectrum[i], &error);
#elif SCHEME == CQUAD
        double error;
        size_t neval;
        gsl_integration_cquad_workspace *w =
            gsl_integration_cquad_workspace_alloc(10000);
        gsl_integration_cquad(&F, rmin, rmax, 0, 1e-6, w, &spectrum[i], &error,
                              &neval);
#elif SCHEME == ROMBERG
        size_t neval;
        gsl_integration_romberg_workspace *w =
            gsl_integration_romberg_alloc(24);
        gsl_integration_romberg(&F, rmin, rmax, 0, 1e-3, &spectrum[i], &neval,
                                w);
#endif
    }
    PrintDuration(start_time);
    for (size_t i = 0; i < spectrum.size(); ++i) {
        spectrum[i] *= FUDGE;
    }
    return spectrum;
}

std::vector<double> Process::GetSpectrum(std::vector<double> &E_bins, int N) {
    std::cout << "Compute spectrum with N_r = " << N << " ..." << std::endl;
    std::chrono::steady_clock::time_point start_time =
        std::chrono::steady_clock::now();
    double rmin = 0, rmax = 0;
    GetBoundaries(&rmin, &rmax);
    std::vector<double> r = CreateVector(rmin, rmax, N);
    std::vector<double> spectrum(E_bins.size(), 0.);
#pragma omp parallel for
    for (size_t i = 0; i < E_bins.size(); ++i) {
        for (int j = 0; j < N; ++j) {
            spectrum[i] += Integrand(r[j], E_bins[i]);
        }
    }
    PrintDuration(start_time);
    return spectrum;
}

void Process::Normalize(std::vector<double> &x, std::vector<double> &y) {
    double sum = std::accumulate(y.begin(), y.end(), 0.);
    double norm = static_cast<double>(x.size()) / (sum * (x.back() - x[0]));
    for (size_t i = 0; i < y.size(); ++i) {
        y[i] *= norm;
    }
}
