#include <cmath>
#include <numeric>
#include <gsl/gsl_integration.h>
#include <chrono>

#include "pbf_process.h"
#include "utils.h"

#define QAGS 1
#define QAGP 2
#define CQUAD 3
#define ROMBERG 4
#define SCHEME ROMBERG

double PbfProcess::GslIntegrand(double r, void* params) {
    GslIntegrationParams* gsl_params =
            static_cast<GslIntegrationParams*>(params);
    return gsl_params->pthis->Integrand(r, gsl_params->E, gsl_params->pthis);
}

std::vector<double> PbfProcess::GetSpectrum(std::vector<double>& E_bins) {
    std::cout << "Compute spectrum with GSL ..." << std::endl;
    std::chrono::steady_clock::time_point start_time =
            std::chrono::steady_clock::now();
    double rmin=0, rmax=0;
    GetBoundaries(&rmin, &rmax);
    std::vector<double> spectrum(E_bins.size(), 0.);
#if SCHEME==QAPG
    nscool->DetermineDeltaTInfty(this);
#endif

#pragma omp parallel for
    for (size_t i = 0; i < E_bins.size(); ++i) {
        GslIntegrationParams gsl_params;
        gsl_params.pthis = this;
        gsl_params.E = E_bins[i];
        gsl_function F;
        F.function = &PbfProcess::GslIntegrand;
        F.params = &gsl_params;

#if SCHEME==QAGS
        double error;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
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
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
        gsl_integration_qagp(&F, &pts[0], pts.size(), 1e-2, 1e-3, 1000, w,
                             &spectrum[i], &error);
#elif SCHEME==CQUAD
        double error;
        size_t neval;
        gsl_integration_cquad_workspace* w =
                gsl_integration_cquad_workspace_alloc(10000);
        gsl_integration_cquad(&F, rmin, rmax, 0, 1e-6, w, &spectrum[i], &error,
                              &neval);
#elif SCHEME==ROMBERG
        size_t neval;
        gsl_integration_romberg_workspace* w =
                gsl_integration_romberg_alloc(19);
        gsl_integration_romberg(&F, rmin, rmax, 0, 1e-3, &spectrum[i], &neval,
                                w);
#endif
    }
    Normalize(E_bins, spectrum);
    PrintDuration(start_time);
    return spectrum;
}

std::vector<double> PbfProcess::GetSpectrum(std::vector<double>& E_bins,
                                            int N) {
    std::cout << "Compute spectrum with N_r = " << N << " ..." << std::endl;
    std::chrono::steady_clock::time_point start_time =
            std::chrono::steady_clock::now();
    double rmin=0, rmax=0;
    GetBoundaries(&rmin, &rmax);
    std::vector<double> r = CreateVector(rmin, rmax, N);
    std::vector<double> spectrum(E_bins.size(), 0.);
#pragma omp parallel for
    for (size_t i = 0; i < E_bins.size(); ++i) {
        for (int j = 0; j < N; ++j) {
            spectrum[i] += Integrand(r[j], E_bins[i], this);
        }
    }
    Normalize(E_bins, spectrum);
    PrintDuration(start_time);
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
    for (size_t i = 0; i < y.size(); ++i) {
        y[i] *= norm;
    }
}
