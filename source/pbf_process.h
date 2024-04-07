#ifndef PBF_PROCESS_H_
#define PBF_PROCESS_H_

#include "nscool.h"

#define keV2K 11605000.

class PbfProcess {
  public:
    explicit PbfProcess(NSCool* pnscool, State this_state)
        : nscool{pnscool}, state{this_state} {};
    std::vector<double> GetSpectrum(std::vector<double>& E_bins, int N);
    std::vector<double> GetSpectrum(std::vector<double>& E_bins);
    virtual double GetDeltaT(double T, double Tc) = 0;
    virtual double Integrand(double r, double E, PbfProcess* pthis) = 0;
    virtual double J(double omega, double T, double DeltaT) = 0;
    NSCool* nscool;

  private:
    static double GslIntegrand(double r, void* params);
    void GetBoundaries(double* rmin, double* rmax);
    void Normalize(std::vector<double>& x, std::vector<double>& y);

    const State state;
};

struct GslIntegrationParams {
    PbfProcess* pthis;
    double E;
};

#endif // PBF_PROCESS_H_
