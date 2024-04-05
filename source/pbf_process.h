#ifndef PBF_PROCESS_H_
#define PBF_PROCESS_H_

#include "nscool.h"

class PbfProcess {
  public:
    explicit PbfProcess(NSCool* pnscool, State this_state)
        : nscool{pnscool}, state{this_state} {};
    std::vector<double> GetSpectrum(std::vector<double>& E_bins, int N);
    std::vector<double> GetSpectrum(std::vector<double>& E_bins);
    virtual double GetDeltaT(double T, double Tc) = 0;

  private:
    virtual double J(double omega, double T, double DeltaT) = 0;
    static double Integrand(double r, double E, PbfProcess* pthis);
    static double GslIntegrand(double r, void* params);
    void GetBoundaries(double* rmin, double* rmax);
    void Normalize(std::vector<double>& x, std::vector<double>& y);

    NSCool* nscool;
    const State state;
};

class PbfProcess_1s0: public PbfProcess {
  public:
    PbfProcess_1s0(NSCool* pnscool) : PbfProcess(pnscool, SF_1s0) {};
    using PbfProcess::PbfProcess;
    virtual double GetDeltaT(double T, double Tc) override;

  private:
    virtual double J(double omega, double T, double DeltaT) override;
};

struct GslIntegrationParams {
    PbfProcess* pthis;
    double E;
};

#endif // PBF_PROCESS_H_
