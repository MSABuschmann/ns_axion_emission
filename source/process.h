#ifndef PROCESS_H_
#define PROCESS_H_

#include <cmath>

#include "nscool.h"

class Process {
  public:
    explicit Process(NSCool *pnscool) : nscool{pnscool} {};
    std::vector<double> GetSpectrumGsl(std::vector<double> &E_bins, int N);
    std::vector<double> GetSpectrum(std::vector<double> &E_bins, int N);
    virtual double GetDeltaT(double T, double Tc) = 0;
    virtual double Integrand(double r, double E) = 0;

  protected:
    virtual void GetBoundaries(double *rmin, double *rmax) = 0;
    double Fermi(double x) { return 1. / (std::exp(x) + 1); };
    NSCool *nscool;

  private:
    static double GslIntegrand(double r, void *params);
    void Normalize(std::vector<double> &x, std::vector<double> &y);
};

struct GslIntegrationParams {
    Process *pthis;
    double E;
};

#endif // PROCESS_H_
