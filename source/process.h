#ifndef PROCESS_H_
#define PROCESS_H_

#include <cmath>

#include "nscool.h"

class Process {
  public:
    explicit Process(NSCool *pnscool, std::string this_source)
        : nscool{pnscool}, source{this_source} {};
    std::vector<double> Get_dI_dE_Gsl(std::vector<double> &E_bins, int N);
    std::vector<double> Get_dI_dE(std::vector<double> &E_bins, int N);
    virtual double GetDeltaT(double T, double Tc) = 0;
    virtual double dI_dE_dr(double r, double E) = 0;

  protected:
    double Fermi(double x) { return 1. / (std::exp(x) + 1); };
    NSCool *nscool;
    std::string source;

  private:
    static double dI_dE_dr_Gsl(double r, void *params);
    void Normalize(std::vector<double> &x, std::vector<double> &y);
};

struct GslIntegrationParams {
    Process *pthis;
    double E;
};

#endif // PROCESS_H_
