#ifndef PBF_1S0_H_
#define PBF_1S0_H_

#include "process.h"

class Pbf_1s0: public Process {
  public:
    Pbf_1s0(NSCool* pnscool, std::string nucleon)
        : Process(pnscool), source{nucleon} {};
    using Process::Process;
    virtual double GetDeltaT(double T, double Tc) override;
    virtual double Integrand(double r, double E) override;

  protected:
    virtual void GetBoundaries(double* rmin, double* rmax) override;

  private:
    std::string source;
    double J(double omega, double T, double DeltaT);
};

#endif // PBF_1S0_H_
