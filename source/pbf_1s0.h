#ifndef PBF_1S0_H_
#define PBF_1S0_H_

#include "process.h"

class Pbf_1s0 : public Process {
  public:
    Pbf_1s0(NSCool *pnscool, std::string nucleon, double this_gaNN)
        : Process(pnscool), source{nucleon}, gaNN{this_gaNN} {};
    using Process::Process;
    virtual double GetDeltaT(double T, double Tc) override;
    virtual double dI_dE_dr(double r, double E) override;

  protected:
    virtual void GetBoundaries(double *rmin, double *rmax) override;

  private:
    double J(double omega, double T, double DeltaT, double I);
    double Epsilon(double r, double T, double I);
    double Ias(double T, double DeltaT);

    std::string source;
    double gaNN;
};

#endif // PBF_1S0_H_
