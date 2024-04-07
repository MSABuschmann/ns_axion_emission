#ifndef PBF_1S0_H_
#define PBF_1S0_H_

#include "pbf_process.h"

class PbfProcess_1s0: public PbfProcess {
  public:
    PbfProcess_1s0(NSCool* pnscool) : PbfProcess(pnscool, SF_1s0) {};
    using PbfProcess::PbfProcess;
    virtual double GetDeltaT(double T, double Tc) override;
    virtual double Integrand(double r, double E, PbfProcess* pthis) override;

  private:
    virtual double J(double omega, double T, double DeltaT) override;
};

#endif // PBF_1S0_H_
