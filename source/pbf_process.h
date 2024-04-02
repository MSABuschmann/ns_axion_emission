#ifndef PBF_PROCESS_H_
#define PBF_PROCESS_H_

#include "nscool.h"

class PbfProcess {
  public:
    explicit PbfProcess(NSCool* pnscool) : nscool(pnscool) {};
    std::vector<double> GetSpectrum(std::vector<double>& E_bins);

  private:
    virtual double J(double omega, double T, double DeltaT) = 0;
    std::vector<double> GetDeltaT(std::vector<double> T,
                                  std::vector<double> Tc);
    virtual double GetDeltaT(double T, double Tc) = 0;
    NSCool* nscool;
};

class PbfProcess_1s0: public PbfProcess {
  public:
    using PbfProcess::PbfProcess;

  private:
    virtual double J(double omega, double T, double DeltaT) override;
    virtual double GetDeltaT(double T, double Tc) override;

    const State state = SF_1s0;
};

#endif // PBF_PROCESS_H_
