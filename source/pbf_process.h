#ifndef PBF_PROCESS_H_
#define PBF_PROCESS_H_

#include "nscool.h"

class PbfProcess {
  public:
    explicit PbfProcess(NSCool* pnscool, State this_state)
        : nscool{pnscool}, state{this_state} {};
    std::vector<double> GetSpectrum(std::vector<double>& E_bins, int N);

  private:
    virtual double J(double omega, double T, double DeltaT) = 0;
    std::vector<double> GetDeltaT(std::vector<double> T,
                                  std::vector<double> Tc);
    virtual double GetDeltaT(double T, double Tc) = 0;
    double Integrand(double E, double r);
    void Normalize(std::vector<double>& x, std::vector<double>& y);

    NSCool* nscool;
    const State state;
};

class PbfProcess_1s0: public PbfProcess {
  public:
    PbfProcess_1s0(NSCool* pnscool) : PbfProcess(pnscool, SF_1s0) {};
    using PbfProcess::PbfProcess;

  private:
    virtual double J(double omega, double T, double DeltaT) override;
    virtual double GetDeltaT(double T, double Tc) override;
};

#endif // PBF_PROCESS_H_
