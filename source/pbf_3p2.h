#ifndef PBF_3P2_H_
#define PBF_3P2_H_

#include "process.h"

class Pbf_3p2 : public Process {
  public:
    Pbf_3p2(NSCool *pnscool, std::string this_source, double this_gann)
        : Process(pnscool, this_source), gann{this_gann} {};
    using Process::Process;
    virtual double dI_dE_dr(double r, double E) override;
    virtual double GetDeltaT(double T, double Tc) override;
    double GetDeltaT(double T, double Tc, double cos_theta);

  private:
    static double dI_dE_dr_dtheta(double r, void *params);
    double J(double omega, double T, double DeltaT, double I0);
    double Epsilon(double r, double T, double I);
    double Ian(double T, double DeltaT);

    double gann;
};

struct GslThetaIntegrationParams {
    Pbf_3p2 *pthis;
    double r;
    double omega;
    double T;
    double Tc;
    double I0;
};

#endif // PBF_3P2__H_
