#ifndef PBF_3P2_H_
#define PBF_3P2_H_

#include "process.h"

class Pbf_3p2 : public Process {
  public:
    Pbf_3p2(NSCool *pnscool, std::string nucleon, double this_gann)
        : Process(pnscool), source{nucleon}, gann{this_gann} {};
    using Process::Process;
    virtual double Integrand(double r, double E) override;
    virtual double GetDeltaT(double T, double Tc) override;
    double GetDeltaT(double T, double Tc, double theta);

  protected:
    virtual void GetBoundaries(double *rmin, double *rmax) override;

  private:
    double J(double omega, double T, double DeltaT, double I);
    double Epsilon(double r, double T, double I);
    double Ian(double T, double DeltaT);

    std::string source;
    double gann;
};

#endif // PBF_3P2__H_
