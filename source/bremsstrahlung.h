#ifndef BREMSSTRAHLUNG_H_
#define BREMSSTRAHLUNG_H_

#include "process.h"

class Bremsstrahlung: public Process {
  public:
    Bremsstrahlung(NSCool* pnscool, std::string nucleon, double this_gann,
                   double this_gapp)
        : Process(pnscool), source{nucleon}, gann{this_gann},
          gapp{this_gapp} {};
    using Process::Process;
    virtual double GetDeltaT(double T, double Tc) override { return 0; };
    virtual double Integrand(double r, double E) override;

  protected:
    virtual void GetBoundaries(double* rmin, double* rmax) override;

  private:
    double F(double x);
    double G(double x);
    std::string source;
    double gann, gapp;
};

#endif // BREMSSTRAHLUNG_H_
