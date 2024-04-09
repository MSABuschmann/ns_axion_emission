#ifndef BREMSSTRAHLUNG_H_
#define BREMSSTRAHLUNG_H_

#include "process.h"

class Bremsstrahlung: public Process {
  public:
    Bremsstrahlung(NSCool* pnscool, std::string nucleon)
        : Process(pnscool), source{nucleon} {};
    using Process::Process;
    virtual double GetDeltaT(double T, double Tc) override { return 0; };
    virtual double Integrand(double r, double E) override;

  protected:
    virtual void GetBoundaries(double* rmin, double* rmax) override;

  private:
    std::string source;
};

#endif // BREMSSTRAHLUNG_H_
