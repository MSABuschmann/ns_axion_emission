#ifndef BREMSSTRAHLUNG_H_
#define BREMSSTRAHLUNG_H_

#include "process.h"

class Bremsstrahlung : public Process {
  public:
    Bremsstrahlung(NSCool *pnscool, std::string this_source, double this_gann,
                   double this_gapp)
        : Process(pnscool, this_source), gann{this_gann}, gapp{this_gapp} {};
    using Process::Process;

    virtual double GetDeltaT(double T, double Tc) override {
        return T * DeltaT0(T, Tc);
    };

    virtual double dI_dE_dr(double r, double E) override;

    double DeltaT0(double T, double Tc);
    double Rnn(double T, double Tc);
    double Rnp(double T, double Tcn, double Tcp);
    double Rnp(double T, double Tc);
    double Rpp(double T, double Tc);

  private:
    double F(double x);
    double G(double x);
    double J(double omega, double T);
    double Epsilon(double r, double T);

    double gann, gapp;
};

#endif // BREMSSTRAHLUNG_H_
