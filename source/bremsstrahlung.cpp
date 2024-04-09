#include <cmath>

#include "bremsstrahlung.h"
#include "utils.h"

double Bremsstrahlung::Integrand(double r, double E) {
    const double T = nscool->GetT(r);
    const double ephi = nscool->GetEphi(r);
    const double omega = E * keV2K / ephi;
    const double z = omega/T;
    if (z <= 0) {
        return 0;
    }

    const double dvdr = nscool->GetDvdr(r);
    const double kfn = nscool->GetKfn(r);
    const double kfp = nscool->GetKfp(r);
    const double mstn = nscool->GetMstn(r);
    const double beta_nn = 0.56;
    const double beta_np = 0.66;
    const double beta_pp = 0.7;
    const double gamma_nn = 0.838;
    const double gamma_np = 0.838;
    const double gamma_pp = 0.838;
    const double gamma = 1/(1. + 1./3. * mstn * (kfn/1.68) );

    const double T8 = T/1e8;
    const double gann10 = gann*1e10;
    const double gapp10 = gapp*1e10;

    const double xn = 0.207 * 1.68 / kfn;
    const double xp = 0.207 * 1.68 / kfp;
    const double Fn = F(xn);
    const double Fp = F(xp);
    const double Gp = G(xp);
    const double Fplus = F(2.*xn*xp / (xp+xn));
    const double Fminus = F(2.*xn*xp / (xp-xn));

    const double Cg = 0.5 *Fp + Fplus + Fminus + xp/xn*(Fplus - Fminus)  + Gp;
    const double Ch = 0.5*(Fp + Fplus + Fminus + xp/xn*(Fplus - Fminus)) + Gp;
    const double geff10 = std::sqrt((gapp10+gann10)*(gapp10+gann10)*Cg
                                  + (gapp10-gann10)*(gapp10-gann10)*Ch);

    const double Rnn = 1;
    const double Rnp = 1;
    const double Rpp = 1;

    double epsilon = 0;
    if (source == "n") {
        epsilon = 7.373e11 * gann10*gann10 * (Fn/0.601566) * (kfn/1.68)
                * T8*T8*T8*T8*T8*T8 * (beta_nn/0.56) * (gamma_nn/0.838)
                * gamma*gamma*gamma*gamma*gamma*gamma * Rnn;
    } else if (source == "np") {
        epsilon = 9.617e11 * geff10*geff10 * (kfp/1.68)
                * T8*T8*T8*T8*T8*T8 * (beta_np/0.56) * (gamma_np/0.838)
                * gamma*gamma*gamma*gamma*gamma*gamma * Rnp;
    } else if (source == "p") {
        epsilon = 9.191e11 * gapp10*gapp10 * (Fp/0.601566) * (kfp/1.68)
                * T8*T8*T8*T8*T8*T8 * (beta_pp/0.56) * (gamma_pp/0.838)
                * gamma*gamma*gamma*gamma*gamma*gamma * Rpp;
    }

    const double N = 315./(124 * M_PI*M_PI*M_PI*M_PI*M_PI*M_PI);
    const double J = N * z*z*z*(z*z + 4.*M_PI)/(std::exp(z)-1.)/T * dvdr * ephi;
    return J * epsilon;
}

double Bremsstrahlung::F(double x) {
    return 1. - 1.5 * x * std::atan(1./x) + 0.5 * x*x / (1 + x*x);
}

double Bremsstrahlung::G(double x) {
    return 1. - x * std::atan(1./x);
}

void Bremsstrahlung::GetBoundaries(double* rmin, double* rmax) {
    if (source == "n") {
        (*rmin) = 0;
        (*rmax) = nscool->GetRMax();
    } else {
        (*rmin) = 0;
        (*rmax) = nscool->Get1s03p2Boundary();
    }
}
