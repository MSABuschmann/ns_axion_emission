#include <cmath>

#include "bremsstrahlung.h"
#include "utils.h"

double Bremsstrahlung::Integrand(double r, double E) {
    if (E <= 0) {
        return 0;
    }
    const double T = nscool->GetT(r);
    const double dvdr = nscool->GetDvdr(r);
    const double ephi = nscool->GetEphi(r);
    const double omega = E * keV2K / ephi;
    return J(r, omega, T) * Epsilon(r, T) * dvdr * ephi;
}

inline double Bremsstrahlung::J(double r, double omega, double T) {
    const double z = omega / T;
    const double N = 315. / (124 * M_PI * M_PI * M_PI * M_PI * M_PI * M_PI);
    return N * z * z * z * (z * z + 4. * M_PI) / (std::exp(z) - 1.);
}

inline double Bremsstrahlung::Epsilon(double r, double T) {
    const double Tcn = nscool->GetTcn(r);
    const double Tcp = nscool->GetTcp(r);
    const double kfn = nscool->GetKfn(r);
    const double kfp = nscool->GetKfp(r);
    const double mstn = nscool->GetMstn(r);
    const double beta_nn = 0.56;
    const double beta_np = 0.66;
    const double beta_pp = 0.7;
    const double gamma_nn = 0.838;
    const double gamma_np = 0.838;
    const double gamma_pp = 0.838;
    const double gamma = 1 / (1. + mstn * (kfn / 1.68) / 3.);
    const double T8 = T / 1e8;
    const double gann10 = gann * 1e10;
    const double gapp10 = gapp * 1e10;
    const double xn = 0.207 * 1.68 / kfn;
    const double xp = 0.207 * 1.68 / kfp;
    const double Fn = F(xn);
    const double Fp = F(xp);
    const double Gp = G(xp);
    const double Fplus = F(2. * xn * xp / (xp + xn));
    const double Fminus = F(2. * xn * xp / (xp - xn));
    const double Cg =
        0.5 * Fp + Fplus + Fminus + xp / xn * (Fplus - Fminus) + Gp;
    const double Ch =
        0.5 * (Fp + Fplus + Fminus + xp / xn * (Fplus - Fminus)) + Gp;
    const double geff10 = std::sqrt((gapp10 + gann10) * (gapp10 + gann10) * Cg +
                                    (gapp10 - gann10) * (gapp10 - gann10) * Ch);
    const double rnn = Rnn(T, Tcn);
    const double rnp = Rnp(T, Tcn, Tcp);
    const double rpp = Rpp(T, Tcn);

    double epsilon = 0;
    if (source == "nn") {
        epsilon = 7.373e11 * gann10 * gann10 * (Fn / 0.601566) * (kfn / 1.68) *
                  T8 * T8 * T8 * T8 * T8 * T8 * (beta_nn / 0.56) *
                  (gamma_nn / 0.838) * gamma * gamma * gamma * gamma * gamma *
                  gamma * rnn;
    } else if (source == "np") {
        epsilon = 9.617e11 * geff10 * geff10 * (kfp / 1.68) * T8 * T8 * T8 *
                  T8 * T8 * T8 * (beta_np / 0.56) * (gamma_np / 0.838) * gamma *
                  gamma * gamma * gamma * gamma * gamma * rnp;
    } else if (source == "pp") {
        epsilon = 9.191e11 * gapp10 * gapp10 * (Fp / 0.601566) * (kfp / 1.68) *
                  T8 * T8 * T8 * T8 * T8 * T8 * (beta_pp / 0.56) *
                  (gamma_pp / 0.838) * gamma * gamma * gamma * gamma * gamma *
                  gamma * rpp;
    }
    return epsilon;
}

inline double Bremsstrahlung::F(double x) {
    return 1. - 1.5 * x * std::atan(1. / x) + 0.5 * x * x / (1 + x * x);
}

inline double Bremsstrahlung::G(double x) { return 1. - x * std::atan(1. / x); }

double Bremsstrahlung::Rnn(double T, double Tc) {
    const double x = T / Tc;
    if (x >= 1) {
        return 1;
    }
    const double v =
        std::sqrt(1. - x) * (1.456 - 0.157 / std::sqrt(x) + 1.764 / x);
    const double ann = 0.1747 + std::sqrt(0.68112009 + 0.0062932489 * v * v);
    const double bnn = 0.7333 + std::sqrt(0.07112889 + 0.02815684 * v * v);
    const double enn = 4.338 - std::sqrt(17.875984 + 4 * v * v);
    const double fnn = 7.762 - std::sqrt(60.248644 + 9 * v * v);
    return 0.5 *
           (ann * ann * std::exp(enn) + std::pow(bnn, 7.5) * std::exp(fnn));
}

double Bremsstrahlung::Rpp(double T, double Tc) { return Rnn(T, Tc); }

double Bremsstrahlung::Rnp(double T, double Tcn, double Tcp) {
    return std::min(Rnp(T, Tcn), Rnp(T, Tcp));
}

double Bremsstrahlung::Rnp(double T, double Tc) {
    const double x = T / Tc;
    if (x >= 1) {
        return 1;
    }
    const double v =
        std::sqrt(1. - x) * (1.456 - 0.157 / std::sqrt(x) + 1.764 / x);
    const double anp = 0.9982 + std::sqrt(0.00000324 + 0.14554225 * v * v);
    const double bnp = 0.3949 + std::sqrt(0.36614601 + 0.07107556 * v * v);
    const double enp = 1.306 - std::sqrt(1.705636 + v * v);
    const double fnp = 3.303 - std::sqrt(10.909809 + 4. * v * v);
    return 1. / 2.732 *
           (anp * std::exp(enp) +
            1.732 * bnp * bnp * bnp * bnp * bnp * bnp * bnp * std::exp(fnp));
}

inline void Bremsstrahlung::GetBoundaries(double *rmin, double *rmax) {
    if (source == "n") {
        (*rmin) = 0;
        (*rmax) = nscool->GetRMax();
    } else {
        (*rmin) = 0;
        (*rmax) = nscool->Get1s03p2Boundary();
    }
}
