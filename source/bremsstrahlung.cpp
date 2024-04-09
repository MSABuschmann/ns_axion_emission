#include <cmath>

#include "bremsstrahlung.h"
#include "utils.h"

double Bremsstrahlung::Integrand(double r, double E) {
    double T = nscool->GetT(r);
    double ephi = nscool->GetEphi(r);
    double dvdr = nscool->GetDvdr(r);
    double omega = E * keV2K / ephi;
    double z = omega/T;
    if (z <= 0) {
        return 0;
    }
    return z*z*z*(z*z + 4.*M_PI)/(std::exp(z)-1.) * dvdr * ephi;
}

void Bremsstrahlung::GetBoundaries(double* rmin, double* rmax) {
    (*rmin) = 0;
    (*rmax) = nscool->GetRMax();
}

