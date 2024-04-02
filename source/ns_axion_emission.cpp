//#include <sstream>
//#include <iostream>
//#include <fstream>
//#include <stdlib.h>
//#include "math.h"
//#include <omp.h>
//#include "hdf5.h"
#include <memory>

#include "nscool.h"
#include "pbf_process.h"

int main(int argc, const char* argv[]) {
    if (argc != 2) {
        std::cout << "No EOS provided!" << std::endl;
        return true;
    }
    std::string eos = argv[1];

    NSCool nscool;
    if (!nscool.LoadEos(eos)) {
        std::cout << "Could not load EOS '" << eos << "'!\n";
        std::cout << "Available EOSs are:\n";
        nscool.PrintEosNames();
    }

    const int NBins = 100;
    const double E_min = 0.;
    const double E_max = 150.;
    std::vector<double> E_bins(NBins);
    for (int i = 0; i < NBins; ++i) {
        E_bins[i] = i/(static_cast<double>(NBins-1.)) * (E_max-E_min) + E_min;
    }

    PbfProcess_1s0 pbf_1s0 = PbfProcess_1s0(&nscool);
    std::vector<double> spectrum_1s0 = pbf_1s0.GetSpectrum(E_bins);
    for (int i = 0; i < NBins; ++i) {
        std::cout << i << ": " << E_bins[i] << " " << spectrum_1s0[i]
                  << std::endl;
    }

    return 0;
}
