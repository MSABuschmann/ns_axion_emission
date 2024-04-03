//#include <sstream>
//#include <iostream>
//#include <fstream>
//#include <stdlib.h>
//#include "math.h"
//#include <omp.h>
#include "hdf5.h"
#include <memory>

#include "nscool.h"
#include "pbf_process.h"
#include "utils.h"

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

    std::vector<double> E_bins = CreateVector(0., 150., 100);
    PbfProcess_1s0 pbf_1s0 = PbfProcess_1s0(&nscool);
    std::vector<double> spectrum_1s0 = pbf_1s0.GetSpectrum(E_bins);

    std::string filename = "output/" + eos + ".hdf5";
    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
    WriteDataset(file_id, "E", E_bins);
    WriteDataset(file_id, "1s0n", spectrum_1s0);

    TestInterpolation(file_id, nscool);
    WriteDataset(file_id, "raw_rTc", nscool.raw_rTc);
    WriteDataset(file_id, "raw_Tcn", nscool.raw_Tcn);
    WriteDataset(file_id, "raw_Tcp", nscool.raw_Tcp);
    WriteDataset(file_id, "raw_rT", nscool.raw_rT);
    WriteDataset(file_id, "raw_T", nscool.raw_T);
    WriteDataset(file_id, "raw_ephi", nscool.raw_ephi);
    H5Fclose(file_id);

    return 0;
}
