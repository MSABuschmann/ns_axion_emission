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

    std::vector<double> E_bins = CreateVector(0., 150., 1000);
    PbfProcess_1s0 pbf_1s0 = PbfProcess_1s0(&nscool);

    std::string filename = "output/" + eos + ".hdf5";
    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
    WriteDataset(file_id, "E", E_bins);
    WriteDataset(file_id, "1s0n_N100", pbf_1s0.GetSpectrum(E_bins, 100));
    WriteDataset(file_id, "1s0n_N1000", pbf_1s0.GetSpectrum(E_bins, 1000));
    WriteDataset(file_id, "1s0n_N10000", pbf_1s0.GetSpectrum(E_bins, 10000));
    WriteDataset(file_id, "1s0n_N100000", pbf_1s0.GetSpectrum(E_bins, 100000));
    TestInterpolation(file_id, nscool);
    nscool.WriteRawData(file_id);
    H5Fclose(file_id);

    return 0;
}
