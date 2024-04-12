#include "hdf5.h"
#include <memory>

#include "bremsstrahlung.h"
#include "nscool.h"
#include "pbf_1s0.h"
#include "utils.h"

int main(int argc, const char *argv[]) {
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

    double gann = 1e-10;
    double gapp = 1e-10;
    std::vector<double> alpha = {0.5, 1., 1.5};
    std::vector<double> E = CreateVector(0., 80., 100);
    // std::vector<double> E = {5,10,15};

    Pbf_1s0 pbf_1s0n(&nscool, "n", gann);
    Pbf_1s0 pbf_1s0p(&nscool, "p", gapp);
    Bremsstrahlung bremsstrahlung_nn(&nscool, "nn", gann, gann);
    Bremsstrahlung bremsstrahlung_np(&nscool, "np", gann, gapp);
    Bremsstrahlung bremsstrahlung_pp(&nscool, "pp", gapp, gapp);

    // *************************************************** Write Header
    std::string filename = "output/" + eos + ".hdf5";
    hid_t file_id =
        H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    std::vector<double> params = {nscool.Get1s03p2Boundary(), nscool.GetRMax(),
                                  nscool.GetUnweightedMeanT()};
    WriteDataset(file_id, "E", E);
    WriteDataset(file_id, "alpha", alpha);
    WriteDataset(file_id, "params", params);
    nscool.WriteRawData(file_id);

    // *************************************************** Write Tests
    WriteIntegrand(file_id, nscool);
    TestInterpolation(file_id, nscool);
    WriteDataset(file_id, "1s0n_N100", pbf_1s0n.GetSpectrum(E, 100));
    WriteDataset(file_id, "1s0n_N1000", pbf_1s0n.GetSpectrum(E, 1000));
    WriteDataset(file_id, "1s0n_N10000", pbf_1s0n.GetSpectrum(E, 10000));
    WriteDataset(file_id, "1s0n_N262145", pbf_1s0n.GetSpectrum(E, 262145));

    // *************************************************** Write Spectra
    for (size_t i = 0; i < alpha.size(); ++i) {
        nscool.SetAlpha(alpha[i]);
        WriteDataset(file_id, "1s0n_a" + std::to_string(i),
                     pbf_1s0n.GetSpectrum(E));
        WriteDataset(file_id, "1s0p_a" + std::to_string(i),
                     pbf_1s0p.GetSpectrum(E));
        WriteDataset(file_id, "bremsstrahlung_nn_a" + std::to_string(i),
                     bremsstrahlung_nn.GetSpectrum(E));
        WriteDataset(file_id, "bremsstrahlung_np_a" + std::to_string(i),
                     bremsstrahlung_np.GetSpectrum(E));
        WriteDataset(file_id, "bremsstrahlung_pp_a" + std::to_string(i),
                     bremsstrahlung_pp.GetSpectrum(E));
    }

    H5Fclose(file_id);
    return 0;
}
