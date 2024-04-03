#ifndef UTILS_H_
#define UTILS_H_

#include "hdf5.h"

#include "nscool.h"

void WriteDataset(hid_t file_id, const char* name,
                  std::vector<double> dataset) {
    hsize_t dims[1] = {dataset.size()};
    hid_t space = H5Screate_simple(1, dims, NULL);
    hid_t dataset_id = H5Dcreate(file_id, name, H5T_IEEE_F64LE, space,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             &dataset[0]);
    H5Dclose(dataset_id);
}

std::vector<double> CreateVector(double min, double max, int N) {
    std::vector<double> v(N);
    for (int i = 0; i < N; ++i) {
        v[i] = static_cast<double>(i)/(static_cast<double>(N-1)) * (max-min)
                + min;
    }
    return v;
}


void TestInterpolation(hid_t file_id, NSCool& nscool) {
    const int N = 1000;
    std::vector<double> r = CreateVector(0, nscool.GetRMax(), N); 
    std::vector<double> T(N);
    std::vector<double> Tcn(N);
    for (int i = 0; i < N; ++i) {
        T[i] = nscool.GetT(r[i]);
        Tcn[i] = nscool.GetTcn(r[i]);
    }
    WriteDataset(file_id, "r_interp", r);
    WriteDataset(file_id, "T_interp", T);
    WriteDataset(file_id, "Tcn_interp", Tcn);
}

#endif // UTILS_H_
