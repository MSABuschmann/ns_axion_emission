#include "utils.h"

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
