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

void PrintDuration(std::chrono::steady_clock::time_point& start_time) {
    std::chrono::steady_clock::time_point stop_time =
            std::chrono::steady_clock::now();
    std::chrono::microseconds duration_micro =
            std::chrono::duration_cast<std::chrono::microseconds>(
                    stop_time - start_time);
    double duration_seconds = static_cast<double>(duration_micro.count())/1e6;
    std::cout << "Time: " << duration_seconds << "s." << std::endl;
}
