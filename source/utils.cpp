#include "utils.h"
#include "bremsstrahlung.h"
#include "pbf_1s0.h"

std::vector<double> CreateVector(double min, double max, int N) {
    std::vector<double> v(N);
    for (int i = 0; i < N; ++i) {
        v[i] = static_cast<double>(i) / (static_cast<double>(N - 1)) *
                   (max - min) +
               min;
    }
    return v;
}

void TestInterpolation(hid_t file_id, NSCool &nscool) {
    const int N = 1000;
    std::vector<double> r = CreateVector(0, nscool.GetRMax(), N);
    std::vector<double> T(N);
    std::vector<double> Tcn(N);
    std::vector<double> Tcp(N);
    for (int i = 0; i < N; ++i) {
        T[i] = nscool.GetT(r[i]);
        Tcn[i] = nscool.GetTcn(r[i]);
        Tcp[i] = nscool.GetTcp(r[i]);
    }
    WriteDataset(file_id, "r_interp", r);
    WriteDataset(file_id, "T_interp", T);
    WriteDataset(file_id, "Tcn_interp", Tcn);
    WriteDataset(file_id, "Tcp_interp", Tcp);

    Bremsstrahlung bremsstrahlung(&nscool, "nn", 1e-10, 1e-10);
    std::vector<double> Rnn(N);
    std::vector<double> Rpp(N);
    std::vector<double> Rnp(N);
    for (int i = 0; i < N; ++i) {
        Rnn[i] = bremsstrahlung.Rnn(T[i], Tcn[i]);
        Rpp[i] = bremsstrahlung.Rpp(T[i], Tcp[i]);
        Rnp[i] = bremsstrahlung.Rnp(T[i], Tcn[i], Tcp[i]);
    }
    WriteDataset(file_id, "Rnn", Rnn);
    WriteDataset(file_id, "Rpp", Rpp);
    WriteDataset(file_id, "Rnp", Rnp);
}

void PrintDuration(std::chrono::steady_clock::time_point &start_time) {
    std::chrono::steady_clock::time_point stop_time =
        std::chrono::steady_clock::now();
    std::chrono::microseconds duration_micro =
        std::chrono::duration_cast<std::chrono::microseconds>(stop_time -
                                                              start_time);
    double duration_seconds = static_cast<double>(duration_micro.count()) / 1e6;
    std::cout << "Time: " << duration_seconds << "s." << std::endl;
}

void WriteIntegrand(hid_t file_id, NSCool &nscool) {
    const int N = 500000;
    std::vector<double> E = {5, 10, 15};
    std::vector<double> r = CreateVector(0, nscool.GetRMax(), N);
    WriteDataset(file_id, "integrand_r", r);
    WriteDataset(file_id, "integrand_E", E);

    Pbf_1s0 pbf_1s0n(&nscool, "n", 1e-10);
    for (size_t i = 0; i < E.size(); ++i) {
        std::vector<double> integrand(N, 0), int_p(N, 0);
        for (size_t j = 0; j < N; ++j) {
            integrand[j] = pbf_1s0n.Integrand(r[j], E[i]);
        }
        std::string str = "integrand_n_" + std::to_string(i);
        WriteDataset(file_id, str.c_str(), integrand);
    }

    Pbf_1s0 pbf_1s0p(&nscool, "p", 1e-10);
    for (size_t i = 0; i < E.size(); ++i) {
        std::vector<double> integrand(N, 0), int_p(N, 0);
        for (size_t j = 0; j < N; ++j) {
            integrand[j] = pbf_1s0p.Integrand(r[j], E[i]);
        }
        std::string str = "integrand_p_" + std::to_string(i);
        WriteDataset(file_id, str.c_str(), integrand);
    }
}
