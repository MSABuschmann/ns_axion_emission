#ifndef NSCOOL_H_
#define NSCOOL_H_

#include "hdf5.h"
#include <algorithm>
#include <iostream>
#include <vector>

class Process;

class NSCool {
  public:
    NSCool();

    bool LoadEos(std::string eos);
    void PrintEosNames();
    void WriteRawData(hid_t file_id);

    double GetRMax() const { return r_max; }
    double GetCoreCrustBoundary() const { return boundary_core_crust; }
    bool GetBoundaries(double *rmin, double *rmax, std::string source);
    double GetT(double r) { return lerp(raw_rT, raw_T, r) * alpha; }
    double GetTcn(double r) { return lerp(raw_rTc, raw_Tcn, r); }
    double GetTcp(double r) { return lerp(raw_rTc, raw_Tcp, r); }
    double GetKfn(double r) { return lerp(raw_rTc, raw_kfn, r); }
    double GetKfp(double r) { return lerp(raw_rTc, raw_kfp, r); }
    double GetMstn(double r) { return lerp(raw_rTc, raw_mstn, r); }
    double GetMstp(double r) { return lerp(raw_rTc, raw_mstp, r); }
    double GetEphi(double r) { return lerp(raw_rT, raw_ephi, r); }
    double GetDvdr(double r) { return lerp(raw_rT, raw_dvdr, r); }
    void SetAlpha(double new_alpha) { alpha = new_alpha; }

    double GetUnweightedMeanT();
    void DetermineDeltaTInfty(Process *process);
    std::vector<double> GetResonanceLayer(double omega);

    template <typename T> double FindDownwards(std::vector<T> &vec, T val);
    template <typename T> double FindUpwards(std::vector<T> &vec, T val);

  private:
    void LocateEos();
    bool LoadCriticalTemps(std::string path);
    bool LoadProfile(std::string path);
    std::vector<std::string> Split(std::string line);
    double lerp(std::vector<double> &x, std::vector<double> &y, double nx);
    int GetState(std::string state);
    void DetermineBoundaries();
    void PrintBoundary(std::string process, std::string source);
    void PrintBoundaries();

    const std::string eos_parent_dir = "./eos/";
    std::vector<std::string> eos_names;
    double alpha = 1;

    double boundary_3p2 = -1;
    double boundary_1s0n_lower = -1;
    double boundary_1s0n_upper = -1;
    double boundary_1s0p_lower = -1;
    double boundary_1s0p_upper = -1;
    double boundary_core_crust = -1;
    double r_max = -1;

    std::vector<double> raw_rTc;
    std::vector<double> raw_Tcn;
    std::vector<double> raw_Tcp;
    std::vector<double> raw_kfp;
    std::vector<double> raw_kfn;
    std::vector<double> raw_mstp;
    std::vector<double> raw_mstn;
    std::vector<double> raw_state;

    std::vector<double> raw_rT;
    std::vector<double> raw_T;
    std::vector<double> raw_ephi;
    std::vector<double> raw_dvol;
    std::vector<double> raw_dvdr;

    std::vector<double> raw_DeltaT_infty;
};

template <typename T> double NSCool::FindDownwards(std::vector<T> &vec, T val) {
    int loc = -1;
    for (size_t i = 1; i < vec.size(); ++i) {
        if (vec[i] < val && vec[i - 1] >= val) {
            loc = i;
            break;
        }
    }

    if (loc <= 0) {
        if (vec.back() > val) {
            return raw_rTc.back();
        }

        return -1;
    }

    return (raw_rTc[loc] + raw_rTc[loc - 1]) / 2.;
}

template <typename T> double NSCool::FindUpwards(std::vector<T> &vec, T val) {
    int loc = -1;
    for (size_t i = 1; i < vec.size(); ++i) {
        if (vec[i] > val && vec[i - 1] <= val) {
            loc = i;
            break;
        }
    }

    if (loc <= 0) {
        return -1;
    }

    return (raw_rTc[loc] + raw_rTc[loc - 1]) / 2.;
}

#endif // NSCOOL_H_
