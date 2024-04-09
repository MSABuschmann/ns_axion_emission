#ifndef NSCOOL_H_
#define NSCOOL_H_

#include <iostream>
#include <vector>
#include "hdf5.h"

class Process;

class NSCool {
  public:
    NSCool();

    bool LoadEos(std::string eos);
    void PrintEosNames();
    void WriteRawData(hid_t file_id);

    double GetRMax() const {
        return std::min(raw_rTc.back(), raw_rT.back()) - 1e-10;
    }

    double Get1s03p2Boundary() const {
        return boundary_1s03p2;
    }

    double GetT(double r) {
        return lerp(raw_rT, raw_T, r) * alpha;
    }

    double GetTcn(double r) {
        return lerp(raw_rTc, raw_Tcn, r);
    }

    double GetTcp(double r) {
        return lerp(raw_rTc, raw_Tcp, r);
    }

    double GetKfn(double r) {
        return lerp(raw_rTc, raw_kfn, r);
    }

    double GetKfp(double r) {
        return lerp(raw_rTc, raw_kfp, r);
    }

    double GetMstn(double r) {
        return lerp(raw_rTc, raw_mstn, r);
    }

    double GetMstp(double r) {
        return lerp(raw_rTc, raw_mstp, r);
    }

    double GetEphi(double r) {
        return lerp(raw_rT, raw_ephi, r);
    }

    double GetDvdr(double r) {
        return lerp(raw_rT, raw_dvdr, r);
    }

    void SetAlpha(double new_alpha) {
        alpha = new_alpha;
    }

    double GetUnweightedMeanT();
    void DetermineDeltaTInfty(Process* process);
    std::vector<double> GetResonanceLayer(double omega);

  private:
    void LocateEos();
    bool LoadCriticalTemps(std::string path);
    bool LoadProfile(std::string path);
    std::vector<std::string> Split(std::string line);
    double lerp(std::vector<double>& x, std::vector<double>& y, double nx);
    void Determine1s03p2Boundary();

    const std::string eos_parent_dir = "./eos/";
    std::vector<std::string> eos_names;
    double boundary_1s03p2 = 0;
    double alpha = 1;

    std::vector<double> raw_rTc;
    std::vector<double> raw_Tcn;
    std::vector<double> raw_Tcp;
    std::vector<double> raw_kfp;
    std::vector<double> raw_kfn;
    std::vector<double> raw_mstp;
    std::vector<double> raw_mstn;
    std::vector<int> raw_state;

    std::vector<double> raw_rT;
    std::vector<double> raw_T;
    std::vector<double> raw_ephi;
    std::vector<double> raw_dvol;
    std::vector<double> raw_dvdr;

    std::vector<double> raw_DeltaT_infty;
};

#endif // NSCOOL_H_
