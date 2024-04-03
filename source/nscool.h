#ifndef NSCOOL_H_
#define NSCOOL_H_

#include <iostream>
#include <vector>

enum State {
    SF_1s0,
    SF_3p2
};

class NSCool {
  public:
    NSCool();

    bool LoadEos(std::string eos);
    void PrintEosNames();

    double GetRMax() const {
        return std::min(raw_rTc.back(), raw_rT.back()) - 1e-10;
    }

    double GetT(double r) {
        return lerp(raw_rT, raw_T, r);
    }

    double GetTcn(double r) {
        return lerp(raw_rTc, raw_Tcn, r);
    }

    std::vector<double> raw_rTc;
    std::vector<double> raw_Tcn;
    std::vector<double> raw_Tcp;
    std::vector<State> raw_state;

    std::vector<double> raw_rT;
    std::vector<double> raw_T;
    std::vector<double> raw_ephi;

  private:
    void LocateEos();
    bool LoadCriticalTemps(std::string path);
    bool LoadProfile(std::string path);
    std::vector<std::string> Split(std::string line);
    double lerp(std::vector<double>& x, std::vector<double>& y, double nx);

    const std::string eos_parent_dir = "./eos/";
    std::vector<std::string> eos_names;
};

#endif // NSCOOL_H_
