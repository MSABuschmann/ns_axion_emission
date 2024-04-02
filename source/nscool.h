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

    std::vector<double> GetT() const {
        return raw_T;
    }

    std::vector<double> GetTc() const {
        return raw_Tcn;
    }

    std::vector<double> GetEphi() const {
        return raw_ephi;
    }

  private:
    void LocateEos();
    bool LoadCriticalTemps(std::string path);
    bool LoadProfile(std::string path);
    std::vector<std::string> Split(std::string line);

    const std::string eos_parent_dir = "./eos/";
    std::vector<std::string> eos_names;

    std::vector<double> raw_rTc;
    std::vector<double> raw_Tcn;
    std::vector<double> raw_Tcp;
    std::vector<State> raw_state;

    std::vector<double> raw_rT;
    std::vector<double> raw_T;
    std::vector<double> raw_ephi;
};

#endif // NSCOOL_H_
