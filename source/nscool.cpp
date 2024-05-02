#include <filesystem>
#include <fstream>
#include <iostream>

#include "nscool.h"
#include "process.h"
#include "utils.h"

NSCool::NSCool() { LocateEos(); }

void NSCool::LocateEos() {
    for (const auto &dir :
         std::filesystem::directory_iterator(eos_parent_dir)) {
        std::string dir_path = dir.path();
        dir_path.erase(0, eos_parent_dir.length());
        eos_names.push_back(dir_path);
    }
}

bool NSCool::LoadEos(std::string eos) {
    std::string eos_path = eos_parent_dir + eos;
    if (std::any_of(eos_names.begin(), eos_names.end(),
                    [eos](std::string name) { return eos == name; })) {
    } else {
        return false;
    }

    if (!LoadCriticalTemps(eos_path)) {
        return false;
    }

    if (!LoadProfile(eos_path)) {
        return false;
    }

    DetermineBoundaries();
    PrintBoundaries();
    return true;
}

bool NSCool::LoadCriticalTemps(std::string path) {
    std::ifstream file;
    file.open(path + "/Star_Try.dat");
    if (!file.is_open()) {
        return false;
    }

    std::string line;
    for (int i = 0; i < 6; ++i)
        std::getline(file, line);

    while (file) {
        std::getline(file, line);
        std::vector<std::string> entries = Split(line);

        switch (entries.size()) {
        case 24:
            raw_rTc.push_back(std::stod(entries[1]));
            raw_kfp.push_back(std::stod(entries[8]));
            raw_kfn.push_back(std::stod(entries[9]));
            raw_Tcn.push_back(std::stod(entries[14]));
            raw_Tcp.push_back(std::stod(entries[15]));
            raw_mstn.push_back(std::stod(entries[22]));
            raw_mstp.push_back(std::stod(entries[23]));
            raw_state.push_back(GetState(entries[18]));
            break;
        case 15:
            raw_rTc.push_back(std::stod(entries[1]));
            raw_kfp.push_back(0);
            raw_kfn.push_back(std::stod(entries[7]));
            raw_Tcn.push_back(std::stod(entries[8]));
            raw_Tcp.push_back(0);
            raw_mstn.push_back(std::stod(entries[13]));
            raw_mstp.push_back(std::stod(entries[14]));
            raw_state.push_back(GetState(entries[9]));
            break;
        default:
            continue;
        }
    }

    return true;
}

bool NSCool::LoadProfile(std::string path) {
    std::ifstream file;
    file.open(path + "/Temp_Try.dat");
    if (!file.is_open()) {
        return false;
    }

    std::string line;
    for (int i = 0; i < 9; ++i)
        std::getline(file, line);

    while (file) {
        std::getline(file, line);
        std::vector<std::string> entries = Split(line);
        if (entries.size() != 7) {
            continue;
        }
        raw_rT.push_back(std::stod(entries[1]));
        raw_T.push_back(std::stod(entries[5]));
        raw_ephi.push_back(std::stod(entries[3]));
        raw_dvol.push_back(std::stod(entries[4]));
    }

    raw_rT.push_back(0);
    raw_T.push_back(raw_T.back());
    raw_ephi.push_back(raw_ephi.back());
    raw_dvol.push_back(raw_dvol.back());

    std::reverse(raw_rT.begin(), raw_rT.end());
    std::reverse(raw_T.begin(), raw_T.end());
    std::reverse(raw_ephi.begin(), raw_ephi.end());
    std::reverse(raw_dvol.begin(), raw_dvol.end());

    raw_dvdr.resize(raw_rT.size());
    raw_dvdr[0] = 0;
    for (size_t i = 1; i < raw_rT.size(); ++i) {
        raw_dvdr[i] = raw_dvol[i] / (raw_rT[i] - raw_rT[i - 1]);
    }

    return true;
}

std::vector<std::string> NSCool::Split(std::string line) {
    std::stringstream ss(line);
    std::string s;
    std::vector<std::string> vs;
    while (getline(ss, s, ' ')) {
        if (s != "") {
            vs.push_back(s);
        }
    }

    return vs;
}

void NSCool::PrintEosNames() {
    for (const std::string &str : eos_names) {
        std::cout << str << "\n";
    }
}

void NSCool::DetermineBoundaries() {
    r_max = std::min(raw_rTc.back(), raw_rT.back()) - 1e-10;
    boundary_3p2 = FindDownwards(raw_state, 1.5);
    boundary_1s0n_lower = FindUpwards(raw_state, 0.5);
    boundary_1s0n_upper = FindDownwards(raw_state, 0.5);
    boundary_1s0p_lower = FindUpwards(raw_Tcp, 2.);
    boundary_1s0p_upper = FindDownwards(raw_Tcp, 2.);
    boundary_core_crust = FindDownwards(raw_Tcp, 0.5);
    if (boundary_3p2 > 0) {
        boundary_1s0n_lower = boundary_3p2;
    }
}

void NSCool::DetermineDeltaTInfty(Process *process) {
    raw_DeltaT_infty.clear();
    for (size_t i = 0; i < raw_Tcn.size(); ++i) {
        double r = raw_rTc[i];
        double Tc = raw_Tcn[i];
        double T = GetT(r);
        double DeltaT = process->GetDeltaT(T, Tc);
        double DeltaT_infty = DeltaT * GetEphi(r);
        raw_DeltaT_infty.push_back(DeltaT_infty);
    }
}

std::vector<double> NSCool::GetResonanceLayer(double omega) {
    std::vector<double> res;
    double target = omega / 2.;
    for (size_t i = 0; i < raw_DeltaT_infty.size() - 1; ++i) {
        if ((raw_DeltaT_infty[i] - target) *
                (raw_DeltaT_infty[i + 1] - target) <=
            0) {
            std::vector<double> x = {raw_DeltaT_infty[i],
                                     raw_DeltaT_infty[i + 1]};
            std::vector<double> y = {raw_rTc[i], raw_rTc[i + 1]};
            if (x[0] > x[1]) {
                std::reverse(x.begin(), x.end());
                std::reverse(y.begin(), y.end());
            }
            res.push_back(lerp(x, y, target));
        }
    }
    return res;
}

double NSCool::lerp(std::vector<double> &x, std::vector<double> &y, double nx) {
    std::vector<double>::iterator it = std::upper_bound(x.begin(), x.end(), nx);
    double t = (*it - nx) / (*it - *(it - 1));
    size_t i = std::distance(x.begin(), it);
    return std::lerp(y[i], y[i - 1], t);
}

void NSCool::WriteRawData(hid_t file_id) {
    WriteDataset(file_id, "raw_rTc", raw_rTc);
    WriteDataset(file_id, "raw_Tcn", raw_Tcn);
    WriteDataset(file_id, "raw_Tcp", raw_Tcp);
    WriteDataset(file_id, "raw_kfn", raw_kfn);
    WriteDataset(file_id, "raw_kfp", raw_kfp);
    WriteDataset(file_id, "raw_mstn", raw_mstn);
    WriteDataset(file_id, "raw_mstp", raw_mstp);
    WriteDataset(file_id, "raw_state", raw_state);
    WriteDataset(file_id, "raw_rT", raw_rT);
    WriteDataset(file_id, "raw_T", raw_T);
    WriteDataset(file_id, "raw_ephi", raw_ephi);
    WriteDataset(file_id, "raw_dvol", raw_dvol);
    WriteDataset(file_id, "raw_dvdr", raw_dvdr);
}

double NSCool::GetUnweightedMeanT() {
    double res = 0;
    for (size_t i = 0; i < raw_rTc.size(); ++i) {
        double T = lerp(raw_rT, raw_T, raw_rTc[i]);
        double ephi = lerp(raw_rT, raw_ephi, raw_rTc[i]);
        res += T * ephi;
    }
    return res / static_cast<double>(raw_rTc.size());
}

int NSCool::GetState(std::string state) {
    if (state == "no") {
        return 0;
    } else if (state == "1s0") {
        return 1;
    } else if (state == "3p2") {
        return 2;
    }
    return -1;
}

bool NSCool::GetBoundaries(double *rmin, double *rmax, std::string source) {
    if (source == "nn") {
        (*rmin) = 0;
        (*rmax) = r_max;
        return true;
    } else if (source == "np") {
        (*rmin) = 0;
        (*rmax) = boundary_core_crust;
        return true;
    } else if (source == "pp") {
        (*rmin) = 0;
        (*rmax) = boundary_core_crust;
        return true;
    } else if (source == "n") {
        if (boundary_1s0n_upper <= 0) {
            return false;
        }
        (*rmin) = boundary_1s0n_lower;
        (*rmax) = boundary_1s0n_upper;
        return true;
    } else if (source == "p") {
        if (boundary_1s0p_upper <= 0) {
            return false;
        }
        (*rmin) = boundary_1s0p_lower;
        (*rmax) = boundary_1s0p_upper;
        return true;
    } else if (source == "A") {
        if (boundary_3p2 <= 0) {
            return false;
        }
        (*rmin) = 0;
        (*rmax) = boundary_3p2;
        return true;
    } else if (source == "B") {
        if (boundary_3p2 <= 0) {
            return false;
        }
        (*rmin) = 0;
        (*rmax) = boundary_3p2;
        return true;
    }
    return false;
}

void NSCool::PrintBoundary(std::string process, std::string source) {
    double rmin = 0, rmax = 0;
    std::cout << "Boundary " << process << ": ";
    if (GetBoundaries(&rmin, &rmax, source)) {
        std::cout << rmin << " - " << rmax << " metres." << std::endl;
    } else {
        std::cout << "N/A." << std::endl;
    }
}

void NSCool::PrintBoundaries() {
    PrintBoundary("Bremsstrahlung nn", "nn");
    PrintBoundary("Bremsstrahlung np", "np");
    PrintBoundary("Bremsstrahlung pp", "pp");
    PrintBoundary("1s0(n)", "n");
    PrintBoundary("1s0(p)", "p");
    PrintBoundary("3p2 A", "A");
    PrintBoundary("3p2 B", "B");
}
