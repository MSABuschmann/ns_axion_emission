#include <filesystem>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cmath>

#include "utils.h"
#include "nscool.h"

NSCool::NSCool() {
    LocateEos();
}

void NSCool::LocateEos() {
    for (const auto& dir: std::filesystem::directory_iterator(eos_parent_dir)) {
        std::string dir_path = dir.path();
        dir_path.erase(0, eos_parent_dir.length());
        eos_names.push_back(dir_path);
    }
}

bool NSCool::LoadEos(std::string eos) {
    std::string eos_path = eos_parent_dir + eos;
    if (std::any_of(eos_names.begin(), eos_names.end(), [eos](std::string name)
        { return eos == name; })) {} else {
        return false;
    }

    if (!LoadCriticalTemps(eos_path)) {
        return false;
    }

    if (!LoadProfile(eos_path)) {
        return false;
    }

    Determine1s03p2Boundary();
    
    std::cout << "3p2: 0m - " << Get1s03p2Boundary() << "m" << std::endl;
    std::cout << "1s0: " << Get1s03p2Boundary() << "m - " 
              << GetRMax() << "m" << std::endl;

    return true;
}

bool NSCool::LoadCriticalTemps(std::string path) {
    std::ifstream file;
    file.open(path+"/Star_Try.dat");
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
            case 22:
                raw_rTc.push_back(std::stod(entries[1]));
                raw_Tcn.push_back(std::stod(entries[14]));
                raw_Tcp.push_back(std::stod(entries[15]));
                raw_state.push_back(static_cast<int>(SF_3p2));
                break;
            case 13:
                raw_rTc.push_back(std::stod(entries[1]));
                raw_Tcn.push_back(std::stod(entries[8]));
                raw_Tcp.push_back(0);
                raw_state.push_back(static_cast<int>(SF_1s0));
                break;
            default: continue;
        }
    }

    return true;
}

bool NSCool::LoadProfile(std::string path) {
    std::ifstream file;
    file.open(path+"/Temp_Try.dat");
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
        size_t last = raw_rT.size() - 1;
        if (last>0) {
            double dr = raw_rT[last-1] - raw_rT[last];
            raw_dvdr.push_back(raw_dvol.back() / dr);
        } else {
            raw_dvdr.push_back(0);
        }
    }

    raw_rT.push_back(0);
    raw_T.push_back(raw_T.back());
    raw_ephi.push_back(raw_ephi.back());
    raw_dvol.push_back(raw_dvol.back());
    raw_dvdr[0] = raw_dvdr[1];
    raw_dvdr.push_back(raw_dvdr.back());

    std::reverse(raw_rT.begin(), raw_rT.end());
    std::reverse(raw_T.begin(), raw_T.end());
    std::reverse(raw_ephi.begin(), raw_ephi.end());
    std::reverse(raw_dvol.begin(), raw_dvol.end());
    std::reverse(raw_dvdr.begin(), raw_dvdr.end());

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
    for (const std::string& str: eos_names) {
        std::cout << str << "\n";
    }
}

void NSCool::Determine1s03p2Boundary() {
    size_t i = std::distance(raw_state.begin(),
                    std::lower_bound(raw_state.begin(), raw_state.end(), 1));
    boundary_1s03p2 = (raw_rTc[i] + raw_rTc[i-1])/2.;
}

double NSCool::lerp(std::vector<double>& x, std::vector<double>& y, double nx) {
    std::vector<double>::iterator it = std::upper_bound(x.begin(), x.end(), nx);
    double t = (*it - nx)/(*it - *(it-1));
    size_t i = std::distance(x.begin(), it);
    return std::lerp(y[i], y[i-1], t);
}

void NSCool::WriteRawData(hid_t file_id) {
    WriteDataset(file_id, "raw_rTc", raw_rTc);
    WriteDataset(file_id, "raw_Tcn", raw_Tcn);
    WriteDataset(file_id, "raw_Tcp", raw_Tcp);
    WriteDataset(file_id, "raw_state", raw_state);
    WriteDataset(file_id, "raw_rT", raw_rT);
    WriteDataset(file_id, "raw_T", raw_T);
    WriteDataset(file_id, "raw_ephi", raw_ephi);
    WriteDataset(file_id, "raw_dvol", raw_dvol);
    WriteDataset(file_id, "raw_dvdr", raw_dvdr);
}
