//#include <sstream>
//#include <iostream>
//#include <fstream>
//#include <stdlib.h>
//#include "math.h"
//#include <omp.h>
//#include "hdf5.h"
#include <memory>

#include "nscool.h"

int main(int argc, const char* argv[]) {
    if (argc != 2) {
        std::cout << "No EOS provided!" << std::endl;
        return true;
    }
    std::string eos = argv[1];

    NSCool nscool;
    if (!nscool.LoadEos(eos)) {
        std::cout << "Could not load EOS '" << eos << "'!\n";
        std::cout << "Available EOSs are:\n";
        nscool.PrintEosNames();
    }

    return 0;
}
