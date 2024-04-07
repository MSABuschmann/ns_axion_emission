#ifndef UTILS_H_
#define UTILS_H_

#include "hdf5.h"
#include <chrono>

#include "nscool.h"

template <typename T>
void WriteDataset(hid_t file_id, const char* name,
                  std::vector<T> dataset) {
    // Identify datatype.
    hid_t mem_type_id, dset_type_id;
    if (std::is_same<T, float>::value) {
        mem_type_id = H5T_NATIVE_FLOAT;
        dset_type_id = H5T_IEEE_F32LE;
    } else if (std::is_same<T, double>::value) {
        mem_type_id = H5T_NATIVE_DOUBLE;
        dset_type_id = H5T_IEEE_F64LE;
    } else if (std::is_same<T, int>::value) {
        mem_type_id = H5T_NATIVE_INT;
        dset_type_id = H5T_IEEE_F64LE;
    } else {
        std::cout << "Unknown datatype, cannot write: " << name << std::endl;
        return;
    }

    hsize_t dims[1] = {dataset.size()};
    hid_t space = H5Screate_simple(1, dims, NULL);
    hid_t dataset_id = H5Dcreate(file_id, name, dset_type_id, space,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             &dataset[0]);
    H5Dclose(dataset_id);
}

extern std::vector<double> CreateVector(double min, double max, int N);
extern void TestInterpolation(hid_t file_id, NSCool& nscool);
extern void PrintDuration(std::chrono::steady_clock::time_point& start_time);

#endif // UTILS_H_
