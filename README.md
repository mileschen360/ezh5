ezh5: Easy HDF5 C++ Library
===========================

An easy to use c++ wrapper for HDF5 library

Features
--------

- Single-header library, just copy into your source code folder
- Very simple API do not require knowledge of HDF5
- Work with std::complex which the original library does not provide
- Easy to mix with original library for customized need

Example
-------

```C++

#include "ezh5.hpp"
#include <complex>

using namespace ezh5;

int main(){
    File fh5 ("test.h5", H5F_ACC_TRUNC);
    fh5["x"]=2;
    fh5["group"]["y"] = 0.1;
    fh5["group"]["subgroup"]["z"] = 0.2;
    auto gr= fh5["group2"];
    gr["x2"] = 0.5;
    gr["y2"] = std::complex<double> (1.2, 3.5);
    return 0;
}



```
