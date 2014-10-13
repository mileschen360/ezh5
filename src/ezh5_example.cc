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
