#include "ezh5.hpp"
#include <complex>
#include <vector>

int main(){
    ezh5::File f ("universe.h5", H5F_ACC_TRUNC);  // create a file object called f
    f["length unit"] = "meter";
    f["time unit"] = "second";
    f["light speed"] = 3.0e8;
    f["solar system"]["mass unit"] = "kg";
    f["solar system"]["earth"]["mass"] = 1e9;
    ezh5::Node myhome = f["solar system"]["earth"];
    myhome["population"] = 157893294;  // why doesn't work here
    // std::vector<int> x(3);
    // fh5["vec"] = x;
    
    // fh5["group"]["y"] = 0.1;
    // fh5["group"]["subgroup"]["z"] = 0.2;
    // ezh5::Node gr= fh5["group2"];
    // gr["x2"] = 0.5;
    // gr["y2"] = std::complex<double> (1.2, 3.5);
    return 0;
}
