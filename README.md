ezh5: Easy HDF5 C++ Library
===========================

An easy to use c++ wrapper for HDF5 library

Features
--------

- Single-header library, just copy into your source code folder
- Very simple API so that knowledge of HDF5 is not required
- Work with std::complex which the original library does not provide
- Easy to mix with original library for customized need
- 

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

Compare with original library
-----------------------------


### Save an integer

#### original library

steals from [tutorial](http://www.hdfgroup.org/HDF5/Tutor/index.html) of HDF5 ([h5crtdat.c](http://www.hdfgroup.org/ftp/HDF5/current/src/unpacked/examples/h5crtdat.c) [h5_rdwt.c](http://www.hdfgroup.org/ftp/HDF5/current/src/unpacked/examples/h5_rdwt.c))

```C++

	hid_t file, dataset_id, dataspace_id;
	herr_t status;
	int data = 7;
	hsize_t dim = 1;
	
	file_id = H5Fcreate("test.h5", H5F_ACC_TRUNC);
	dataspace_id = H5Screate_simple(1, &dim, NULL);
	dataset_id = H5Dcreate2(file_id, "/dset", H5T_STD_I32BE, dataspace_id,
	                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
	
	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);
	status = H5Fclose(file_id);
	
```

#### ezh5

```C++
	File fh5 ("test.h5", H5F_ACC_TRUNC);
	fh5["dset"] = 7;
	// not need to close, resource is released in destructors of intermediate classes
```


Implementation Details for Other Developers
-------------------------------------------

Two teniques (template and lazy evaluation) to simplify this implementation.

I use template to provide type traits and construct composited datatype such as complex number.

Lazy evaluation is used to improve the performance, and it happens to greatly simplify creating groups.
