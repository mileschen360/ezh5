#ifndef EZH5_H
#define EZH5_H

#include <hdf5.h>

#ifndef CANNOT_FIND_BOOST
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#endif

#include <complex>
#include <string>
#include <vector>



namespace ezh5{

    using namespace std;
    
    template<typename T>
    struct DataType{
	static hid_t id;
    };

    
    // define complex number datatype
    template<typename T>
    struct DataType<std::complex<T> >{
	static hid_t id;
	DataType(){
	    H5Tinsert(id, "r", 0, DataType<T>::id);  // the name 'r' is to be compatible with h5py
	    H5Tinsert(id, "i", sizeof(T), DataType<T>::id); // TODO: use HOFFSET
	}
    };

    
    /// char* dsname
    template<typename T>
    hid_t write(hid_t loc_id, const char* dsname, const T& buf){
	hid_t dataspace_id = H5Screate(H5S_SCALAR);
	hid_t dataset_id = H5Dcreate(loc_id, dsname, DataType<T>::id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t error_id = H5Dwrite(dataset_id, DataType<T>::id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buf);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);
	return error_id;
    }
    /// string dsname
    template<typename T>
    hid_t write(hid_t loc_id, const std::string& dsname, const T& buf){
	return write(loc_id, dsname.c_str(), buf);
    }

    // /// char* dsname     // TODO: use enable_if to distinguish this from above
    // template<typename VEC>
    // hid_t write(hid_t loc_id, const char* dsname, const VEC& vec){
    // 	hsize_t dims[1];
    // 	dims[0] = vec.size();
    // 	hid_t dp_id = H5Screate_simple(1, dims, NULL);
    // 	hid_t ds_id = H5Dcreate(loc_id, dsname, DataType<VEC::value_type>::id, dp_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // 	hid_t err_id = H5Dwrite(ds_id, DataType<VEC::value_type>::id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vec[0]);
    // 	H5Dclose(ds_id);
    // 	H5Sclose(dp_id);
    // 	return err_id;
    // }
    // /// string dsname
    // template<typename VEC>
    // hid_t write(hid_t loc_id, const std::string& dsname, const VEC& vec){
    // 	return write(loc_id, dsname.c_str(), vec);
    // }


    /// write std::vector to hdf5
    /// char* dsname
    template<typename T>
    hid_t write(hid_t loc_id, const char* dsname, const std::vector<T>& vec){
	hsize_t dims[1];
	dims[0] = vec.size();
	std::cout<<dsname<<std::endl;
	hid_t dp_id = H5Screate_simple(1, dims, NULL);
	assert(dp_id>=0);
	hid_t ds_id = H5Dcreate(loc_id, dsname, DataType<T>::id, dp_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	assert(ds_id>=0);
	hid_t err_id = H5Dwrite(ds_id, DataType<T>::id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vec[0]);
	H5Dclose(ds_id);
	H5Sclose(dp_id);
	return err_id;
    }
    /// string dsname
    template<typename T>
    hid_t write(hid_t loc_id, const std::string& dsname, const std::vector<T>& vec){
	return write(loc_id, dsname.c_str(), vec);
    }
    
#ifdef _BOOST_UBLAS_VECTOR_
    /// char* dsname
    template<typename T>
    hid_t write(hid_t loc_id, const char* dsname, const boost::numeric::ublas::vector<T>& vec){
	hsize_t dims[1];
	dims[0] = vec.size();
	hid_t dp_id = H5Screate_simple(1, dims, NULL);
	hid_t ds_id = H5Dcreate(loc_id, dsname, DataType<T>::id, dp_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t err_id = H5Dwrite(ds_id, DataType<T>::id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vec(0));
	H5Dclose(ds_id);
	H5Sclose(dp_id);
	return err_id;
    }

    /// string dsname
    template<typename T>
    hid_t write(hid_t loc_id, const std::string& dsname, const boost::numeric::ublas::vector<T>& vec){
	return write(loc_id, dsname.c_str(), vec);
    }
#endif
    

#ifdef _BOOST_UBLAS_MATRIX_
    template<typename T>
    hid_t write(hid_t loc_id, const char* dsname, const boost::numeric::ublas::matrix<T>& mat){
	hsize_t dims[2];
	dims[0] = mat.size1();
	dims[1] = mat.size2();
	hid_t dp_id = H5Screate_simple(2, dims, NULL);
	htri_t is_exist = H5Lexists(loc_id, dsname, H5P_DEFAULT);
	hid_t ds_id = -1;
	switch (is_exist){
	case true:
	    ds_id = H5Dopen1(loc_id, dsname);
	    break;
	case false:
	    ds_id = H5Dcreate(loc_id, dsname, DataType<T>::id, dp_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	    break;
	default:
	    break;
	}
	assert(ds_id >=0);
	hid_t err_id = H5Dwrite(ds_id, DataType<T>::id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &mat(0,0));
	assert(err_id>=0);
	H5Dclose(ds_id);
	H5Sclose(dp_id);
	return err_id;
    }
    /// string dsname
    template<typename T>
    hid_t write(hid_t loc_id, const std::string& dsname, const boost::numeric::ublas::matrix<T>& mat){
	return write(loc_id, dsname.c_str(), mat);
    }
#endif



//---------- read function    
    template<typename T>
    hid_t read(hid_t loc_id, const char* dsname, T* p_buf){
	hid_t dataset_id = H5Dopen2(loc_id, dsname, H5P_DEFAULT);
	hid_t datatype_id = H5Dget_type(dataset_id);
	hid_t error_id = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, p_buf);
	H5Dclose(dataset_id);
	return error_id;
    }
    
    
    // hid_t write_complex_vec(hid_t loc_id, const char* name, const std::size_t vec_size, const std::complex<double>* buf){
    // 	hsize_t dim=vec_size;
    // 	hid_t dataspace_id = H5Screate_simple(1, &dim, NULL);
    // 	hid_t dataset_id = H5Dcreate(loc_id, name, complex_tid, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // 	hid_t error_id = H5Dwrite(dataset_id, complex_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
    // 	H5Dclose(dataset_id);
    // 	H5Sclose(dataspace_id);
    // 	return error_id;
    // }

    // hid_t read_complex(hid_t loc_id, const char* name, std::complex<double>* buf){
    // 	hid_t dataset_id = H5Dopen2(loc_id, name, H5P_DEFAULT);
    // 	hid_t datatype_id = H5Dget_type(dataset_id);
    // 	hid_t error_id = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
    // 	H5Dclose(dataset_id);
    // 	return error_id;
    // }
    
}

namespace ezh5{

    template<> hid_t DataType<double>::id = H5T_NATIVE_DOUBLE;
    template<> hid_t DataType<float>::id = H5T_NATIVE_FLOAT;
    template<> hid_t DataType<int>::id = H5T_NATIVE_INT;
    template<> hid_t DataType<long>::id = H5T_NATIVE_LONG;
    template<> hid_t DataType<unsigned int>::id = H5T_NATIVE_UINT;
    template<> hid_t DataType<unsigned long>::id = H5T_NATIVE_ULONG;

    template<> hid_t DataType<std::complex<float> >::id = H5Tcreate(H5T_COMPOUND, sizeof(std::complex<float>));  // create compound datatype
    template<> hid_t DataType<std::complex<double> >::id = H5Tcreate(H5T_COMPOUND, sizeof(std::complex<double>));
    // TODO: why I need specification, instead of defining follow
    //template<typename T>  
    //hid_t DataType<std::complex<T> >::id = H5Tcreate(H5T_COMPOUND, sizeof(std::complex<T>));


    DataType<std::complex<float> > obj_to_run_constructor_float;  // to add LAYOUT (STRUCTURE) to compound datatype
    DataType<std::complex<double> > obj_to_run_constructor_double;

}

namespace ezh5{

    class ID {
    public:
	hid_t id;
	ID(): id(-1){}
	ID(hid_t id_in) : id(id_in){}

	~ID(){
	}
    
    };


    class Node :public ID{
    public:

	Node(){}
    
	Node(hid_t pid_in, const string& path_in)
	    : ID(-1),
	      pid(pid_in),
	      path(path_in){
	    cout<<"creating "<<path<<endl;
	}

    
	Node& operator()(const string& path){
	    return *this;
	}

    
	Node operator[](const string& path_more){
	    if(this->id==-1){// in lazy
		htri_t is_exist = H5Lexists(pid, path.c_str(), H5P_DEFAULT);
		if (is_exist<0){
		    assert(false);
		}else if (is_exist==false){
		    this->id = H5Gcreate2(pid, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		}else{
		    this->id = H5Gopen(pid, path.c_str(), H5P_DEFAULT);
		}
	    }
	    assert(this->id>0);
	    return Node(this->id, path_more);
	
	    // std::size_t i_slash = path.find('/');
	    // this->id = H5Gcreate2(pid, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	    // return Node(this->id, path_more);
	    // }else{
	    //     hid_t gid;
	    //     const string& grpname = path.substr(0,i_slash);
	    //     htri_t is_exist = H5Lexists(this->id, grpname.c_str(), H5P_DEFAULT);
	    //     if (is_exist<0){
	    // 	assert(false);
	    //     }else if (is_exist==false){
	    // 	gid = H5Gcreate2(this->pid, path.substr(0,i_slash).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	    //     }else{
	    // 	gid = H5Gopen(this->id, grpname.c_str(), H5P_DEFAULT);
	    //     }
	    //     assert(gid>=0);
	    //     return Node(gid, path.substr(i_slash+1, string::npos));
	    // }
	}

    
	template<typename T>
	Node& operator=(T val){
	    if(id==-1){
		hid_t dataspace_id = H5Screate(H5S_SCALAR);
		this->id = H5Dcreate(pid, path.c_str(), DataType<T>::id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	    }
	    hid_t error_id = H5Dwrite(id, DataType<T>::id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &val);
	    H5Dclose(this->id);
	    this->id = -1;
	    return *this;

	
	    // std::size_t i_slash;
	    // if((i_slash= path.find('/'))!=string::npos){
	    //     this->id = H5Gcreate2(pid, path.substr(0,i_slash).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	    //     Node(id, path.substr(i_slash+1, string::npos)) = val;
	    // }else{//leaf
	    //     if(id==-1){
	    // 	hid_t dataspace_id = H5Screate(H5S_SCALAR);
	    // 	this->id = H5Dcreate(pid, path.c_str(), DataType<T>::id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	    //     }
	    //     hid_t error_id = H5Dwrite(id, DataType<T>::id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &val);
	    // }
	}

	~Node(){
	    cout<<"closing "<<path<<endl;
	    if(this->id>0){
		H5Gclose(this->id);
		this->id = -1;
	    }
	}
    
    public:
	hid_t pid;
	string path;
    };

    class Dataset : public Node{


    };


    class File : public Node{
    public:
	// TODO: open an opened file
	// TODO: implement open if exists, create if not
	File(const string& path, unsigned flags){
	    if(flags==H5F_ACC_RDWR || flags==H5F_ACC_RDONLY){
		this->id = H5Fopen(path.c_str(), flags, H5P_DEFAULT);
	    }else if (flags==H5F_ACC_TRUNC || flags==H5F_ACC_EXCL){
		this->id = H5Fcreate(path.c_str(), flags, H5P_DEFAULT, H5P_DEFAULT);
	    }else{
		assert(false && "unknow file access mode");
	    }
	}


	~File(){
	    H5Fclose(this->id);
	    this->id = -1;
	}

	// Node operator[](const string& path){
	// 	std::size_t i_slash = path.find('/');
	// 	cout<<i_slash<<endl;
	// 	if (i_slash==string::npos)
	// 	    return Node(this->id, path);
	// 	else{
	// 	    hid_t gid;
	// 	    const string& grpname = path.substr(0,i_slash);
	// 	    htri_t is_exist = H5Lexists(this->id, grpname.c_str(), H5P_DEFAULT);
	// 	    if (is_exist<0){
	// 		assert(false);
	// 	    }else if (is_exist==false){
	// 		gid = H5Gcreate2(this->id, path.substr(0,i_slash).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// 	    }else{
	// 		gid = H5Gopen(this->id, grpname.c_str(), H5P_DEFAULT);
	// 	    }
	// 	    assert(gid>=0);
	// 	    return Node(gid, path.substr(i_slash+1, string::npos));
	// 	}
	// }
    };


}


#endif
