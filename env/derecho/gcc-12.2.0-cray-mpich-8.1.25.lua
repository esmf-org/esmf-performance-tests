-- -*- lua -*-
-- Module file created by ESMF (https://earthsystemmodeling.org/) on 2024-12-03
--
-- gcc@13.2.0
-- cray-mpich@8.1.27
--

whatis([[Name : esmf-environment]])
whatis([[Compiler : gcc-12.2.0]])
whatis([[MPI : cray-mpich-8.1.25]])
whatis([[Python : python3-3.11.6]])
whatis([[Short description : Sets up ESMF build environment for Derecho.]])

help([[ESMF environment on derecho. Include intel-classic, cray-mpich, and
python.]])

-- Services provided by the package
family("esmfenv")

always_load("ncarenv/23.09")
always_load("craype/2.7.23")
always_load("gcc/13.2.0")
always_load("ncarcompilers/1.0.0")
always_load("hdf5/1.14.3")
always_load("netcdf/4.9.2")
always_load("cray-mpich/8.1.27")
always_load("cmake/3.26.3")

setenv("CC", "mpicc")
setenv("CXX", "mpicxx")
setenv("FC", "mpif90")
