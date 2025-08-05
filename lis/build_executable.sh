#!/bin/bash

# Load required modules
module purge
module load openmpi/gcc/9.4.0
module load openblas
module load netcdf 
module load netcdf-fortran
module load hdf4 
module load hdf5 
module load hdf-eos2 
module load jasper 
module load xerces-c
module load esmf eccodes gdal
module load python

# Current directory should be the LIS directory
cd /home/alirezam/scratch/LIS_library/LISF_GPS_added/lis

# Set up necessary environment variables (from lis-build-zaratan)
export LIS_ARCH=linux_gfortran
export LIS_FC=mpifort
export LIS_CC=mpicc
export LIS_MODESMF=$ESMF_INCDIR
export LIS_LIBESMF=$ESMF_LIBDIR
export LIS_ECCODES=$ECCODES_ROOT
export LIS_NETCDF=$NETCDF_FORTRAN_ROOT
export LIS_HDF4=$HDF4_ROOT
export LIS_HDF5=$HDF5_ROOT
export LIS_HDFEOS=$HDFEOS2_ROOT
export LIS_JASPER=$JASPER_ROOT
export LIS_LAPACK=$OPENBLAS_ROOT
export LIS_GDAL=$GDAL_ROOT
export LIS_LIBGEOTIFF=$GEOTIFF_ROOT

# Check if configure.lis exists, if not, run configure script
if [ ! -f "make/configure.lis" ]; then
    echo "Running configure script..."
    # Run configure with all the options that were used in lis-build-zaratan
    # Redirecting input from the here-document
    ./configure <<EOF
1
-1
2
2
0
1
2
1
1
9
1
1
1
0
0
0
0
0
0
EOF

    # Apply the same sed modifications from lis-build-zaratan
    sed -i "s/-L\$(LIB_NETCDF)/-L\$(LIB_NETCDF) -L\$(NETCDF_LIBDIR)/g" make/configure.lis
    sed -i "s/-lmpi/-lmpi -L\$(OPENBLAS_LIBDIR) -lopenblas/g" make/configure.lis
    sed -i "s/-I\$(INC_LIBGEOTIFF)/-I\$(INC_LIBGEOTIFF) -I\$(TIFF_INCDIR)/g" make/configure.lis
    sed -i "s/-L\$(LIB_LIBGEOTIFF)/-L\$(LIB_LIBGEOTIFF) -L\$(TIFF_LIBDIR)/g" make/configure.lis
    sed -i '/^LIB_ECCODES/s/lib\//lib64\//' make/configure.lis

    alldirs=(`echo OPENBLAS_LIBDIR LIB_NETCDF NETCDF_LIBDIR LIB_ESMF LIB_JASPER LIB_ECCODES LIB_HDF4 LIB_HDF5 LIB_HDFEOS LIB_GDAL LIB_FORTRANGIS LIB_LIBGEOTIFF TIFF_LIBDIR`)
    for dir in "${alldirs[@]}"; do
        sed -i "s/-L\$($dir)/-L\$($dir) -Wl,-rpath,\$($dir)/g" make/configure.lis
    done

    sed -i "s/fallow-argument-mismatch/Wargument-mismatch/g" make/Makefile
fi

# Option 1: Just link the main executable (fastest if you've already compiled all the needed .o files)
echo "Linking the LIS executable..."
make -f make/Makefile LIS

# Alternative options (uncomment if the first option doesn't work):

# Option 2: Only rebuild if source files changed, but don't clean
# echo "Building LIS without cleaning..."
# ./compile

# Option 3: If you need to rebuild a specific object file first
# echo "Rebuilding specific object file..."
# cd make
# make GPSdispObs_Mod.o
# cd ..
# make -f make/Makefile LIS

echo "Build completed!" 