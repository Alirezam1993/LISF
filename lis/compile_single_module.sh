#!/bin/bash

if [ $# -lt 1 ]; then
  echo "Usage: $0 <module_name>"
  echo "Example: $0 GPSdispObs_Mod"
  exit 1
fi

MODULE_NAME=$1

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

# Go to the make directory
cd make

# First, rebuild the specific module
echo "Rebuilding module: ${MODULE_NAME}.o"
make ${MODULE_NAME}.o

# Now link the full executable
echo "Linking the LIS executable..."
make -f Makefile LIS

# Copy the executable to the main directory if needed
if [ -f LIS ]; then
  echo "Copying executable to parent directory..."
  cp LIS ..
fi

cd ..
echo "Build completed!" 