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

# First, rebuild the specific module(s) we modified
echo "Rebuilding GPS module..."
cd make
make GPSdispObs_Mod.o
make read_GPSdispObs.o
make enksgps_general.o
make enksgps_Mod.o
make read_GPSdispObs.o
make noahmp401_getgpsdisppred.o
make noahmp401_gpsdisp_DAlogMod.o
make LIS_lsmda_pluginMod.o
make noahmp401_setgpsdispvars.o
make retrospective_runMod.o
make module_sf_noahmplsm_401.o
make noahmp401_updategpsdisp.o
make LIS_DAobservationsMod.o
make merra2_forcingMod.o

# Now link the full executable with the updated objects (stay in the make directory)
echo "Linking the LIS executable..."
make -f Makefile LIS

# Copy the executable to the main directory if needed
if [ -f LIS ]; then
  echo "Copying executable to parent directory..."
  cp LIS ..
fi

cd ..
echo "Build completed!" 
noahmp401_getgpsdisppred.o
noahmp401_gpsdisp_DAlogMod.o
surfacemodels/land/noahmp.4.0.1/da_gpsdisp/noahmp401_updategpsdisp.o
surfacemodels/land/noahmp.4.0.1/da_gpsdisp/noahmp401_getgpsdisppred.o
dataassim/obs/GPS/read_GPSdispObs.o
make LIS_DAobservationsMod.o
dataassim/obs/GPS/GPSdispObs_Mod.o dataassim/obs/GPS/read_GPSdispObs.o
core/LIS_DAobservationsMod.o
core/LIS_DAobservationsMod.o
dataassim/obs/GPS/read_GPSdispObs.o
