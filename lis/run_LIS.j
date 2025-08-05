#!/bin/csh
#SBATCH -J LIS_v1  #job name
#SBATCH -o slurm_LIS.o%j
#SBATCH -A forman-prj-aac
#SBATCH -t 7:00:00
#SBATCH -N 1
#SBATCH -n 1
####SBATCH --ntasks=112  try this as well!!
####SBATCH --mem=60G
#SBATCH --mail-user=alirezam@umd.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


. ~/.bash_profile

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
module list

cd make 
make realclean
cd ..

sh ./lis-build-zaratan
# sh ./makebuild.glue_mpi
# Echo job information
# echo "JobID is $SLURM_JOBID"
# scontrol show job $SLURM_JOBID
# End of script
# exit

