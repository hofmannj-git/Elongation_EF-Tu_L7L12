#!/bin/bash
#
#SBATCH --job-name=run_ssf
#
#SBATCH -p rzia
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=hofmannj@stanford.edu
#SBATCH --mail-type=FAIL,END

module load gcc/10.1.0
export LD_LIBRARY_PATH=/home/users/hofmannj/bin/hdf5/lib:$LD_LIBRARY_PATH

~/bin/fpostproc/src/plot_ssf ./ 51000000 200000 10000 8.0 802 0.1 10.0
