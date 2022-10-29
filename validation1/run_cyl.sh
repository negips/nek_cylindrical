#!/bin/bash
#SBATCH -N 8
#SBATCH --exclusive 
#SBATCH --time=16:00:00
# SBATCH --time=00:02:00
#SBATCH --exclusive
#SBATCH -A snic2022-5-297
#SBATCH --mail-type=ALL
#SBATCH --mail-user=prabal.negi@su.se 
#SBATCH --output=job.%J.out

#SBATCH -J cyl

casename=cyl

echo $casename > SESSION.NAME

echo $PWD/ >> SESSION.NAME

echo $SLURM_JOB_NODELIST > node

##module add openmpi/intel
rm -f $casename.sch

mpprun ./nek5000 > out.256

rm -f $casename.sch



