#!/bin/bash
#SBATCH -N 8
#SBATCH --exclusive 
#SBATCH --time=8:00:00
# #SBATCH --time=00:10:00
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
echo $SLURM_NPROCS > nprocs
echo $SLURM_NNODES > nnodes

nproc=$(expr '32' '*' "$SLURM_NNODES")
##module add openmpi/intel
rm -f $casename.sch

mpprun ./nek5000 > out.$nproc

rm -f $casename.sch



