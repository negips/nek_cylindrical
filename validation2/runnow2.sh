#!/bin/bash
# #SBATCH -N 1 
# #SBATCH --exclusive 
# # #SBATCH --time=72:00:00
# #SBATCH --time=00:15:00
# #SBATCH --exclusive
# #SBATCH -A snic2021-22-728
# #SBATCH --mail-type=ALL
# #SBATCH --mail-user=negi@mech.kth.se 
# #SBATCH --output=job.%J.out
# 
# #SBATCH -J impulse

casename=cyl

echo $casename > SESSION.NAME

echo $PWD/ >> SESSION.NAME

echo $SLURM_JOB_NODELIST > node
echo $SLURM_NNODES > nnodes

nproc=$(expr '32' '*' "$SLURM_NNODES")

echo $nproc > nprocs
##module add openmpi/intel
rm -f $casename.sch

mpprun -n 256 ./nek5000 > out.$nproc

rm -f $casename.sch

