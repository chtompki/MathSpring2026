#!/bin/bash
#SBATCH --job-name=rk4_convergence
#SBATCH --output=rk4_convergence.out
#SBATCH --error=rk4_convergence.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

module load matlab/R2023b

matlab -batch "rk4_challenge"
