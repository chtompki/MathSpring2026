#!/bin/bash
#SBATCH --job-name=rk4_ivp
#SBATCH --output=rk4_%j.out
#SBATCH --error=rk4_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
set -e
cd "$SLURM_SUBMIT_DIR"
module load matlab/R2024a
export RK4_QUICK=0
matlab -singleCompThread -batch "rk4_ivp_solution"

