#!/bin/bash
#SBATCH -o gself_array.o%j
#SBATCH --time=36:00:00
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --mem=96gb

module load matlab

cd $SLURM_SUBMIT_DIR
matlab -batch "run('metabolic_challenge_autapse_ReLU.m')"
