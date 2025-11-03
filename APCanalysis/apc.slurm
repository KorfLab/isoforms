#!/bin/bash
#SBATCH --job-name=example1
#SBATCH --account=publicgrp
#SBATCH --partition=low
#SBATCH --output=jobs/example1.stdout # relative path starts in $HOME
#SBATCH --error=jobs/example1.stderr # relative path starts in $HOME
##SBATCH --nodes=1
##SBATCH --ntasks-per-nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10M
#SBATCH --time=1:05 # days-hours:min:sec

# optional #
# module load conda
# conda activate something

# Various commands
ls $HOME
ls $HOME | wcount
