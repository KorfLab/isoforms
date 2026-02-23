#!/bin/bash
#SBATCH --job-name=dif_genes            # jobs should have unique names
#SBATCH --account=publicgrp            # use publicgrp by default
#SBATCH --partition=low                # use low by default
#SBATCH --output=jobs/dif_genes.stdout  # relative path starts in $HOME
#SBATCH --error=jobs/dif_genes.stderr   # relative path starts in $HOME
#SBATCH --nodes=1                      # generally we request 1
#SBATCH --ntasks-per-node=1            # 1 for simple tasks
#SBATCH --cpus-per-task=7             # depends on job
#SBATCH --mem=14G                      # depends on job
#SBATCH --time=5:00                    # days-hours:min:sec

#--- You will probably do these most of the time ---#
# module load conda
# conda activate something

# Various commands
./../../bin/optiso genes/mod1_750bp.fa genes/mod1_750bp.gff3 ../../models/worm.splicemodel --name mod1_750bp.optiso
./../../bin/optiso genes/dyn1_750bp.fa genes/dyn1_750bp.gff3 ../../models/worm.splicemodel --name dyn1_750bp.optiso
./../../bin/optiso genes/erp1_750bp.fa genes/erp1_750bp.gff3 ../../models/worm.splicemodel --name erp1_750bp.optiso
./../../bin/optiso genes/pdr1_783bp.fa genes/pdr1_783bp.gff3 ../../models/worm.splicemodel --name pdr1_783bp.optiso
./../../bin/optiso genes/rbf1_750bp.fa genes/rbf1_750bp.gff3 ../../models/worm.splicemodel --name rbf1_750bp.optiso
./../../bin/optiso genes/sad1_750bp.fa genes/sad1_750bp.gff3 ../../models/worm.splicemodel --name sad1_750bp.optiso
./../../bin/optiso genes/unc16_750bp.fa genes/unc16_750bp.gff3 ../../models/worm.splicemodel --name unc16_750bp.optiso
