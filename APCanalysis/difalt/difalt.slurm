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
./../../bin/optiso finder_out/mod_1.fa search_out/mod_1.gff3 ../../models/worm.splicemodel > mod_1.optiso
./../../bin/optiso finder_out/dyn_1.fa search_out/dyn_1.gff3 ../../models/worm.splicemodel > dyn_1.optiso
./../../bin/optiso finder_out/erp_1.fa search_out/erp_1.gff3 ../../models/worm.splicemodel > erp_1.optiso
./../../bin/optiso finder_out/pdr_1.fa search_out/pdr_1.gff3 ../../models/worm.splicemodel > pdr_1.optiso
./../../bin/optiso finder_out/rbf_1.fa search_out/rbf_1.gff3 ../../models/worm.splicemodel > rbf_1.optiso
./../../bin/optiso finder_out/sad_1.fa search_out/sad_1.gff3 ../../models/worm.splicemodel > sad_1.optiso
./../../bin/optiso finder_out/unc_16.fa search_out/unc_16.gff3 ../../models/worm.splicemodel > unc_16.optiso
