#!/bin/bash
#SBATCH -c 4 # number of cores
#SBATCH --mem 5G # memory pool for all cores

module load SequenceAnalysis/SequenceSearch/diamond/2.0.9
module load Emboss/EMBOSS/6.6.0
module load R

sh ./running_diamond.shell ./data/Species.info.xls . diamond 4

Rscript prepare_i-ADHoRe.v2.R \
	-i ./data/Species.info.xls \
	-o . \
	-c run_iadhore.sh
