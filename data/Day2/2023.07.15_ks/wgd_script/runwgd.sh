#!/bin/bash
#SBATCH -c 8 # number of cores
#SBATCH --mem 5G # memory pool for all cores

export OPENBLAS_NUM_THREADS=2
export GOTO_NUM_THREADS=2
export OMP_NUM_THREADS=2

module load SequenceAnalysis/SequenceSearch/diamond/2.0.9
module load SequenceAnalysis/StructurePrediction/mcl/14.137
module load SequenceAnalysis/MultipleSequenceAlignment/mafft/7.475
module load Phylogeny/paml/4.9j
module load Phylogeny/FastTree/2.1.10 
module load Conda/miniconda/latest

source activate wgd

wgd dmd -I 3 Elaeis.fa -o 01.wgd_dmd
wgd ksd 01.wgd_dmd/Elaeis.fa.mcl Elaeis.fa -o 02.wgd_ksd -n 8
wgd syn -f mRNA -a ID -ks 02.wgd_ksd/Elaeis.fa.ks.tsv Elaeis.gff 01.wgd_dmd/Elaeis.fa.mcl -o 03.wgd_syn
wgd mix -ni 100 --method bgmm -n 1 5 02.wgd_ksd/Elaeis.fa.ks.tsv -o 04.wgd_mix
