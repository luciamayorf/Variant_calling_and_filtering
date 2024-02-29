#!/bin/bash

#SBATCH --job-name=repeatmasker_lp_lib
#SBATCH --output=/mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/repetitive_regions/repeatmasker_lowcomplex.out
#SBATCH --error=/mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/repetitive_regions/repeatmasker_lowcomplex.err
#SBATCH --time=1-00:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=36

module load repeatmasker

RepeatMasker -pa 9 -gff -xsmall -noint -dir /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynPar1.2/repetitive_regions/low_complex/ mLynPar1.2.scaffolds.revcomp.scaffolds.fa
