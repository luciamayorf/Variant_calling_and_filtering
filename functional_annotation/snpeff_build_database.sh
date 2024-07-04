#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/functional_annotation/snpeff/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/functional_annotation/snpeff/slurm-%j.out

# This scripts build the database for the snpEff tool. All the necessary files have already been allocated in the adequate directories (see markdown).

# Load the snpEff module
module load snpeff/5.0

# Build the database
java -jar $EBROOTSNPEFF/snpEff.jar build -d -gff3 -v LYPA1_2A -c /mnt/netapp1/Store_CSIC/home/csic/eye/lmf/snpEff/snpEff.config -dataDir /mnt/netapp1/Store_CSIC/home/csic/eye/lmf/snpEff/data
