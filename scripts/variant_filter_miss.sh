#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/variant_filtering/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/variant_filtering/slurm-%j.out

# This script remove the variants missing in more than 70% of the individuals.

# Usage example: sbatch variant_filter_miss.sh <input_vcf> 

# Load the samtools module
module load samtools/1.19

# Define input vcf
vcf=${1}

# Basename of input vcf
vcf_basename=$(basename ${vcf} .vcf)

# Path of input vcf
vcf_path=$(dirname ${vcf})

# Run bcftools
bcftools filter -e 'F_MISSING > 0.7' -Ob -o ${vcf_path}/${vcf_basename}.miss.vcf ${vcf}