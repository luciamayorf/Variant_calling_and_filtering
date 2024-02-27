#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/variant_filtering/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/variant_filtering/slurm-%j.out

# This script remove the variants from a VCF file that are found in a BED file containing windows to exclude.

# Usage example: sbatch variant_filter_rd.sh <input_vcf> <input_bed> 

# Load the bedtools module
module load bedtools/2.31.0

# Define input vcf
vcf=${1}

# Define input bed
bed=${2}

# Basename of input vcf
vcf_basename=$(basename ${vcf} .vcf)

# Path of input vcf
vcf_path=$(dirname ${vcf})

# Run bedtools
bedtools subtract -header -a ${vcf} -b ${bed} > ${vcf_path}/${vcf_basename}_rd.vcf