#!/bin/bash

#SBATCH --job-name=variant_filtering
#SBATCH --output=/mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/variant_filtering/slurm-%j.out
#SBATCH --error=/mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/variant_filtering/slurm-%j.err

# Purpose: Same filters as 1 to 5, but keeping SNPs with QUAL<20 instead of QUAL>30
# Input:
#   - $1: reference genome
#   - $2: input vcf
#   - $3: mask bed file (repetitive and low mappability regions)
# Output:
#   - $2.filter5.vcf: filtered vcf

# Usage:
#   - bash variant_filter_1to5_QUAL20.sh <ref.fa> <in.vcf> 


# calling the modules:
module load gatk/4.2.0.0
module load bedtools/2.31.0
module load samtools/1.14


#### 0. Variables definition ####

# fasta reference genome
ref=${1}

# input vcf
invcf=${2}

# basename of input vcf
vcf_basename=$(basename ${invcf} .vcf.gz)

# path of input vcf
vcf_dir=$(dirname ${invcf})

# prefix of vcfs
vcf_pre="${vcf_dir}/${vcf_basename}"



#### 4. Removing sites with low quality #####

echo "starting step 4: filtering out sites with low genotype quality (QUAL>=20)"

# Filter for QUAL field:
bcftools view -i 'QUAL>=20' -Ov -o ${vcf_pre}.filter4_QUAL20.vcf ${vcf_pre}.filter3.vcf



#### 5. Removing indels #####

echo "starting step 5: filtering out INDELS"

# Filter to end up with only SNPs:
gatk SelectVariants \
  -select-type SNP \
  -R ${ref} \
  -V ${vcf_pre}.filter4_QUAL20.vcf \
  -O ${vcf_pre}.filter5_QUAL20.vcf