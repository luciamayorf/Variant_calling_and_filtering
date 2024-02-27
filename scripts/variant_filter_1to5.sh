#!/bin/bash

#SBATCH --job-name=variant_filtering
#SBATCH --output=/mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/variant_filtering/slurm-%j.out
#SBATCH --error=/mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/variant_filtering/slurm-%j.err

# Purpose: Filter out repetitive regions, indels, non-biallelic sites, invariant sites, and sites with low quality.
# Input:
#   - $1: reference genome
#   - $2: input vcf
#   - $3: mask bed file (repetitive and low mappability regions)
# Output:
#   - $2.filter5.vcf: filtered vcf

# Usage:
#   - bash variant_filter_1to5.sh <ref.fa> <in.vcf> <masked_regions.bed>


# calling the modules:
module load gatk/4.2.0.0
module load bedtools/2.31.0
module load samtools/1.14


#### 0. Variables definition ####

# fasta reference genome
ref=${1}

# input vcf
invcf=${2}

# mask bed file (repetitive and low complexity regions)
mask=${3}

# basename of input vcf
vcf_basename=$(basename ${invcf} .vcf.gz)

# path of input vcf
vcf_dir=$(dirname ${invcf})

# prefix of vcfs
vcf_pre="${vcf_dir}/${vcf_basename}"



#### 1. Removing variants in repetitive/low_complexity regions ####

echo "starting step 1: filtering out repetitive regions"

# Apply the filter with BedTools subtract
bedtools subtract -a ${invcf} -b ${mask} -header | uniq > ${vcf_pre}.filter1.vcf


#### 2. Removing non-biallelic sites ####

echo "starting step 2: filtering out non-biallelic sites"

# Apply the filter with GATK SelectVariants
gatk SelectVariants \
  --restrict-alleles-to BIALLELIC \
  -R ${ref} \
  -V ${vcf_pre}.filter1.vcf \
  -O ${vcf_pre}.filter2.vcf

		# An alternative : opt/bcftools-1.6/bcftools view \
		# -m2 -M2 -Ov -o ${vcf_pre}.filter2.vcf ${vcf_pre}.filter1.vcf


#### 3. Removing invariant sites #####

echo "starting step 3: filtering out invariant sites"

# Apply the filter with BCFtools view
bcftools view \
  -e 'INFO/AF=1.00' ${vcf_pre}.filter2.vcf \
  > ${vcf_pre}.filter3.vcf
  
  

#### 4. Removing sites with low quality #####

echo "starting step 4: filtering out sites with low genotype quality (QUAL>=30)"

# Filter for QUAL field:
bcftools view -i 'QUAL>=30' -Ov -o ${vcf_pre}.filter4.vcf ${vcf_pre}.filter3.vcf



#### 5. Removing indels #####

echo "starting step 5: filtering out INDELS"

# Filter to end up with only SNPs:
gatk SelectVariants \
  -select-type SNP \
  -R ${ref} \
  -V ${vcf_pre}.filter4.vcf \
  -O ${vcf_pre}.filter5.vcf