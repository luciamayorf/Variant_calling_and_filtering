#!/bin/bash
#SBATCH --output=/mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/variantqc/slurm-%j.out
#SBATCH --error=/mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/variantqc/slurm-%j.err

# This script calculates different statitics for a VCF file with vcftools and plink.

# Usage of the script: bash genotypes_qc_vcftools.sh <input_vcf> <output_directory>

# input vcf
input_vcf=${1}

# base name of the input vcf
base_name=$(basename ${input_vcf} .vcf)

# output directory
output_dir=${2}
mkdir -p ${output_dir}

# Load the modules
module load samtools
module load plink

# calculate allele frequency
vcftools --vcf ${input_vcf} --freq2 --out ${output_dir}/${base_name} --max-alleles 2

# calculate mean depth per individual
vcftools --vcf ${input_vcf} --out ${output_dir}/${base_name} --depth

# mean depth per site
vcftools --vcf ${input_vcf} --out ${output_dir}/${base_name} --site-mean-depth

# calculate site quality
vcftools --vcf ${input_vcf} --out ${output_dir}/${base_name} --site-quality

# calculate proportion of missing data per individual
vcftools --vcf ${input_vcf} --out ${output_dir}/${base_name} --missing-indv

# calculate proportion of missing data per site
vcftools --vcf ${input_vcf} --out ${output_dir}/${base_name} --missing-site

# calculate heterozygosity and inbreeding coefficient per individual
vcftools --vcf ${input_vcf} --out ${output_dir}/${base_name} --het

# calculate number of singletons 
vcftools --vcf ${input_vcf} --singletons --out ${output_dir}/${base_name}

# calculate number of singletons per individual
cut -f5 ${output_dir}/${base_name}.singletons | sort -r | uniq -c | sed 's/ \+ /\t/g' | cut -f2 | sed 's/ /\t/' > ${output_dir}/${base_name}.singletons_per_ind

# perform a PCA (setting half-calls to missing)
plink --vcf ${input_vcf} --double-id --allow-extra-chr --set-missing-var-ids @:# --vcf-half-call m --pca --out ${output_dir}/${base_name}


