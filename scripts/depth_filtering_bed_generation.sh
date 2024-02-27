#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/mosdepth/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/mosdepth/slurm-%j.out

# This script generates a bed file that contains the mean depth of consecutive 10kbp windows from a set of bam files with the tool mosdepth.

# Example of the mosdepth command:
#   mosdepth -t 4 -n --by 10000 <output_prefix> <input_bam>

# Usage of the script simultaneously for several bam files:
#    for input_bam in $(ls /bams/directory/*.bam); do 
#        job_id=(sbatch -c 16 -mem=20GB -t 00:30:00 depth_filtering_bed_generation.sh <${input_bam}> <output_directory> | awk '{print $4}')
#        echo "${job_id} ${input_bam}" >> /logs/directory/job_ids_depth_filtering_bed_generation.txt
#    done


# Load the mosdepth module
module load mosdepth

# Define the input bam file basename:
basename_bam=$(basename ${1} .bam)

# Define the output directory:
output_dir=${2}


# Run mosdepth:
mosdepth -t 4 -n --by 10000 ${output_dir}/${basename_bam}_mosdepth ${1}

