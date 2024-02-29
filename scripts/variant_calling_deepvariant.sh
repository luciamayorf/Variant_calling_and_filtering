#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/novogene_lp_sept2023/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/novogene_lp_sept2023/slurm-%j.out

# This script calls variants from a set of bam files to obtain gVCFs with DeepVariant, using the model WGS.

# Example of the deepvariant command:
#   deepvariant_gpu --model_type=WGS --ref=/path/to/reference/genome \
#   --reads=/path/to/bam --output_gvcf=/path/to/output/output.g.vcf.gz --num_shards=32 --output_vcf=/path/to/output/output.vcf.gz

# Usage of the script simultaneously for several bam files:
#    for input_bam in $(ls /bams/directory/*.bam); do 
#        job_id=(sbatch -c 32 --mem=50GB -t 06:00:00 --cpus-per-task=32 --gres=gpu:a100:1 variant_calling_deepvariant.sh \
#               <${input_bam}> <ref_genome> <files_list> <output_directory> | awk '{print $4}')
#        echo "${job_id} ${input_bam}" >> /logs/directory/job_ids_deepvariant.txt
#    done

# Load the deepvariant module:
module load cesga/2020 deepvariant/1.6.0

# Define the input bam file:
bam=${1}

# Define reference genome:
ref_genome=${2}

# Define list of samples (FASTQ files):
list_file=${3}

# Define the output directory:
output_dir=${4}

# Define basename of bam:
basename_bam=$(basename ${bam} _mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner.bam)
echo "basename_bam: ${basename_bam}"

# Define final sample name:
sample_name=$(grep ${basename_bam} ${list_file} | cut -f2 | sort -u)
echo "sample_name: ${sample_name}"

# Define sample sex:
sex=$(grep ${basename_bam} ${list_file} | cut -f3 | sort -u )
echo "sample_name: ${sex}"

# Run deepvariant (changinf the command according to the individuals sex):
if [[ ${sex} == "male" ]]; then
    deepvariant_gpu --model_type=WGS \
                    --ref=${ref_genome} \
                    --reads=${bam} \
                    --output_gvcf=${output_dir}/${sample_name}_mLynPar1.2_ref.g.vcf.gz \
                    --output_vcf=${output_dir}/${sample_name}_mLynPar1.2_ref.vcf.gz \
                    --par_regions_bed="/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.PAR1_sexChr.bed" \
                    --haploid_contigs="mLynPar1.2_ChrX,mLynPar1.2_ChrY,mLynPar1.2_ChrY_unloc_1,mLynPar1.2_ChrY_unloc_2,mLynPar1.2_ChrY_unloc_3,mLynPar1.2_ChrY_unloc_4,mLynPar1.2_ChrY_unloc_5,mLynPar1.2_ChrY_unloc_6,mLynPar1.2_ChrY_unloc_7,mLynPar1.2_ChrY_unloc_8,mLynPar1.2_ChrY_unloc_9,mLynPar1.2_ChrY_unloc_10,mLynPar1.2_ChrY_unloc_11,mLynPar1.2_ChrY_unloc_12" \
                    --num_shards=32
else
    deepvariant_gpu --model_type=WGS \
                    --ref=${ref_genome} \
                    --reads=${bam} \
                    --output_gvcf=${output_dir}/${sample_name}_mLynPar1.2_ref.g.vcf.gz \
                    --output_vcf=${output_dir}/${sample_name}_mLynPar1.2_ref.vcf.gz \
                    --num_shards=32
fi
