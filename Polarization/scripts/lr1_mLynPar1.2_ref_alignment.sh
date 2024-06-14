#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/alignments/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/alignments/slurm-%j.out

set -e
module load bwa
module load samtools
module load picard
module load gatk/3.7-0-gcfedb67

# align reads from files 'LR1_R1.fastp.fastq.gz' and 'LR1_R2.fastp.fastq.gz':
bwa mem /mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.scaffolds.revcomp.scaffolds.fa /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/FASTQ_files/MAGROGEN/fastp/LR1_R1.fastp.fastq.gz /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/FASTQ_files/MAGROGEN/fastp/LR1_R2.fastp.fastq.gz -t 20 | samtools view -hbS -@ 20 - -o /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/polarization/lr1_mLynPar1.2_ref.bam

# sort the resulting bam file:
samtools sort -@ 20 /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/polarization/lr1_mLynPar1.2_ref.bam -o /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/polarization/lr1_mLynPar1.2_ref_sorted.bam

# add read groups to the sorted bam file:
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/polarization/lr1_mLynPar1.2_ref_sorted.bam O=/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/polarization/lr1_mLynPar1.2_ref_sorted_rg.bam RGID=macrogen RGLB=lr1_lib RGPL=Illumina RGPU=macrogen RGSM=lr1 VALIDATION_STRINGENCY=SILENT

# lr1 only has one fastqid!
# renaming its only sorted_rg BAM to merged_sorted:
mv /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/polarization/lr1_mLynPar1.2_ref_sorted_rg.bam /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/polarization/lr1_mLynPar1.2_ref_sorted_rg_merged_sorted.bam

# marking duplicates of merged bam:
java -jar $EBROOTPICARD/picard.jar MarkDuplicates METRICS_FILE=/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/polarization/lr1_rmdup.txt I=/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/polarization/lr1_mLynPar1.2_ref_sorted_rg_merged_sorted.bam O=/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/polarization/lr1_mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup.bam MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800

# index the duplicate marked bam for gatk:
samtools index /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/polarization/lr1_mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup.bam

# identify indel realignment targets:
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 20 -R /mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.scaffolds.revcomp.scaffolds.fa -I /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/polarization/lr1_mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup.bam -o /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/polarization/lr1_realignertargetcreator.intervals

# realign indels:
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -nt 20 -R /mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.scaffolds.revcomp.scaffolds.fa -targetIntervals /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/polarization/lr1_realignertargetcreator.intervals -I /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/polarization/lr1_mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup.bam -o /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/polarization/lr1_mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner.bam

# index the indel realigned bam for future analyses:
samtools index /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/polarization/lr1_mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner.bam
