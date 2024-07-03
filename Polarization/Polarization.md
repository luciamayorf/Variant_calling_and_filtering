
# Synteny

## Raw reads download

We decided to use 3 species as outgroups: cat-marbled_cat-bobcat.

We already had data available for lynx rufus, but we had to download raw reads from the NCBI SRA for the other two species. 
I chose Illumina paired-end read sequencing data used to build the reference genomes for the cat (Felis catus): [SRR5055389](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR5055389&display=metadata) for the cat. For the marbled cat (Pardofelis marmorata), I downloaded the data generated [by Li et al. (2019)](https://academic.oup.com/mbe/article/36/10/2111/5518928): [SRR6071637](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR6071637&display=metadata).

We first downloading the data to EBD genomics server, extract the fastq files and perform a fastQC:
```bash
# For cat
cd /GRUPOS/grupolince/rawdata/felis_catus_SRA
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5055389/SRR5055389
fastq-dump --gzip --split-files SRR5055389
fastqc -o fastqc -t 16 SRR5055389_1.fastq.gz SRR5055389_2.fastq.gz

# For marbled cat
cd /GRUPOS/grupolince/rawdata/pardofelis_marmorata_SRA
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR6071637/SRR6071637
fastq-dump --gzip --split-files SRR6071637
fastqc -o fastqc -t 16 SRR6071637_1.fastq.gz SRR6071637_2.fastq.gz
```

Then I run fastp for both fastq files, as in [Data_preprocessing_alignment](https://github.com/luciamayorf/Data_preprocessing_alignment_v2) repository, adding the flag "--qualified_quality_phred 25" to the fastp command to remove lower quality reads. We performed another FastQC of the trimmed data and everything seemed correct.


## Alignments

I aligned the cat and marbled cata reads to the new Iberian lynx reference genome in the EBD genomics server and the lynx rufus data in CESGA ft3, running scripts [fc_SRA_mLynPar1.2_ref_alignment.sh](https://github.com/luciamayorf/Variant_calling_and_filtering/blob/main/Polarization/scripts/fc_SRA_mLynPar1.2_ref_alignment.sh), [pm_SRA_mLynPar1.2_ref_alignment.sh](https://github.com/luciamayorf/Variant_calling_and_filtering/blob/main/Polarization/scripts/pm_SRA_mLynPar1.2_ref_alignment.sh) and [lr1_mLynPar1.2_ref_alignment.sh](https://github.com/luciamayorf/Variant_calling_and_filtering/blob/main/Polarization/scripts/lr1_mLynPar1.2_ref_alignment.sh).

We also performed a quality control of the alignments using Qualimap bamqc, as in [Data_preprocessing_alignment](https://github.com/luciamayorf/Data_preprocessing_alignment_v2) repository.

## Base synteny file

The next step consists of generating a fasta file for each species with the Iberian lynx reference genome, but containing the alleles for each species in the variants positions of the pardinus VCF. For that, we use samtools pileup and put2fa.







