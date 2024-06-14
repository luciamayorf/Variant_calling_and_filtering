
# Synteny

## Raw reads download

We decided to use 3 species as outgroups: tiger-cat-lynx rufus.

We already had data available for lynx rufus, but we had to download raw reads from the NCBI SRA for the other two species. 
I chose Illumina paired-end read sequencing data used to build the reference genomes in both species: [SRR5055389](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR5055389&display=metadata) for cat, and [SRR17297154](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR17297154&display=metadata) for tiger.

Downloading the data to EBD genomics server:
```bash
# For cat
cd /GRUPOS/grupolince/rawdata/felis_catus_SRA
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5055389/SRR5055389
fastq-dump --gzip --split-files SRR5055389 

# For tiger
cd /GRUPOS/grupolince/rawdata/panthera_tigris_SRA
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR17297154/SRR17297154
fastq-dump --gzip --split-files SRR17297154
```

Then I run a FastQC and fastp for both fastq files, as in [Data_preprocessing_alignment](https://github.com/luciamayorf/Data_preprocessing_alignment_v2) repository, adding the flag "--qualified_quality_phred 25" to the fastp command to remove lower quality reads.

## Alignments

I aligned cat and tiger data to the new Iberian lynx reference genome in the EBD genomics server and rufus data in CESGA ft3, running scripts [fc_SRA_alignment.sh](), [pt_SRA_alignment.sh]() and [lr1_mLynPar1.2_ref_alignment.sh](https://github.com/luciamayorf/Variant_calling_and_filtering/blob/main/Polarization/scripts/lr1_mLynPar1.2_ref_alignment.sh).

SUBIR LOS SCRIPST DE GATO Y TIGRE CUANDO COMPRUEBE QUE FUNCIONAN

## Base synteny file







