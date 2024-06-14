# Raw reads download

We decided to use 3 species as outgroups: tiger-cat-lynx rufus.

We already had data available for lynx rufus, but we had to download raw reads from the NCBI SRA for the other two species. 
I chose Illumina paired-end read data used to build the reference genomes of both species: [SRR5055389](https://www.ncbi.nlm.nih.gov/sra?LinkName=biosample_sra&from_uid=5980360) for cat, and [SRR17297154](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR17297154&display=metadata) for tiger.

I run a FastQC and fastp for both fastq files.

## Alignments

