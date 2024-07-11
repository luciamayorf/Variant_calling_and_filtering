# 0. Building SnpEff database

The Iberian lynx reference genome is not found in the SnpEff database. We will build the database using the reference genome and the annotation file. 
For that, we first need to modify the snpEff.config file to include the new database. As I don't have the permissions to write the file, I need to create a copy in my $STORE.


```bash
module load snpeff/5.0
cp /opt/cesga/2020/software/Core/snpeff/5.0/snpEff.config /mnt/netapp1/Store_CSIC/home/csic/eye/lmf/snpEff/
```

I need to add these lines at the end of the config file:
```bash
# Lynx_pardinus, version mLynPar1.2
LYPA1_2A.genome : Iberian lynx        # from now on, LYPA1_2A is the code for the Lynx pardinus reference genome (in snpEff)
```

### Longest isoform selection

To build the database, SnpEff needs all the files to be in the same folder (named after the reference genome code in the configuration file). It needs the following files: the reference genome FASTA, the GFF file, and a coding and protein sequences fasta).

To get one effect per variant, we need to select one transcript and one protein per gene. We will base out decision on the longest peptide that it's transcribed (and its correspondent transcript).

We have the FASTA with the protein sequence of the longest peptide available (provided by CNAG with the new reference genomes). We will start with this file to build the longest transcript file, after which we will generate the new gff3 that contains only the longest isoform.

```bash
cd /mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_pardinus_mLynPar1.2
# For the transcripts file
for i in $(zgrep ">" LYPA1_2A.longestpeptide.fa.gz | sed '/^>/ s/\(.*\)P/\1T/;s/>//'); do
    awk -v id="${i}" -v RS='>' -v ORS='' '$1 == id {print ">"$0}' LYPA1_2A.transcripts.fa
done > snpeff_annotation/LYPA1_2A.transcripts_longest_CESGA.fa

# To check that they match:
diff <(zgrep ">" LYPA1_2A.longestpeptide.fa.gz | sed '/^>/ s/\(.*\)P/\1T/' | sed 's/>//' | sort) <(grep ">" snpeff_annotation/LYPA1_2A.transcripts_longest_CESGA.fa | sed 's/>//' | sort)
```

To generate the new GFF3 file, I use the custom script [gff_longest_generation](https://github.com/luciamayorf/Variant_calling_and_filtering/blob/main/functional_annotation/gff_longest_generation.sh):
```bash
sbatch -t 04:00:00 --mem=500MB /home/csic/eye/lmf/scripts/functional_annotation/snpeff/gff_longest_generation.sh
```

Now that we have all the files ready, we move them to the same folder and rename them so that SnpEff can build the database
```bash
# create a directory inside the software's dependencies whose name matches the code
mkdir /mnt/netapp1/Store_CSIC/home/csic/eye/lmf/snpEff/data/LYPA1_2A

# copy the reference genome FASTA to the new folder and rename it following the manuals instructions
cp /mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.scaffolds.revcomp.scaffolds.fa /mnt/netapp1/Store_CSIC/home/csic/eye/lmf/snpEff/data/LYPA1_2A/
mv mLynPar1.2.scaffolds.revcomp.scaffolds.fa sequences.fa

# copy and rename the gff3 file
cp /mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_pardinus_mLynPar1.2/snpeff_annotation/LYPA1_2A.FA_longestisoform_CESGA.gff3 /mnt/netapp1/Store_CSIC/home/csic/eye/lmf/snpEff/data/LYPA1_2A/
mv /mnt/netapp1/Store_CSIC/home/csic/eye/lmf/snpEff/data/LYPA1_2A/LYPA1_2A.FA_longestisoform_CESGA.gff3 /mnt/netapp1/Store_CSIC/home/csic/eye/lmf/snpEff/data/LYPA1_2A/genes.gff

# copy and rename the transcripts file 
cp /mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_pardinus_mLynPar1.2/snpeff_annotation/LYPA1_2A.transcripts_longest_CESGA.fa /mnt/netapp1/Store_CSIC/home/csic/eye/lmf/snpEff/data/LYPA1_2A/
mv /mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_pardinus_mLynPar1.2/snpeff_annotation/LYPA1_2A.transcripts_longest_CESGA.fa /mnt/netapp1/Store_CSIC/home/csic/eye/lmf/snpEff/data/LYPA1_2A/cds.fa

# move and rename the protein file (the code for the protein needs to be the same as the code for the transcript, so we need to change the last "P" for a "T").
zcat LYPA1_2A.longestpeptide.fa.gz | sed '/^>/ s/\(.*\)P/\1T/' > /mnt/netapp1/Store_CSIC/home/csic/eye/lmf/snpEff/data/LYPA1_2A/protein.fa
```

Now we can build the database with the script [snpeff_build_database.sh](https://github.com/luciamayorf/Variant_calling_and_filtering/blob/main/functional_annotation/snpeff_build_database.sh):
```{bash}
sbatch -t 00:05:00 --mem 10G /home/csic/eye/lmf/scripts/functional_annotation/snpeff/snpeff_build_database.sh     # job ID: 7509429
```

CDS percentage error: 0.27%. Protein percentage error: 4.64%.

## Running SnpEff

Dani runs SnpEff defining intervals in a BED file. This gives the variants intersecting those regions the annotation provided by the BED file. They have a very detailed annotation file, containing intergenic regions (and promoters, upstream/downstream, UTRs, etc.), but we don't have that information (only genes, transcripts, cds and exons).

```{bash}
java -Xmx16g -jar $EBROOTSNPEFF/snpEff.jar LYPA1_2A -v -c /mnt/netapp1/Store_CSIC/home/csic/eye/lmf/snpEff/snpEff.config -s "/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/functional_annotation/snpeff/snpeff_c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.html" -csvStats "/mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/functional_annotation/snpeff/snpeff_c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.csv" /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.vcf > /mnt/netapp2/Store_csebdjgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/functional_annotation/snpeff/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss_annotated.vcf.gz
```

