# Variant_calling_and_filtering

In this repository I will perform the variant calling and filtering of a high coverage dataset of WGS samples of 50 individuals at ~30X. All the analyses were performed at CESGA server.


## Variant calling

We will use the WGS model from [DeepVariant](https://github.com/google/deepvariant) to call variants in our dataset.

For the pseudoatosomal regions, we will establish a standard PAR1 region of 7 Mb in the beginning and at the end of the X chromosome, according to the existing bibliography: between [6 Mb](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5155386/)-[6.5Mb](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5522595/) and [10 Mb](https://onlinelibrary.wiley.com/doi/full/10.1111/eva.13302).

We will use the script [variant_calling_deepvariant.sh](https://github.com/luciamayorf/Variant_calling_and_filtering/blob/main/scripts/variant_calling_deepvariant.sh) <{input_bam}> <ref_genome> <files_list> <output_directory>:

```bash
for input_bam in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/*_sorted_rg_merged_sorted_rmdup_indelrealigner.bam); do 
  job_id=$(sbatch -c 32 --mem=50GB -t 06:00:00 --gres=gpu:a100:1 /home/csic/eye/lmf/scripts/deepvariant/variant_calling_deepvariant.sh ${input_bam} /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.scaffolds.revcomp.scaffolds.fa /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/FASTQ_files/novogene_lp_sept2023/fastq_samples_list.txt /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_gvcfs/novogene_lp_sept2023 | awk '{print $4}')
    echo "${job_id} ${input_bam}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/deepvariant/job_ids_deepvariant_novogene_lp_sept2023.txt
done
```

We later perform a small quality check with [bcftools stats](https://samtools.github.io/bcftools/bcftools.html#stats) to make sure that everything went correctly and that there are not weird samples. It doesn't take long to perform the operation, so we can just run it interactively in a loop.

#### Variant QC

```bash
module load samtools

for i in /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_gvcfs/novogene_lp_sept2023/*_mLynPar1.2_ref.vcf.gz; do
  bcftools stats ${i} > /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_gvcfs/novogene_lp_sept2023/bcftools_stats/$(basename "${i}" .vcf.gz)_stats.txt
done

module load multiqc 
multiqc *_stats.txt
```

---

## gVCFs merging

We will use [GLnexus](https://github.com/dnanexus-rnd/GLnexus), which converts multiple VCF files to a single BCF file, which then needs to be converted to a VCF. We will use the configuration "DeepVariantWGS", which already applies some soft quality filters (AQ >10) to decrease the false positive rate and the VCF file size.

First, I need to create a list of all the gVCFS I want to merge:
```bash
ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_gvcfs/novogene_lp_sept2023/*g.vcf.gz > /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_gvcfs/novogene_lp_sept2023/novogene_lp_sep2023_gvcfs_list.txt
```

CAUTION: GLnexus creates a directory "GLnexus.DB" in the folder where the job is launched. If that folder already exists, the program won't run.

Now we can run the script [glnexus_script.sh](https://github.com/luciamayorf/Variant_calling_and_filtering/blob/main/scripts/glnexus_script.sh) <gvcfs_list> <output_file>: 
```bash
sbatch -c 30 --mem=100G -t 03:00:00 /home/csic/eye/lmf/scripts/glnexus_script.sh /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_gvcfs/novogene_lp_sept2023/novogene_lp_sep2023_gvcfs_list.txt /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/c_lp_all_novogene_sept23_mLynPar1.2_ref_originalnames.vcf.gz
```

I need to take an additional step to change the name of the samples to include their origin population:
```bash
module load samtools

bcftools reheader -s <(cut -f2 /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/FASTQ_files/novogene_lp_sept2023/fastq_samples_list.txt | uniq) -o c_lp_all_novogene_sept23_mLynPar1.2_ref.vcf.gz c_lp_all_novogene_sept23_mLynPar1.2_ref_originalnames.vcf.gz
```


#### VCF QC

We also perform a small QC with [bcftools stats](https://samtools.github.io/bcftools/bcftools.html#stats) to make sure that the gVCFs were merged correctly.
```bash
module load samtools
module load multiqc

bcftools index -t c_lp_all_novogene_sept23_mLynPar1.2_ref.vcf.gz
bcftools stats c_lp_all_novogene_sept23_mLynPar1.2_ref.vcf.gz > c_lp_all_novogene_sept23_mLynPar1.2_ref_stats.txt

multiqc c_lp_all_novogene_sept23_mLynPar1.2_ref_stats.txt
```

---

## Variant filtering

### 0. Identification of repetitive regions

Before filtering, we need to identify the repetitive regions in our reference genome. 

I ran [RepeatModeler](https://www.repeatmasker.org/RepeatModeler/) to identify transposable elements, but that was previously done during the reference genome assembly by CNAG collaborators. These repetitive regions are found in the file /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynPar1.2/Repeats.4jb.gff3.gz

After testing different options, we decided to mask the genome with a combination of the Repeats.4jb.gff3.gz file of intersperse repeats and the low complexity regions calculated by [Repeatmasker](https://www.repeatmasker.org/) using the script [repeatmasker_lowcomplex.sh](https://github.com/luciamayorf/Variant_calling_and_filtering/blob/main/scripts/repeatmasker_lowcomplex.sh)

```bash
sbatch /home/csic/eye/lmf/scripts/repetitive_regions/repeatmasker_lowcomplex.sh
```

#### Checking percentage of the genome masked

```{bash}
# generating the merged GFF
zcat Repeats.4jb.gff3.gz > repeats_temp.gff
cat repeats_temp.gff repetitive_regions/low_complex/mLynPar1.2.scaffolds.revcomp.scaffolds.fa.out.gff | sort -k1,1 -k4,4n > repeats_lowcomplexity_mLynPar1.2.scaffolds.revcomp.gff3
rm repeats_temp.gff

# generating the merged BED:
awk -F'\t' '{OFS="\t"; print $1, $4-1, $5}' repeats_lowcomplexity_mLynPar1.2.scaffolds.revcomp.gff3 > repeats_lowcomplexity_mLynPar1.2.scaffolds.revcomp.bed

# merging the repetitive from both files in the BED:
bedtools merge -i repeats_lowcomplexity_mLynPar1.2.scaffolds.revcomp.bed > repeats_lowcomplexity_mLynPar1.2.scaffolds.revcomp_merged.bed

# Calculate the total length of intersected regions:
awk '{sum+=$3-$2} END {print sum}' repeats_lowcomplexity_mLynPar1.2.scaffolds.revcomp_merged.bed     # = 1042291408

# Calculating the % of the genome that would be masked (taking into account that the total mLynPar1.2genome size is: 2465344228 <-- mLynPar1.2.scaffolds.revcomp.scaffolds.fa.fai)
1042291408/2465344228
```
42,28% of the genome is masked with this library. 

IMPORTANT: the BED file containing the windows with repetitive regions to filter out is /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynPar1.2/**repeats_lowcomplexity_mLynPar1.2.scaffolds.revcomp_merged.bed**


---

### 1. First round filtering (filters 1 to 5)

We are going to apply the following filters to the VCF file:

1.  Filter 1: removing variants in **repetitive/low complexity** regions.
2.  Filter 2: removing **non-biallic** sites.
3.  Filter 3: removing **invariant**sites (they would be singletons of the reference genome).
4.  Filter 4: removing variants with a **low quality score**. I think QUAL >= 20 should be enough (99% confident that the genotype is real)
5.  Filter 5: removing **indels**.

We apply the filters by running the script [variant_filter_1to5.sh](https://github.com/luciamayorf/Variant_calling_and_filtering/blob/main/scripts/variant_filter_1to5.sh) <ref.fa> <in.vcf> <masked_regions.bed>
```bash
sbatch -t 00:30:00 -c 5 --mem 10GB /home/csic/eye/lmf/scripts/variant_filtering/variant_filter_1to5.sh /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.scaffolds.revcomp.scaffolds.fa /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/c_lp_all_novogene_sept23_mLynPar1.2_ref.vcf.gz /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynPar1.2/repeats_lowcomplexity_mLynPar1.2.scaffolds.revcomp_merged.bed
```

These are the number of SNPs removed in each of the filters:
| Filter       | N of variants  | Filtered variants |
|:-------------|:--------------:|:-----------------:|
| No filter    |     4288231    |   0               |
| Filter1      |     2062793    |   2225438         |
| Filter2      |     1816919    |   245874          |
| Filter3      |     1810275    |   6644            |
| Filter4      |     1663124    |   147151          |
| Filter5      |     1325874    |   337250          |


To obtain the VCFs with the QUAL20 filter (I will probably use this VCF from now on), I run the script [variant_filter_1to5_QUAL20.sh]() <ref.fa> <in.vcf>:
```{bash}
sbatch -t 00:10:00 -c 5 --mem 10GB /home/csic/eye/lmf/scripts/variant_filtering/variant_filter_1to5_QUAL20.sh /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.scaffolds.revcomp.scaffolds.fa /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/c_lp_all_novogene_sept23_mLynPar1.2_ref.vcf.gz #job ID: 5990353
```

| Filter       | N of variants  | Filtered variants |
|:-------------|:--------------:|:-----------------:|
| No filter    |     4288231    |   0               |
| Filter1      |     2062793    |   2225438         |
| Filter2      |     1816919    |   245874          |
| Filter3      |     1810275    |   6644            |
| Filter4      |     1721036    |   89239           |
| Filter5      |     1354440    |   366596          |

---

### 2. Mean read depth filtering

After testing how to filter each individual per chromosome mean read depth (to exclude possible paralogs found in each individual), we saw that most windows were shared across individuals, so we decided to apply this filter for the whole population.

#### Calculating the mean read depth in 10kbp windows

First we will calculate mean read depth in consecutive 10kbp windows along the genome using [mosdepth](https://github.com/brentp/mosdepth), which generates a bed file with the mean coverage per 10kbp window.

I run the script [mean_rd_filtering_bed_generation.sh](https://github.com/luciamayorf/Variant_calling_and_filtering/blob/main/scripts/depth_filtering_bed_generation.sh) <input_bam> <output_path>:
```bash
for input_bam in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/*_sorted_rg_merged_sorted_rmdup_indelrealigner.bam); do 
  job_id=$(sbatch -c 16 --mem=10GB -t 00:15:00 /home/csic/eye/lmf/scripts/variant_filtering/depth_filtering_bed_generation.sh ${input_bam} /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/depth_filtering | awk '{print $4}')
    echo "${job_id} ${input_bam}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/mosdepth/job_ids_depth_filtering_bed_generation.txt
done
```

IMPORTANT: from now on, we will also **exclude the scaffolds**, as they usually present very high read depth values (and very few defined 10kbp windows, as they are quite small). If I ever need to recover this information, I would need to apply a different depth coverage filter for these scaffolds (I have not removed the SNPs yet).


#### Obtaining the bed of excluded windows

To obtain the bed file with the excluded windows and a plot of the mean read depth distribution, I run the script [mean_rd_filter_population_bed.R](https://github.com/luciamayorf/Variant_calling_and_filtering/blob/main/scripts/mean_rd_filter_population_bed.R) <input_directory> <samples_list> <output_directory>

```bash
Rscript /home/csic/eye/lmf/scripts/variant_filtering/mean_rd_filter_population_bed.R /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/depth_filtering/ <(cut -f1,2 -d'_' /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/FASTQ_files/novogene_lp_sept2023/fastq_samples_list.txt | sort -u) /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/depth_filtering/
```

The output of the script is a BED file with the excluded windows, and a mean read depth distribution plot with the different cutoff values.

![mean_rd_cutoffs_pop](https://github.com/luciamayorf/Variant_calling_and_filtering/assets/96131721/88d68b5d-6995-4200-96f8-583e662465f1)

IMPORTANT: the BED file containing the windows with a higher mean read depth to filter out is /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/depth_filtering/**excluded_windows_mean_rd_pop_filter.bed** 

1182 windows will be excluded (out of a total of 242578 windows: ~0.49% of windows excluded), with the mean+0.5*sd cutoff.

#### Applying the filter

We decide to apply the mean+0.5*sd as the read-depth cutoff for the population analysis, using the script [variant_filter_rd.sh](https://github.com/luciamayorf/Variant_calling_and_filtering/blob/main/scripts/variant_filter_rd.sh). 

```bash
sbatch -c 5 --mem=5GB -t 00:10:00 /home/csic/eye/lmf/scripts/variant_filtering/variant_filter_rd.sh /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20.vcf /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/depth_filtering/excluded_windows_mean_rd_pop_filter.bed 
```

After this the read depth filter, we have the following number of SNPs:

| Filter       | N of variants  | Filtered variants |
|:-------------|:--------------:|:-----------------:|
| No filter    |     4288231    |   0               |
| Filter1      |     2062793    |   2225438         |
| Filter2      |     1816919    |   245874          |
| Filter3      |     1810275    |   6644            |
| Filter4      |     1721036    |   89239           |
| Filter5      |     1354440    |   366596          |
| Mean_rd      |     1327386    |   27054           |

---

#### Genotypes QC

Before applying the missingness, I will apply a fast genotypes QC with vcftools and [plink](https://www.cog-genomics.org/plink/)

I will calculate the allele frequencies, the mean depth per individual and the mean depth per site, the site quality, the proportion of missing data per individual and the number of singletons to check the quality of the genotypes and the heterozigosity and inbreeding coefficients.

```bash
sbatch -c 5 --mem 5GB -t 00:15:00 /home/csic/eye/lmf/scripts/genotypes_QC/genotypes_qc_vcftools.sh /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.vcf /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/vcf_stats # job ID: 5995880
```

To generate the corresponding plots, we use the script [genotypes_qc_plots.R](https://github.com/luciamayorf/Variant_calling_and_filtering/blob/main/scripts/genotypes_qc_plots.R) <input_vcf_prefix> </path/to/tables/>
```bash
Rscript /home/csic/eye/lmf/scripts/genotypes_QC/genotypes_qc_plots.R c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/vcf_stats/
```

---

### 3. Missing genotypes filtering

We will apply a genotype missingness filtering to avoid analyzing variants that are not informative enough and might introduce noise into our results. First, we will calculate the number of missing genotypes for each SNP in order to draw a distribution of data missingness across the entire genome.

```{bash}
bcftools query -f '%CHROM:%POS\t[%GT\t]\n' c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.vcf | awk '{missing=0; for(i=2; i<=NF; i++) if($i=="./.") missing++; print $1, missing}' > ./missingness_filtering/missing_gt_count_c_lp_novogene_sept23.txt
```

Then, we will generate a plot with the % of SNPs that would be included when the proportion of missing data allowed is increased, using the script [missingness_plot.R](https://github.com/luciamayorf/Variant_calling_and_filtering/blob/main/scripts/missingness_plot.R) <input_data> <output_directory>
```{bash}
Rscript /home/csic/eye/lmf/scripts/variant_filtering/missingness_plot.R /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/missingness_filtering/missing_gt_count_c_lp_novogene_sept23.txt /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/missingness_filtering/
```

![cumulative_miss_plot](https://github.com/luciamayorf/Variant_calling_and_filtering/assets/96131721/31ba6f0c-178c-41e3-b071-2cd5e29e68b1)

#### Applying the filter

For now, we decide to apply a very lax filter: we will filter out SNPs missing in more than 70% of the individuals. Later, we can apply a more strict filter if needed. I use the script [variant_filter_miss.sh](https://github.com/luciamayorf/Variant_calling_and_filtering/blob/main/scripts/variant_filter_miss.sh)

```bash
sbatch -c 5 --mem=5GB -t 00:10:00 /home/csic/eye/lmf/scripts/variant_filtering/variant_filter_miss.sh /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.vcf
```

| Filter       | N of variants  | Filtered variants |
|:-------------|:--------------:|:-----------------:|
| No filter    |     4288231    |   0               |
| Filter1      |     2062793    |   2225438         |
| Filter2      |     1816919    |   245874          |
| Filter3      |     1810275    |   6644            |
| Filter4      |     1721036    |   89239           |
| Filter5      |     1354440    |   366596          |
| Mean_rd      |     1327386    |   27054           |
| Miss         |     1327251    |   135             |
