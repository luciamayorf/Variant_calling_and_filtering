# Variant_calling_and_filtering

In this repository I will perform the variant calling and filtering of a high coverage dataset of WGS samples of 50 individuals at ~30X.


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

We later perform a small quality check with [bcftools stats](https://samtools.github.io/bcftools/bcftools.html#stats) to make sure everything went correctly and there are not weird samples. It doesn't take long to perform the operation, so I can just run it interactively in a loop.

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

---

## Variant filtering

### 0. Identification of repetitive regions

Before filtering, we need to identify the repetitive regions in our reference genome. 

REDACTAR

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


### 2. Mean read depth filtering

After testing how to filter each individual mean read depth, we saw that most windows were shared across individuals, so we decided to apply this filter for the whole population.

#### Calculating the mean read depth in 10kbp windows

First we will calculate mean read depth in consecutive 10kbp windows along the genome using [mosdepth](https://github.com/brentp/mosdepth), which generates a bed file with the mean coverage per 10kbp window.

I run the script [mean_rd_filtering_bed_generation.sh](https://github.com/luciamayorf/Variant_calling_and_filtering/blob/main/scripts/depth_filtering_bed_generation.sh) <input_bam> <output_path>:
```bash
for input_bam in $(ls /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_bams/novogene_lp_sept23/*_sorted_rg_merged_sorted_rmdup_indelrealigner.bam); do 
  job_id=$(sbatch -c 16 --mem=10GB -t 00:15:00 /home/csic/eye/lmf/scripts/variant_filtering/depth_filtering_bed_generation.sh ${input_bam} /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/depth_filtering | awk '{print $4}')
    echo "${job_id} ${input_bam}" >> /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/mosdepth/job_ids_depth_filtering_bed_generation.txt
done
```

CAUTION: from now on, we will also **exclude the scaffolds**, as they usually present very high read depth values (and very few defined 10kbp windows, as they are quite small). If I ever need to recover this information, I would need to apply a different depth coverage filter for these scaffolds.

#### Obtaining the bed of excluded windows

To obtain the bed file with the excluded windows and a plot of the mean read depth distribution, I run the script [mean_rd_filter_population_bed.R](https://github.com/luciamayorf/Variant_calling_and_filtering/blob/main/scripts/mean_rd_filter_population_bed.R) <input_directory> <samples_list> <output_directory>

```bash
Rscript /home/csic/eye/lmf/scripts/variant_filtering/mean_rd_filter_population_bed.R /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/depth_filtering/ <(cut -f1,2 -d'_' /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/FASTQ_files/novogene_lp_sept2023/fastq_samples_list.txt | sort -u) /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/depth_filtering/
```

#### Applying the filter

We decide to apply the mean+0.5*sd as the read-depth cutoff for the population analysis, using the script [variant_filter_rd.sh](https://github.com/luciamayorf/Variant_calling_and_filtering/blob/main/scripts/variant_filter_rd.sh). 1182 windows will be excluded (out of a total of 242578 windows: ~0.49% of windows excluded).
```bash
sbatch -c 5 --mem=5GB -t 00:10:00 /home/csic/eye/lmf/scripts/variant_filtering/variant_filter_rd.sh /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20.vcf /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/depth_filtering/excluded_windows_mean_rd_pop_filter.bed 
```
Now we have the following number of SNPs:
| Filter       | N of variants  | Filtered variants |
|:-------------|:--------------:|:-----------------:|
| No filter    |     4288231    |   0               |
| Filter1      |     2062793    |   2225438         |
| Filter2      |     1816919    |   245874          |
| Filter3      |     1810275    |   6644            |
| Filter4      |     1721036    |   89239           |
| Filter5      |     1354440    |   366596          |
| Mean_rd      |     1327386    |   27054           |

### 3. Missing genotypes filtering

Firstly, we will generate a plot with the % of SNPs that would be included when the proportion of missing data allowed is increased, using the script [missingness_plot.R](https://github.com/luciamayorf/Variant_calling_and_filtering/blob/main/scripts/missingness_plot.R) <input_data> <output_directory>
```{bash}
Rscript /home/csic/eye/lmf/scripts/variant_filtering/missingness_plot.R /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/missingness_filtering/missing_gt_count_c_lp_novogene_sept23.txt /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/mLynPar1.2_ref_vcfs/novogene_lp_sept23/missingness_filtering/
```
[cumulative_miss_plot.pdf](https://github.com/luciamayorf/Variant_calling_and_filtering/files/14422118/cumulative_miss_plot.pdf)


