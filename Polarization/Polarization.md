
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
I move all the BAM files to EBD genomics server.

## Base synteny file

The next step consists of generating a fasta file for each species with the Iberian lynx reference genome, but containing the alleles for each species in the variants positions of the pardinus VCF. For that, we use [samtools-mpileup](https://www.htslib.org/doc/samtools-mpileup.html) and [put2fa](https://github.com/Paleogenomics/Chrom-Compare/tree/master).

```bash
# Define general paths
REF=/GRUPOS/grupolince/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.scaffolds.revcomp.scaffolds.fa
CHR=($(cat /GRUPOS/grupolince/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.scaffolds.revcomp.scaffolds.fa.fai | cut -f 1 | uniq))
VCF=/GRUPOS/grupolince/mLynPar1.2_ref_vcfs/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.vcf.gz

############################ RUFUS ############################
# Define paths
INBAM=/GRUPOS/grupolince/mLynPar1.2_ref_bams/polarization/lr1_mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner.bam
OUTFASTA=/GRUPOS/grupolince/lucia/polarization/fasta/lr1_mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner_pu2fa.fa
INTERSECTFASTA=/GRUPOS/grupolince/lucia/polarization/lr1_mLynPar1.2_ref_intersect.fa
INTERSECTBED=/GRUPOS/grupolince/mLynPar1.2_ref_bams/polarization/lr1_mLynPar1.2_ref_intersect.bed

# get fasta from bam per chr
rm $OUTFASTAls

for i in ${CHR[@]:0:20}; do
  echo "getting FASTA for chromosome $i"
  samtools mpileup -s -q30 -f $REF $INBAM -r $i | /GRUPOS/grupolince/Chrom-Compare/pu2fa -c $i -C 100 >> $OUTFASTA
done

echo "getting intersection FASTAs"

bedtools getfasta -fi $OUTFASTA -bed $VCF -fo $INTERSECTFASTA # we only get 16257 N (Lorena had 222076)
bedtools getfasta -fi $OUTFASTA -bed $VCF -fo $INTERSECTBED -tab


############################ MARBLED CAT ############################
# Define paths
INBAM=/GRUPOS/grupolince/mLynPar1.2_ref_bams/polarization/pm_mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner.bam
OUTFASTA=/GRUPOS/grupolince/lucia/polarization/fasta/pm_mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner_pu2fa.fa
INTERSECTFASTA=/GRUPOS/grupolince/lucia/polarization/pm_mLynPar1.2_ref_intersect.fa
INTERSECTBED=/GRUPOS/grupolince/lucia/polarization/pm_mLynPar1.2_ref_intersect.bed

# get fasta from bam per chr
rm $OUTFASTA

for i in ${CHR[@]:0:20}; do
  echo "$i"
  samtools mpileup -s -q30 -f $REF $INBAM -r $i | /GRUPOS/grupolince/Chrom-Compare/pu2fa -c $i -C 100 >> $OUTFASTA
done

bedtools getfasta -fi $OUTFASTA -bed $VCF -fo $INTERSECTFASTA    # we only get 30736 N (Lorena had 222076)
bedtools getfasta -fi $OUTFASTA -bed $VCF -fo $INTERSECTBED -tab


############################ CAT ############################
# Define paths
INBAM=/GRUPOS/grupolince/mLynPar1.2_ref_bams/polarization/fc_mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner.bam
OUTFASTA=/GRUPOS/grupolince/lucia/polarization/fasta/fc_mLynPar1.2_ref_sorted_rg_merged_sorted_rmdup_indelrealigner_pu2fa.fa
INTERSECTFASTA=/GRUPOS/grupolince/lucia/polarization/fc_mLynPar1.2_ref_intersect.fa
INTERSECTBED=/GRUPOS/grupolince/lucia/polarization/fc_mLynPar1.2_ref_intersect.bed

# get fasta from bam per chr
rm $OUTFASTA

for i in ${CHR[@]:0:20}; do
  echo "$i"
  samtools mpileup -s -q30 -f $REF $INBAM -r $i | /GRUPOS/grupolince/Chrom-Compare/pu2fa -c $i -C 100 >> $OUTFASTA
done

bedtools getfasta -fi $OUTFASTA -bed $VCF -fo $INTERSECTFASTA # we only get 26685 N
bedtools getfasta -fi $OUTFASTA -bed $VCF -fo $INTERSECTBED -tab
```

Then, we generate the synteny BED file, which contains the nucleotide present in each SNP for 3 ougroup species and the reference allele for the L. pardinus.

```bash
# There's a problem with the cat intersect file. The last 13 SNPs are missing in the C1 chromosome. I have checked and haven't notice anything weird in that chromosome (it doesn't have more "Ns" than the others), there are simply no reads mapping to that region. I add them (with an "N") manually:
sed -i '226132r fc_missing_SNPs_ChrC1.txt' fc_mLynPar1.2_ref_intersect.bed

# generating cat-marbled_cat intersection
cut -f2 pm_mLynPar1.2_ref_intersect.bed | paste fc_mLynPar1.2_ref_intersect.bed - > fc_pm_sinteny_intersect.bed

# generating cat-marbled_cat-rufus-pardinus synteny file (removing scaffolds and ChrY_unloc from the VCF, as they are not in the fasta files)
cut -f2 lr1_mLynPar1.2_ref_intersect.bed | paste fc_pm_sinteny_intersect.bed - > temp_fc_pm_lr.bed
zgrep -v "#" /GRUPOS/grupolince/mLynPar1.2_ref_vcfs/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.vcf.gz | grep -v "scaffold\|unloc" | cut -f4 | paste -d'\0' temp_fc_pm_lr.bed - | awk -F"\t|:|-" '{printf ("%s\t%s\t%s\t%s=%s%s%s\n", $1,$2,$3,"fc_pm_lr_lp",$4,$5,$6)}' > fc_pm_lr_lp_sinteny_intersect.bed
rm temp_fc_pm_lr.bed
```

Apply the parsimonia criteria to add the ancestral column
```bash
awk '{                                       
split($0,a,":");
split(a[1],b,"=");
split(b[2],c,"");
if (c[1]==c[2] && c[1]==c[3]) printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,c[3],c[4]); # all equal
else if (c[1]==c[2] && c[1]!=c[3]) printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,"N",c[4]); #fc pm vs lr N
else if (c[1]==c[3] && c[1]!=c[2]) printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,c[3],c[4]); #fc lr* vs pm
else if (c[2]==c[3] && c[1]!=c[3]) printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,c[3],c[4]); #pm lr * vs fc

else if ((c[1]=="N") && c[2]==c[3]) printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,c[3],c[4]); # fc lr * all equal
else if ((c[1]=="N") && c[2]!=c[3]) printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,"N",c[4]); # lr vs fc N

else printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,"N",c[4]); #others null
}' fc_pm_lr_lp_sinteny_intersect.bed > outgroup_parsimony_ancestral_state_fc_pm_lr_lp_variants.bed
```

The 5th column in this file represents the ancestral allele, and the 6th column contains the allele set as reference in the pardinus VCF.

Now, we will generate different files to polarize the vcf and generate the ancestral fasta reference genome.

```bash
# Generate file with unpolarizable sites (can't determine the ancestral state)
awk '$5=="N" {print $0}' outgroup_parsimony_ancestral_state_fc_pm_lr_lp_variants.bed >  unpolarizable_fc_pm_lr_lp_variants.bed

    ##### 70561/1284590 = 0.0549 --> 5,49% of the sites are not polarizable
    
# Generate file with polarizable sites 
awk '$5!="N" {print $0}' outgroup_parsimony_ancestral_state_fc_pm_lr_lp_variants.bed > polarizable_ancestral_state_fc_pm_lr_lp_variants.bed

    ##### 1214029/1284590: 94,51% of the sites are polarizable. 1214029 total SNPS.

# Generate file with inconsistent sites (wrongly polarised - ancestral state needs to be changed - or unpolarizable) --> for the ancestral reference fasta generation
awk '$5!=$6 {print $0}' outgroup_parsimony_ancestral_state_fc_pm_lr_lp_variants.bed > inconsistent_ancestral_state_fc_pm_lr_lp_variants.bed         
    
    ##### 549422 SNPs, of which 70561 are not polarizable, so the remaning 478861 SNPs with the reference allele will be changed (39,44% of all polarizable SNPs)

# Create the tsv file to annotate the VCF and index it
awk 'BEGIN {OFS="\t"} {print $1, $3, $5}' polarizable_ancestral_state_fc_pm_lr_lp_variants.bed | bgzip -c > aa_annotation_polarizable_snps_fc_pm_lr_lp.tsv.gz
tabix -s1 -b2 -e2 aa_annotation_polarizable_snps_fc_pm_lr_lp.tsv.gz            # -b2 and -e2 to indicate that it is 1-based
```

## VCF polarization

For the VCF polarization, we will first add a column in the INFO field containing the ancestral state (AA).
```bash
cd /GRUPOS/grupolince/mLynPar1.2_ref_vcfs

# I generate a VCF file excluding the SNPs in scaffolds and ChrY_unloc (I should've done this from the beginning)
zgrep -v "_scaffold_\|_unloc_" c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.vcf.gz > c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss_bigChr.vcf
bgzip c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss_bigChr.vcf

# First, we are going to keep only the sites in my VCF that are polarizable.
bedtools intersect -a c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss_bigChr.vcf.gz -b /GRUPOS/grupolince/lucia/polarization/polarizable_ancestral_state_fc_pm_lr_lp_variants.bed -header > c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss_bigChr_polarizable.vcf
bgzip c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss_bigChr_polarizable.vcf

# Then we add the ancestral allele column to the INFO field
/opt/bcftools-1.6/bcftools annotate -a /GRUPOS/grupolince/lucia/polarization/aa_annotation_polarizable_snps_fc_pm_lr_lp.tsv.gz -c CHROM,POS,INFO/AA c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss_bigChr_polarizable.vcf -h <(echo '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">') -Oz -o c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss_bigChr_polarizable_aafilled.vcf.gz
```

To polarize the SNPs in the VCF (alleles will be switched whenever the ancestral allele matches the alternative one, and genotypes will be properly recoded as well, based on the new INFO/AA column). I take the code from [Dani's improved_polarization.Rmd, line 4372](https://github.com/Arynio/genome-annotation/blob/master/improved_polarisation.Rmd)

```bash
# The following code was originally provided by Pierre Lindenbaum and modified by José Luis Castro.
java -jar /opt/jvarkit/dist/vcffilterjdk.jar -e 'if(variant.getNAlleles()!=2 || !variant.hasAttribute("AA")) return true; 
final String aa = variant.getAttributeAsString("AA",""); 
if(!variant.getAlleles().get(1).getDisplayString().equalsIgnoreCase(aa)) return true; 
VariantContextBuilder vb=new VariantContextBuilder(variant); 

Allele oldalt = variant.getAlleles().get(1);
Allele oldref = variant.getAlleles().get(0); 
Allele ref= Allele.create(oldalt.getDisplayString(),true); 
Allele alt= Allele.create(oldref.getDisplayString(),false);

vb.alleles(Arrays.asList(ref,alt)); 

List genotypes= new ArrayList<>(); 
for(Genotype g: variant.getGenotypes()) 
  { 
  if(!g.isCalled()) 
  { genotypes.add(g); continue;} 
  GenotypeBuilder gb = new GenotypeBuilder(g); 
  List alleles = new ArrayList<>(); 
  for(Allele a:g.getAlleles()) { 
    if(a.equals(oldalt)) { a=ref;} 
    else if(a.equals(oldref)) { a=alt;} 
    alleles.add(a); 
    } 
  if(g.hasPL()) { 
    int pl[] = g.getPL(); 
    int pl2[] = new int[pl.length]; 
    for(int i=0;i< pl.length;i++) pl2[i]=pl[(pl.length-1)-i]; 
    gb.PL(pl2); 
    } 
  if(g.hasAD()) 
    { int ad[] = g.getAD(); 
    int ad2[] = new int[ad.length]; 
    for(int i=0;i< ad.length;i++) ad2[i]=ad[(ad.length-1)-i];
    gb.AD(ad2); 
  } 
  genotypes.add(gb.alleles(alleles).make()); 
  }

vb.attribute("AF",1.0d - Double.parseDouble(variant.getAttributeAsString("AF",""))); vb.attribute("AC",variant.getGenotypes().stream().flatMap(G->G.getAlleles().stream()).filter(A->A.equals(oldref)).count()); 
vb.genotypes(genotypes); 
return vb.make();' -o c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss_bigChr_polarizable_aafilled_polarized.vcf.gz c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss_bigChr_polarizable_aafilled.vcf.gz
```

### Detect ancestral state allele incompatibilities

Later on, I noticed that there were some SNPs whose neither of the two alleles (reference or alternative) matched the ancentral state, and they were not changed during the VCF polarization.

```bashcut -f2 lr1_mLynPar1.2_ref_intersect.bed | paste fc_pm_sinteny_intersect.bed - > temp_fc_pm_lr.bed
zgrep -v "#" /GRUPOS/grupolince/mLynPar1.2_ref_vcfs/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss.vcf.gz | grep -v "scaffold\|unloc" | cut -f4,5 | sed 's/\t//' | paste -d'\0' temp_fc_pm_lr.bed - | awk -F"\t|:|-" '{printf ("%s\t%s\t%s\t%s=%s%s%s\n", $1,$2,$3,"fc_pm_lr_lpRefAlt",$4,$5,$6)}' > temp_incompatibilities.bed

# Adding the ancestral state column to remove non-polarizable variants
awk '{                                       
split($0,a,":");
split(a[1],b,"=");
split(b[2],c,"");
if (c[1]==c[2] && c[1]==c[3]) printf ("%s\t%s\t%s\t%s %s\n", $1,$2,$3,$4,c[3]); # all equal
else if (c[1]==c[2] && c[1]!=c[3]) printf ("%s\t%s\t%s\t%s %s\n", $1,$2,$3,$4,"N"); #fc pm vs lr N
else if (c[1]==c[3] && c[1]!=c[2]) printf ("%s\t%s\t%s\t%s %s\n", $1,$2,$3,$4,c[3]); #fc lr* vs pm
else if (c[2]==c[3] && c[1]!=c[3]) printf ("%s\t%s\t%s\t%s %s\n", $1,$2,$3,$4,c[3]); #pm lr * vs fc

else if ((c[1]=="N") && c[2]==c[3]) printf ("%s\t%s\t%s\t%s %s\n", $1,$2,$3,$4,c[3]); # fc lr * all equal
else if ((c[1]=="N") && c[2]!=c[3]) printf ("%s\t%s\t%s\t%s %s\n", $1,$2,$3,$4,"N"); # lr vs fc N

else printf ("%s\t%s\t%s\t%s %s\n", $1,$2,$3,$4,"N"); #others null
}' temp_incompatibilities.bed | awk -F ' ' '$5!="N" {print $0}' > temp_incompatibilities_polarizable.bed


# Compare both of the lp alleles with the ancestral state to see if there's any incompatibility
sed 's/ //g' temp_incompatibilities_polarizable.bed |  
awk '{                                       
split($0,a,":");
split(a[1],b,"=");
split(b[2],c,"");
if (c[4]==c[6]) printf ("%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,"ok_ref"); # ref allele matches any of the other species alleles
else if (c[5]==c[6]) printf ("%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,"ok_alt"); # alt allele matches any of the other species alleles
else printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$4,"incompatible","N"); # none of the alleles matches any of the other species alleles
}' | grep "incompatible" > allele_incompatibilities_polarization.bed
rm temp*.bed
```
There are 3317 SNPs whose none of two alleles match the ancestral state. Those variants are already included in the inconsistencies file, but they have the wrong ancestral state (should be an N, as that allele is not found in the reference). To change the ancestral state for an N in that file:

```bash
awk 'NR==FNR{a[$1,$2,$3]=$0;next}{if(($1,$2,$3) in a){print $1"\t"$2"\t"$3"\t"$4"\t""N""\t"$6}else{print $0}}' allele_incompatibities_polarization.bed inconsistent_ancestral_state_fc_pm_lr_lp_variants.bed > ancestral_state_changes_fc_pm_lr_lp_variants.bed
```

I also need to remove those variants from the polarized VCF
```bash
bedtools intersect -a /GRUPOS/grupolince/mLynPar1.2_ref_vcfs/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss_bigChr_polarizable_aafilled_polarized.vcf.gz -b allele_incompatibities_polarization.bed -v -header > /GRUPOS/grupolince/mLynPar1.2_ref_vcfs/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss_bigChr_polarized.vcf
bgzip /GRUPOS/grupolince/mLynPar1.2_ref_vcfs/c_lp_all_novogene_sept23_mLynPar1.2_ref.filter5_QUAL20_rd.miss_bigChr_polarized.vcf
```

## Ancestral reference fasta generation

[GATK FastaAlternateReferenceMaker](https://gatk.broadinstitute.org/hc/en-us/articles/360037594571-FastaAlternateReferenceMaker) takes a VCF and the reference genome as inputs and replaces the reference variant for the alternative one for the positions in the VCF.

For that, I first need to generate a "pseudovcf" (it will only have useful information about the alleles and the SNP position) file from the BED file containing the positions that need to be changed (the alternative allele column containing the ancestral state).

```bash
# Generate the VCF
awk 'BEGIN {printf("##fileformat=VCFv4.2\n##INFO=<ID=END,Number=1,Type=Integer>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");} {printf("%s\t%s\t.\t%s\t%s\t.\t.\tEND=%d\n",$1,int($2)+1,$6,$5,$3);}' ancestral_state_changes_fc_pm_lr_lp_variants.bed > ancestral_state_changes_fc_pm_lr_lp_variants.vcf

# Generate the ancestral fasta
java -jar /opt/GATK/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker \
-R /GRUPOS/grupolince/reference_genomes/lynx_pardinus_mLynPar1.2/mLynPar1.2.scaffolds.revcomp.scaffolds.fa \
-V ancestral_state_changes_fc_pm_lr_lp_variants.vcf \
-o ancestral_genome/ancestral_mLynPar1.2.scaffolds.revcomp.scaffolds.fa

# Check that everything went well
rm kaka.borrar
seqkit fx2tab ancestral_genome/ancestral_mLynPar1.2.scaffolds.revcomp.scaffolds.fa > ancestral_genome/ancestral_mLynPar1.2.scaffolds.revcomp.scaffolds.tab
while read -r SCAFFOLD START STOP SYNTENY ANCESTRAL PARDINUS; do
  OLD=$(grep "$SCAFFOLD" mLynPar1.2.scaffolds.revcomp.scaffolds.tab | awk '{printf ("%s\n", $2)}' | cut -c$STOP)
  NEW=$(grep "$SCAFFOLD" ancestral_genome/ancestral_mLynPar1.2.scaffolds.revcomp.scaffolds.tab | awk '{printf ("%s\n", $3)}' | cut -c$STOP)
  echo -e "$SCAFFOLD\t$STOP\t$PARDINUS\t$OLD\t$ANCESTRAL\t$NEW" >> kaka.borrar
done < ancestral_state_changes_fc_pm_lr_lp_variants.bed
```
