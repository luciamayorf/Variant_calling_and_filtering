#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/functional_annotation/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/functional_annotation/slurm-%j.out

output_dir="/mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_pardinus_mLynPar1.2/snpeff_annotation/new_version_Paulina_june2025"
ref_dir="/mnt/netapp2/Store_csebdjgl/reference_genomes/lynx_pardinus_mLynPar1.2"

while IFS= read -r line; do
    # Check if the third column matches the word "gene"
    if [[ $(echo "$line" | awk '{print $3}') == "gene" ]]; then
        echo "$line" >> ${output_dir}/LYPA1_2A.longest.isoform_single_transcript.gff3
        GENE=$(echo "$line" | cut -f9 | awk -F'[=;]' '{print $2}')
	    PRODUCT=$(grep "${GENE}" ${output_dir}/LYPA1_2A.longest.cds.fa | sed 's/>//' | cut -f1 -d' ')
        TRANSCRIPT=($(grep -w "${PRODUCT}" ${ref_dir}/LYPA1_2A.FA.gff3 | cut -f3 -d';' | sed 's/Name=//'))
        grep -w "${TRANSCRIPT[0]}" ${ref_dir}/LYPA1_2A.FA.gff3 >> ${output_dir}/LYPA1_2A.longest.isoform_single_transcript.gff3
        
    fi
done < ${ref_dir}/LYPA1_2A.FA.gff3
