#!/bin/bash

while IFS= read -r line; do
    # Check if the third column matches the word "gene"
    if [[ $(echo "$line" | awk '{print $3}') == "gene" ]]; then
        echo "$line" >> snpeff_annotation/LYPA1_2A.FA_longestisoform_CESGA.gff3
        GENE=$(echo "$line" | cut -f1 -d';' | cut -f9 | sed 's/ID=//')
        TRANSCRIPT=$(grep "${GENE}" snpeff_annotation/LYPA1_2A.transcripts_longest_CESGA.fa | sed 's/>//')
        grep -w "${TRANSCRIPT}" LYPA1_2A.FA.gff3 >> snpeff_annotation/LYPA1_2A.FA_longestisoform_CESGA.gff3
    else
        echo "Skypping this line, not a gene"
    fi
done < LYPA1_2A.FA.gff3

