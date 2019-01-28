#!/bin/bash

# ref: https://www.biostars.org/p/112251/

if [ $# -eq 0 ]
  then
    echo "${0##*/} [chrom_size.tsv] [gene_gff] (both file are in the current directory)"
    exit
fi


chrom_size=$1
gene_gff=$2
b_name=$(basename ${gene_gff%.*}) # or ${0##*/}

# module load bedtools

mkdir sub_beds
cd sub_beds

# sort gff
cat ../$gene_gff | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' > "${b_name}.sorted.gff3"

# filter gene.gff
cat "${b_name}.sorted.gff3" | awk '$1 ~ /^#/ {print $0;next}$3=="gene"{print}' > "${b_name}.gene.gff3"

# intergenic gff (genome - gene)
bedtools complement -i "${b_name}.gene.gff3" -g ../$chrom_size > "${b_name}.intergenic.bed"

# generate exon.gff (gff - gene <cds, utr, mrna>); bed is 0-base, gff is 1-base
cat "${b_name}.sorted.gff3" | awk 'BEGIN{OFS="\t"} NF==9 && $3=="CDS"{print $1,$4-1,$5}' > "${b_name}.exon.bed"

# generate intron.gff (genome - intergenic - exon)
bedtools complement -i <(cat "${b_name}.intergenic.bed" "${b_name}.exon.bed"| sort -k1,1 -k2,2n) -g ../$chrom_size > "${b_name}.intron.bed"

