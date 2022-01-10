#!/bin/bash

# updated 9/20/2020

########################################
###	     Adjust Input Variables	     ###
########################################

# path to repeat masker track (with "chr" deleted from chromosome names)
repmasker="/path/to/hg19RepMasker_removed_chr.gtf"

# path(s) to the directory with newly generated gtf files - should have cds_regions, utr3/5_annotation
declare -a arr=(
	"/path/to/gtf_files/"
)


# generate cds GTF, utr3 GTF and utr5 GTF (remove overlaps and repeats)

for direc in "${arr[@]}"
do
  cd "$direc" || { echo "cd Failure"; exit 1; }
  bedtools subtract -a cds_regions.gtf -b $repmasker > cds_regions_cleaned.gtf
  bedtools subtract -a utr3_annotation.gtf -b cds_regions.gtf > utr3_annotation_removed_cds.gtf #removing residual CDS
  bedtools subtract -a utr3_annotation_removed_cds.gtf -b $repmasker > utr3_annotation_cleaned.gtf #getting rid of SINE/LINE
  bedtools subtract -a utr5_annotation.gtf -b cds_regions.gtf > utr5_annotation_removed_cds.gtf
  bedtools subtract -a utr5_annotation_removed_cds.gtf -b $repmasker > utr5_annotation_cleaned.gtf
done
