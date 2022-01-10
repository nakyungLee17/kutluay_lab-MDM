#!/bin/bash

# REQUIRES plot_readlength_dist_mdm.R


# directory to mapping output folder; this folder contains several subfolders,
# each containing mapping outputs of one single sample (ex. one subfolder for each barcode).
# the file used to calculate cellular read length distribution is a transcriptome alignment BAM file.
declare -a arr=(
  "/path/to/mapping_output_folder/" # <<< keep slash at the end
)

## cellular transcriptome

for direc in "${arr[@]}" 
do

  cd "$direc" || { echo "cd Failure"; exit 1; }

  mkdir readlen_distributions
  outdir=$direc'/readlen_distributions'

  ## find transcriptome bam file and count the read length distributions
  for i in $(find . -name '*Trans*.bam' -exec readlink -f {} \;)
  do
    IFS='/' read -r -a array <<< "$(dirname $i)"
    echo ">>>>>>> processing: ${array[10]}"    # <<< this array index is used to distinguish each sample
    samtools view -F 4 "$i" | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"'| sort | uniq -c | tee "$outdir/${array[10]}_cellular.txt" # <<< change array index accordingly
  done


  # plot distributions
  # make sure the R script is in correct path

  cd readlen_distributions || { echo "cd Failure"; exit 1; }

  for i in $(find . -name '*.txt' -exec readlink -f {} \;)
  do
    IFS='.' read -r -a array <<< "$(basename $i)"
    Rscript /path/to/plot_readlength_dist_mdm.R "$i" "$outdir" "${array[0]}"
  done

done
