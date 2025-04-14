#!/bin/bash

# Read the CSV file line by line, skipping the header
tail -n +2 to_correct_flag.csv | while IFS=$'\t' read -r library probeset; do
  bam="${library}_uniqL35MQ25_MDALL_NALL_deam.bam"
  output="${library}_uniqL35MQ25_MDALL_NALL_deam_clean.bam"

  # Construct the subdirectory path dynamically for each combination
  subdir="mappedbams/${library}/${probeset}/rmdupL35MQ25/target/Mam_div_score_ALL/N_score_ALL/deam/"

  # Debug: Print the directory path to check if it's correct
  echo "Checking directory: $subdir"

  # Check if the directory exists
  if test -d "$subdir"; then
    echo "Directory exists: $subdir"
    
    # Change into the directory
    cd "$subdir" || exit 1

    # Run the samtools command to process the BAM file and adjust flag
    samtools view -h "$bam" | \
      awk '{if (!($0 ~ /^@/) && and($2,0x200)==512) {printf $1"\t"$2-512; for (i=3; i<=NF; i++) printf "\t"$i; printf "\n"} else {print $0}}' | \
      samtools view -hb > "$output"

    # Remove specific files
    rm "${library}_uniqL35MQ25_MDALL_NALL_deam_depth"
    rm "${library}_uniqL35MQ25_MDALL_NALL_deam_average"
    rm "${library}_uniqL35MQ25_MDALL_NALL_deam.cov_burden.txt"
    rm "${library}_uniqL35MQ25_MDALL_NALL_deam.bam_average"

  else
    echo "Directory not found for ${library} and ${probeset} in ${subdir}. Skipping..."
  fi
done
