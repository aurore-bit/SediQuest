#!/bin/bash

#for old ben run
# Read the input CSV or TSV file
input_file="overview_output/new_run.csv"  # Replace with the actual file path

# Loop through each line in the input file
while IFS=';' read -r aa_code lane run old_id new_id; do
    new_id_cleaned=$(echo "$old_id" | tr -d '\r' | tr -d '\t')   
    capture_id=$(echo "$new_id" | tr -d '\r' | tr -d '\t')

    # Create directories based on AA code and Library ID
   mkdir -p "overview_output/split/$new_id_cleaned/$aa_code"


    # Copy files from the source to the destination folder
    # Modify the source path based on the file_path and other variables in your line
    cp "$run$lane$capture_id.bam" "overview_output/split/$new_id_cleaned/$aa_code/$new_id_cleaned.bam"

    # Optional: You can add additional logic if the file copy involves more specific subdirectories or handling
    echo "Copied split bam file"
done < "$input_file"


#id=("Lib.J.8599" "Lib.L.5206" "Lib.L.5225" "Lib.L.5349" "Lib.L.5384" "Lib.L.5239" "Lib.L.5241" "Lib.F.9795" "Lib.G.1893" "Lib.G.1929")

#probe=("TW1")

#for i in "${id[@]}"; do
 #   mkdir -p "$i/$probe"
#done

