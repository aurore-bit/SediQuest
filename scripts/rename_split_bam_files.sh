#!/bin/bash

# For old Ben run
# Read the input CSV or TSV file
#input_file="overview_output/ben_old_run.csv"

# Loop through each line in the input file
#while IFS=';' read -r aa_code lane slash old_path file_path old_id new_id; do
    # Clean the new_id to remove unwanted characters like \r and spaces
 #   new_id_cleaned=$(echo "$new_id" | tr -d '\r' | tr -d '\t' | tr -d ';' | xargs)

    # Check if the directory exists
  #  if [ -d "overview_output/split/$new_id_cleaned/$aa_code" ]; then
        # Change to the correct directory
   #     cd "overview_output/split/$new_id_cleaned/$aa_code" || exit

        # Check if the old file exists
    #    if [ -f "$old_id.bam" ]; then
            # Move the file
     #       mv "$old_id.bam" "$new_id_cleaned.bam"
      #      echo "Moved $old_id.bam to $new_id_cleaned.bam in $(pwd)"
       # else
        #    echo "Error: $old_id.bam not found in $(pwd)"
       # fi

        # Change back to the original directory
        #cd - || exit
    #else
     #   echo "Error: Directory overview_output/split/$new_id_cleaned/$aa_code does not exist"
   # fi

#done < "$input_file"

#For bam I had to split
# Read the input CSV or TSV file
input_file="overview_output/bam_to_split.csv"


# Loop through each line in the input file
tail -n +2 "$input_file" | while IFS=';' read -r aa_code indexlibid new_id bam p7 p5; do
    # Clean the new_id to remove unwanted characters like \r and spaces
    new_id_cleaned=$(echo "$new_id" | tr -d '\r' | tr -d '\t' | tr -d ';' | xargs)

    # Check if the directory exists
    if [ -d "overview_output/split/$new_id_cleaned/$aa_code" ]; then
        # Change to the correct directory
        cd "overview_output/split/$new_id_cleaned/$aa_code" || exit

        # Check if the old file exists
        if [ -f "$indexlibid.bam" ]; then
            # Move the file
            mv "$indexlibid.bam" "$new_id_cleaned.bam"
            echo "Moved $indexlibid.bam to $new_id_cleaned.bam in $(pwd)"
        else
            echo "Error: $indexlibid.bam not found in $(pwd)"
        fi

        # Change back to the original directory
        cd - || exit
    else
        echo "Error: Directory overview_output/split/$new_id_cleaned/$aa_code does not exist"
    fi

done < "$input_file"


