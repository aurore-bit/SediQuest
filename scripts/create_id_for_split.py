import pandas as pd
import os

# Read the CSV file using pandas
indexlibid_df = pd.read_csv("../bam_to_split.csv", sep=";")

# Group the DataFrame by 'indexlibid' and 'probeset'
grouped = indexlibid_df.groupby(['indexlibid','new_id', 'probeset'])

# Iterate over each group
for (indexlibid, new_id, probeset), group in grouped:
    # Define the directory path
    dir_path = os.path.join(str(new_id), probeset)
    os.makedirs(dir_path, exist_ok=True)
    # Define the file path
    file_path = os.path.join(dir_path, 'id.txt')

    # Open the file for writing
    with open(file_path, 'w') as file:
        # Write the header
        file.write('#index\tp7\tp5\n')

        # Write the data rows
        for _, row in group.iterrows():
            file.write("{0}\t{1}\t{2}\n".format(row['indexlibid'], row['p7'], row['p5']))
