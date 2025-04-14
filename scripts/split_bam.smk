import pandas as pd
import os

# Read the CSV file using pandas
indexlibid_df = pd.read_csv("bam_to_split.csv", sep=";")

# Get the unique 'indexlibid' values
#INDEXLIBID = indexlibid_df['indexlibid'].unique()

# Create a dictionary to map indexlibid to probeset
INDEXLIBID_TO_PROBESETS = indexlibid_df.groupby('new_id')['probeset'].apply(list).to_dict()

INDEXLIBID_PROBESET_PAIRS = list(indexlibid_df[['new_id', 'probeset']].itertuples(index=False, name=None))

INDEXLIBID_OLD_PAIRS = list(indexlibid_df[['new_id', 'indexlibid']].itertuples(index=False, name=None))

# Retrieve BAM file paths based on 'indexlibid'
def get_bam(wildcards):
    return indexlibid_df.loc[(indexlibid_df['new_id'] == wildcards.indexlibid) & 
                             (indexlibid_df['probeset'] == wildcards.probeset), "bam"].iloc[0]


rule all:
    input:
        expand(["/mnt/expressions/Aurore/sediment_pipeline_v0/overview_output/split/{indexlibid}/{probeset}/{old_indexlibid}.bam"],
               zip, indexlibid=[pair[0] for pair in INDEXLIBID_PROBESET_PAIRS],
                   probeset=[pair[1] for pair in INDEXLIBID_PROBESET_PAIRS],
 old_indexlibid=[pair[1] for pair in INDEXLIBID_OLD_PAIRS])

#split bam 
rule split_bam:
    input:
        ogbam = get_bam,
        ids = "/mnt/expressions/Aurore/sediment_pipeline_v0/overview_output/split/{indexlibid}/{probeset}/id.txt"
    output:
        old_split = "/mnt/expressions/Aurore/sediment_pipeline_v0/overview_output/split/{indexlibid}/{probeset}/{old_indexlibid}.bam",
    params: 
        splitbams = "{indexlibid}.bam",
        old="{old_indexlibid}.bam"
    shell:
        """
        echo 'Splitting the BAM file'
        outdir=$(dirname {output.old_split})
        cd $outdir
        time /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/splitBAM.pl -byfile {input.ids} {input.ogbam} > split_stats.txt
        """
