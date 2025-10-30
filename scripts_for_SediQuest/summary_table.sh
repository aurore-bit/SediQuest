#!/usr/bin/env bash

set -euo pipefail

splitbam="${13}"

if [ -e $splitbam ] ; then
  sc=$(samtools view -c $splitbam)
else
  sc="NA"
fi


map_bam="$1"
rmdup_bam="$2"
target_bam="$3"
deam_bam="$4"
split_kraken="$5"
split_kraken_deam="$6"
deam_stats="$7"
summary_unique="$8"
contam="$9"
count="${10}"
count_deam="${11}"
summary_annotated="${12}"


indexlibid="${14}"
score_n="${15}"
score_b="${16}"
probeset="${17}"


# Get read counts for each input BAM file
map_bam_reads=$(/home/bioinf/usr/bin/samtools view -c "$map_bam")
rmdup_bam_reads=$(/home/bioinf/usr/bin/samtools view -c "$rmdup_bam")
target_bam_reads=$(/home/bioinf/usr/bin/samtools view -c "$target_bam")
deam_bam_reads=$(/home/bioinf/usr/bin/samtools view -c "$deam_bam")
primates=$(/home/bioinf/usr/bin/samtools view -c "$split_kraken")
primates_deam=$(/home/bioinf/usr/bin/samtools view -c "$split_kraken_deam")

        # Get deam stats (5'CT_95CI, 3'CT_95CI, cond5'CT_95CI, cond3'CT_95CI) from deam_stats file
deam_5CT=$(awk 'NR==2 {{print $2}}' "$deam_stats")
deam_3CT=$(awk 'NR==2 {{print $3}}' "$deam_stats")
deam_5CT_95CI=$(awk 'NR==2 {{print $4}}' "$deam_stats")
deam_3CT_95CI=$(awk 'NR==2 {{print $5}}' "$deam_stats")
deam_cond5CT=$(awk 'NR==2 {{print $8}}' "$deam_stats")
deam_cond3CT=$(awk 'NR==2 {{print $9}}' "$deam_stats")
deam_cond5CT_95CI=$(awk 'NR==2 {{print $10}}' "$deam_stats")
deam_cond3CT_95CI=$(awk 'NR==2 {{print $11}}' "$deam_stats")
    
     # Get %unique and %exhausted from exhaustion file
average_dup=$(awk 'NR==2 {{print $11}}' "$summary_unique")

        #get contamination estimate
contamination=$(awk 'NR==8 {{print $2}}' "$contam")
err_estimate=$(awk 'NR==8 {{print $3}}' "$contam")

tmp="$summary_annotated"

        # Write header to the summary file
echo IndexLibID N_score MD_score probeset split mapped unique target kraken_target deam kraken_deam 5'CT 3'CT 5'CT_95CI 3'CT_95CI cond5'CT cond3'CT cond5'CT_95CI cond3'CT_95CI average_dup SNPs_count_target SNPs_count_deam Contamination Contamination_err_estimate | tr ' ' '\t' > $tmp

        # Append the values for each field to the summary file
echo "$indexlibid" "$score_n" "$score_b" "$probeset" \
$sc \
$map_bam_reads \
$rmdup_bam_reads \
$target_bam_reads \
$primates \
$deam_bam_reads \
$primates_deam \
"$deam_5CT" \
"$deam_3CT" \
"$deam_5CT_95CI" \
"$deam_3CT_95CI" \
"$deam_cond5CT" \
"$deam_cond3CT" \
"$deam_cond5CT_95CI" \
"$deam_cond3CT_95CI" \
"$average_dup" \
$(cat "$count") \
$(cat "$count_deam") \
"$contamination" \
"$err_estimate" \
tr ' ' '\t' >> $tmp
