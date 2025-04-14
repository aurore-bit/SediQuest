# Define file and paths
sample_info=/mnt/expressions/Aurore/sediment_pipeline_test/test.txt
mapped_bam_path=/mnt/expressions/Aurore/sediment_pipeline_test/test_pigment_alba/mappedbams
split_bam_path=/mnt/expressions/Aurore/sediment_pipeline_test/test_pigment_alba/splitbams

tail -n +2 "$sample_info" | while read -r line; do
  indexlibid=$(echo "$line" | awk '{print $2}')
  lane=$(echo "$line" | awk '{print $11}')
  seqrun=$(echo "$line" | awk '{print $9}')
  probeset=$(echo "$line" | awk '{print $1}')

echo "Processing indexlibid: $indexlibid with lane: $lane, seqrun: $seqrun, and probeset: $probeset"


# Make folders in the mapped_bams directory
cd "/mnt/expressions/Aurore/sediment_pipeline_v0/output_v0/mappedbams/"
  mkdir -p "$indexlibid/$probeset"

# Make folders in the splitbams directory
cd "/mnt/expressions/Aurore/sediment_pipeline_v0/output_v0/split/"
mkdir -p "$indexlibid/$probeset"

# Copy mapped bam files

cp "$mapped_bam_path/$indexlibid/$probeset/$seqrun/proc_1/$lane/split_exact/map_twist_1240k_third/$indexlibid.bam" "/mnt/expressions/Aurore/sediment_pipeline_v0/output_v0/mappedbams/$indexlibid/$probeset/"

  # Copy split bam files
cp "$split_bam_path/$seqrun/proc_1/$lane/split_exact/$indexlibid.bam" "/mnt/expressions/Aurore/sediment_pipeline_v0/output_v0/split/$indexlibid/$probeset/"
done
