################################################################################
# Draft nuclear sediment pipeline
#
# Aurore GALTIER, 25/11/2024

#FOR NOW YOU NEED TO:
#-have a config folder with config, probeset and sample files
#-copy split bam (raw file) in /output_v0/split/{indexlibid}/{probeset}/
#-copy mapped bam in /output_v0/mappedbams/{indexlibid}/{probeset}/
################################################################################

from snakemake.utils import R
import pandas as pd
import glob

#where to run the script
workdir: "/mnt/expressions/Aurore/sediment_pipeline_v0/"

#to have errors output
if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")

PREFIX=config["prefix"]

# Define project folder
project = config["project"]

#define b_score
score_b = config["score_b"]

# Read the indexlibid information
indexlibid_df = pd.read_table(config["samples_info"], comment="#").set_index(["indexlibid"], drop=False).sort_index()
INDEXLIBID = indexlibid_df.index.unique()

# Retrieve BAM file paths from the indexlibid dataframe
def get_bam(wildcards):
    return indexlibid_df.loc[wildcards.indexlibid, "bam"]

# Read probeset information
probeset_df = pd.read_table(config["probeset_info"], comment="#").set_index(["probeset"], drop=False).sort_index()

# Retrieve reference paths from the probeset dataframe
def get_ref(wildcards):
    return probeset_df.loc[wildcards.probeset, "path_to_ref"]

# Retrieve reference bed file path from the probeset dataframe
def get_bed(wildcards):
    return probeset_df.loc[wildcards.probeset, "path_to_bed"]

def get_control(wildcards):
    return probeset_df.loc[wildcards.probeset, "path_to_control"]

def get_control_info(wildcards):
    return probeset_df.loc[wildcards.probeset, "path_to_control_info"]


##############################################
# map_all Rule
##############################################
sorted_bams = expand('{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.sorted.bam',
    indexlibid=INDEXLIBID,
    probeset=indexlibid_df.loc[INDEXLIBID, "probeset"],
    project=project)

rule map_all:
    input:
        bam_files = sorted_bams,
    # Use the list of files here
    run:
        print("Hey, mapping is done")
        pass

#rule map:
 #   input:
  #      bam = get_bam,
   # output:
    #    mapped = "{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.bam",
     #   temp = temp("{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.bam.tmp")
   # threads: 8
   # params: ref = get_ref
    #shell:
     #   """
      #  echo 'Hey! Mapping BAM files...'
       # bam-fixpair -o {output.temp} {input.bam}
        #bwa bam2bam -t {threads} -g {params.ref} -n 0.01 -o 2 -l 16500 --only-aligned -f {output.mapped} {output.temp}
       # """

rule samtools:
    input:
        mapped = "{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.bam"
    output:
        sorted = "{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.sorted.bam"
    #temp = temp("{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.bam.tmp")
    # sorted = "{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.sorted.bam"
    threads: 8
    shell:
        """
        /home/visagie/.local/bin/samtools sort --threads {threads} -@30 -o {output.sorted} {input.mapped}
        """



##############################################
#processing
##############################################

rule process_all:
    input: 
        expand("{project}/mappedbams/{indexlibid}/{probeset}/target/{score_b}_filter/rmdupL35MQ25/deam/{indexlibid}.bam",
            indexlibid=INDEXLIBID,
            probeset=indexlibid_df.loc[INDEXLIBID, "probeset"],
            project=project,
            score_b=score_b)
    run:
        print('hey! Running target filtering, deamination, quality filtering and duplicate removal :)')
        pass

#keeping the . for now as a 0 should we?
rule filter_control_sites:
    input:
        control_sites = get_bed,
        control= get_control,
    output:
        sites_filtered="{project}/mappedbams/{indexlibid}/{probeset}/target/{score_b}_filter/sites_{b_score}.filtered.txt",
    params:
        b_score=score_b,
    threads: 1
    #conda: "envs/processing.yaml"
    shell: """
        Rscript scripts/filter_SNPs.R {score_b} {input.control_sites} {input.control} {output.sites_filtered}
    """

rule filter_bam_by_control_sites:
    input:
        sites_filtered="{project}/mappedbams/{indexlibid}/{probeset}/target/{score_b}_filter/sites_{score_b}.filtered.txt",
        bam = "{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.bam",
    output:
        sites_bam="{project}/mappedbams/{indexlibid}/{probeset}/target/{score_b}_filter/{indexlibid}.bam",
    threads: 1
    #conda: "envs/processing.yaml"
    shell: """
        bedtools intersect -a {input.bam} -b {input.sites_filtered} > {output.sites_bam} 
    """

rule uniq_q25_l35:
    input:
        bam="{project}/mappedbams/{indexlibid}/{probeset}/target/{score_b}_filter/{indexlibid}.bam"
    output:
        bam="{project}/mappedbams/{indexlibid}/{probeset}/target/{score_b}_filter/rmdupL35MQ25/{indexlibid}.bam",
        bai="{project}/mappedbams/{indexlibid}/{probeset}/target/{score_b}_filter/rmdupL35MQ25/{indexlibid}.bai",
        exhaust = "{project}/mappedbams/{indexlibid}/{probeset}/target/{score_b}_filter/rmdupL35MQ25/{indexlibid}.bam.exhaustion"
    params:
        unique_bam = "{indexlibid}.uniq.L35MQ25.bam",
        unique_bai= "{indexlibid}.uniq.L35MQ25.bam.bai",
        new_bam="{indexlibid}.bam",
        new_bai="{indexlibid}.bai",
        unique_sum="summary_stats_L35MQ25.txt"
    shell: """
    /home/mmeyer/perlscripts/solexa/analysis/analyzeBAM.pl -nof -minlength 35 -qual 25 {input.bam} > {output.exhaust}
    mv {params.unique_bam} {params.new_bam}
    mv {params.unique_bai} {params.new_bai}
    mv {params.new_bam}  $(dirname {output.bam})
    mv {params.new_bai}  $(dirname {output.bam})
    mv {params.unique_sum}  $(dirname {output.bam})
           """
 # bam-rmdup -o {output.bam} -q 25 -l 35 {input.bam}

rule deam_filter:
    input: bam="{project}/mappedbams/{indexlibid}/{probeset}/target/{score_b}_filter/rmdupL35MQ25/{indexlibid}.bam"
    output:
            bam="{project}/mappedbams/{indexlibid}/{probeset}/target/{score_b}_filter/rmdupL35MQ25/deam/{indexlibid}.bam",
            stats="{project}/mappedbams/{indexlibid}/{probeset}/target/{score_b}_filter/rmdupL35MQ25/deam/{indexlibid}.deam53x3.stats"
    #conda: "envs/processing.yaml"
    params:
        unique_bam="{indexlibid}.deam53x3.bam",
        new_bam="{indexlibid}.bam"
    shell: """
    /home/mmeyer/perlscripts/solexa/analysis/filterBAM.pl  -p5 0,1,2 -p3 0,-1,-2 -suffix deam53x3 {input.bam} &> {output.stats}
    mv {params.unique_bam}  {params.new_bam}
    mv {params.new_bam} $(dirname {output.bam})
    """



##############################################
#kraken
##############################################

rule kraken:
    input:
        expand(["{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/deam/{indexlibid}.kraken_spc",
            "{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/deam/{indexlibid}.byread",
            "{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/{indexlibid}.kraken_spc",
            "{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/{indexlibid}.byread",
            "{project}/mappedbams/{indexlibid}/{probeset}/target/{indexlibid}.kraken_spc",
            "{project}/mappedbams/{indexlibid}/{probeset}/target/{indexlibid}.byread",
            "{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.kraken_spc",
            "{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.byread",
            "{project}/split/{indexlibid}/{probeset}/{indexlibid}.kraken_spc",
            "{project}/split/{indexlibid}/{probeset}/{indexlibid}.byread"],
            indexlibid=INDEXLIBID,
            probeset=indexlibid_df.loc[INDEXLIBID, "probeset"],
            project=project)
    run:
        print('hey! Running Kraken')
        pass


rule bam_to_fasta:
    input:
        bam="{indexlibid}.bam"
    output:
        fa="{indexlibid}.fa.gz"
    threads: 1
   # conda: "envs/kraken.yaml"
    shell: """
     bam2fastx -a -A -Q {input.bam} | /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/lenfilter.pl 35 | gzip -c > {output.fa}
     """

rule generic_kraken:
    input: fa="{indexlibid}.fa.gz"
    output: kraken="{indexlibid}.kraken"
    threads: 1
   # conda: "envs/kraken.yaml"
    shell: """

     if [ $(hostname) != "bionc13" ] ; then 
         echo KRAKEN MUST BE RUN ON BIONC13. SLEEPING 30 SECONDS.
         sleep 30
         stop
     fi

    echo THIS FAILS IF THE KRAKEN DB IS NOT COPIED - fix it

     nproc=$(time ~frederic_romagne/kraken/install/kraken --threads {threads} --db /mnt/ramdisk/refseqReleaseKraken \
     --output {output.kraken} {input.fa} 2>&1 | grep 'processed' | cut -f1 -d' ')
     echo "count seqs {input.fa}"
     nseq=$(gunzip -c {input.fa} | grep -c -e '>' -e '@') || echo "no seqs found $nseq"

     if [ ! $nseq -eq $nproc ] ; then echo KRAKEN RUN FAILED; stop; fi
    """

rule generic_kraken_summary:
    input: kraken="{indexlibid}.kraken"
    output: phylo="{indexlibid}.kraken_phylo"
    threads: 1
   # conda: "envs/kraken.yaml"
    shell: """

    if [ $(hostname) != "bionc13" ] ; then 
         echo KRAKEN MUST BE RUN ON BIONC13. SLEEPING 30 SECONDS.
         sleep 30
         stop
    fi

    echo 'making report..'
    ## translate file has info on each read (but not in an easy-to-digest way)
    ## phylo file has a summary for each level in the taxonomy
    python3 /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/kraken_report.py --db /mnt/ramdisk/refseqReleaseKraken {input.kraken} > {output.phylo}
    """

rule generic_kraken_spc_summary:
    input: spc_krak="{indexlibid}.kraken_phylo",
        spc_groups="/mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/data/major_species_groups.txt"
    output: spc_summary="{indexlibid}.kraken_spc"
  #  conda: "envs/kraken.yaml"
    shell: """

    grep -f {input.spc_groups} {input.spc_krak} | awk '{{print $6,$2}}' > {output.spc_summary}

    n=$(awk '{{s += $2}} END {{print s}}' {output.spc_summary})
    echo "$n : " {input.spc_krak}

    m=$(grep '  Mammalia$' {input.spc_krak} | awk '{{print $2-'$n'}}')
    echo "$n : $m : " {input.spc_krak}

    b=$((grep '  Bacteria$' {input.spc_krak} || echo 0 0) | awk '{{print $2}}')
    echo "$n : $m : $b" {input.spc_krak}

    s=$((grep '  Sauropsida$' {input.spc_krak} || echo 0 0) | awk '{{print $2}}')
    echo "$n : $m : $b : $s :" {input.spc_krak}

    r=$(grep '	root$' {input.spc_krak} | awk '{{print $2-'$n'-'$m'-'$b'-'$s'}}')
    echo "$n : $m : $b : $s : $r :" {input.spc_krak}

    u=$(grep '	unclassified$' {input.spc_krak} | awk '{{print $2}}')
    echo "$n : $m : $b : $s : $r : $u : " {input.spc_krak}

    echo Mammalia $m >> {output.spc_summary}
    echo Bacteria $b >> {output.spc_summary}
    echo Sauropsida $s >> {output.spc_summary}
    echo root $r >> {output.spc_summary}
    echo unclassified $u >> {output.spc_summary}

"""

rule generic_kraken_translate:
    input: kraken="{indexlibid}.kraken"
    output: translate="{indexlibid}.translate"
    threads: 1
    shell: """

    if [ $(hostname) != "bionc13" ] ; then 
         echo KRAKEN MUST BE RUN ON BIONC13. SLEEPING 30 SECONDS.
         sleep 30
         stop
    fi

    echo 'making report..'
    ## translate file has info on each read (but not in an easy-to-digest way)
    ## phylo file has a summary for each level in the taxonomy
    python3 /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/kraken_report.py --db /mnt/ramdisk/refseqReleaseKraken {input.kraken} --translate {output.translate} > /dev/null

    """

rule generic_kraken_byread:
    input: translate="{indexlibid}.translate"
    output: byread="{indexlibid}.byread"
    threads: 1
    shell: """

     echo 'parsing translate file..'
     python3 /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/parse_kraken_translate_file.py --translate-file {input.translate} > {output.byread}

    """

##############################################
#summary
##############################################
rule summaries:
    input:
        expand(["{project}_summary/{indexlibid}/{probeset}/_burden_primates_plot.pdf",
        "{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/deam/deam_stats/{indexlibid}.L30-1000_MQ0.plots.pdf",
        "{project}_summary/{indexlibid}/{probeset}/{indexlibid}.pipeline_summary.txt",
        "{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/deam/deam_stats/{indexlibid}.anno_damage.txt",
        "{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/{indexlibid}.read_summary.txt.gz"],
        indexlibid=INDEXLIBID,
        probeset=indexlibid_df.loc[INDEXLIBID, "probeset"],
        project=project)
    run:
        print('hey! Creating summaries')
        pass

#create coverage file
rule cov_bed_files:
    input: "{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/{indexlibid}.bam"
    output:
            cov="{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/{indexlibid}.cov"
    shell:"""
            samtools mpileup -B {input} > {output.cov}
            """

#just count number of reads
rule pipeline_summary:
    input:
        map_bam="{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.bam",
        target_bam="{project}/mappedbams/{indexlibid}/{probeset}/target/{indexlibid}.bam",
        rmdup_bam="{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/{indexlibid}.bam",
        rmdup_stats="{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/summary_stats_L35MQ25.txt"
    output:
        summary_annotated="{project}_summary/{indexlibid}/{probeset}/{indexlibid}.pipeline_summary.txt"
    shell:
        """
        splitbam="{wildcards.project}/split/{wildcards.probeset}/{wildcards.indexlibid}.bam"
        #ls $splitbam

        if [ -e $splitbam ] ; then
          sc=$(samtools view -c $splitbam)
        else
          sc="NA"
        fi

        #tmp=$(mktemp)
        tmp={output.summary_annotated}

        echo IndexLibID probeset ref split mapped \
        $(head -n1 {input.rmdup_stats} | cut -f2- | sed 's/&//g' | sed 's/%/pct_/g' | sed 's/raw/target/g') \
        primate | tr ' ' '\t' > $tmp

        echo {wildcards.indexlibid} {wildcards.probeset} \
        $(samtools view -c {input.map_bam}) \
        $(tail -n1 {input.rmdup_stats} | cut -f2-) \
        | tr ' ' '\t' >> $tmp
        """

#create a file indicating every info for each read
rule bam_read_summary:
     input: bam = "{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/{indexlibid}.bam",
     	    control=get_control_info
     output: summary ="{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/{indexlibid}.read_summary.txt.gz"
     shell: """
            ## bam_basic_stats_pysam2.py is much faster - it goes through the bam read by read instead of looking for every position in the control file. should produce identical output to bam_basic_stats_pysam.py
     	    time /home/benjamin_vernot/miniconda3/bin/python /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/bam_basic_stats_pysam2.py \
	    	 --control {input.control} \
		    --bam {input.bam} \
		    --tags lib --tags-fill {wildcards.indexlibid} \
		    --control-header /mnt/expressions/benjamin_vernot/soil_capture_2017/site_categories_for_capture/soil_probe_designs/probes_CONTROL_HEADER.txt \
		    | gzip -c > {output.summary}
            """



#deam summaries files
rule deam_plot:
    input: bam="{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/deam/{indexlibid}.bam",
    output: plots="{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/deam/deam_stats/{indexlibid}.L30-1000_MQ0.plots.pdf",
    params: bam="../{indexlibid}.bam",
    shell: """
            cd $(dirname {output.plots})
            /home/mmeyer/perlscripts/solexa/analysis/substitution_patterns.pl {params.bam} 
            """

rule deam_summaries:
    input: bam="{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/deam/{indexlibid}.bam",
    output: summary="{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/deam/deam_stats/{indexlibid}.summary_damage.txt"
    shell: """    
            /home/mmeyer/perlscripts/solexa/analysis/summarize_CT_frequencies.pl -screen > {output.summary}
           """

rule annotate_deam_summary:
    input: damage="{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/deam/deam_stats/{indexlibid}.summary_damage.txt"
    output: damage_annotated="{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/deam/deam_stats/{indexlibid}.anno_damage.txt",
    shell:
        """
        sed 1d {input.damage} | tr '-' '\t' | cut -f 3- | awk  'BEGIN {{OFS="\t"}} {{print "{wildcards.indexlibid}", "{wildcards.probeset}", $0}}' > {output.damage_annotated}
        """

#create necessary plots
rule plot:
    input: "{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/{indexlibid}.cov",
    output:
            plot="{project}_summary/{indexlibid}/{probeset}/coverage_plot.png",
    run:
        shell(f"""
			Rscript scripts/primates_coverage_kraken_SNPs_burden_plot.R
		""")
