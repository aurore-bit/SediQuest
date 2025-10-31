################################################################################
#
# Aurore GALTIER, SediQuest pipeline 
#For human nuclear capture in sediments
#
################################################################################

from snakemake.utils import R
import pandas as pd
import glob
import os

#to have errors output
if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")

PREFIX=config["prefix"]

# Define project folder
project = config["project"]

#define score_b
score_b = config["score_b"]

#define n_score
score_n = config["score_n"]

kraken_db = config["kraken_db"]

#define group to extract from kraken output
kraken_group = config["kraken_group"]

#define the name of the kraken group
kraken_group_name = config["kraken_group_name"]

# Read the indexlibid information
indexlibid_df = pd.read_csv(config["samples_info"], comment="#", sep=";")
INDEXLIBID = indexlibid_df.groupby('indexlibid')['probeset_to_lib'].apply(list).to_dict()


# Retrieve BAM file paths from the indexlibid dataframe
def get_bam(wildcards):
    return indexlibid_df.loc[wildcards.indexlibid, "bam"]

# Read probeset information
probeset_df = pd.read_csv(config["probeset_info"], comment="#", sep="\t")
PROBESET = probeset_df.set_index('probeset').to_dict()


# Retrieve reference paths from the probeset dataframe
def get_ref(wildcards):
    return probeset_df.loc[probeset_df['probeset'] == wildcards.probeset, "path_to_ref"].values[0]

# Retrieve reference bed file path from the probeset dataframe
def get_bed(wildcards):
    return probeset_df.loc[probeset_df['probeset'] == wildcards.probeset, "path_to_bed"]

def get_control(wildcards):
    return probeset_df.loc[probeset_df['probeset'] == wildcards.probeset, "path_to_control"]


##############################################
#A rule for everything
##############################################

wildcard_constraints:
    score_n="[^/]+"


Final_bam = expand([f"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{{score_b}}/N_score_{score_n}/deam/{indexlibid}.uniq.L35MQ25_MD{{score_b}}_N{score_n}.deam.bam"
         for indexlibid, probesets in INDEXLIBID.items()
         for probeset in probesets],
         score_b=score_b)

bam_filterON = expand([f"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/combined_filter/{indexlibid}.uniq.L35MQ25_MD{{score_b}}_N{score_n}_K{kraken_group_name}_filterON.bam"
         for indexlibid, probesets in INDEXLIBID.items()
         for probeset in probesets],
        score_b=score_b)



fa_bam_mapped_target = expand([f"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{{score_b}}/N_score_{score_n}/kraken/{indexlibid}.fa.gz"
                     for indexlibid, probesets in INDEXLIBID.items()
                     for probeset in probesets],
                     score_b=score_b)


fa_bam_mapped_target_deam = expand([f"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{{score_b}}/N_score_{score_n}/deam/kraken/{indexlibid}.fa.gz"
                         for indexlibid, probesets in INDEXLIBID.items()
                         for probeset in probesets],
                         score_b=score_b)


Kraken_bam_mapped_target = expand([f"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{{score_b}}/N_score_{score_n}/kraken/{indexlibid}.kraken_spc"
         for indexlibid, probesets in INDEXLIBID.items()
         for probeset in probesets],
         score_b=score_b)


Kraken_byread_bam_mapped_target = expand([f"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{{score_b}}/N_score_{score_n}/kraken/{indexlibid}.byread"
         for indexlibid, probesets in INDEXLIBID.items()
         for probeset in probesets],
         score_b=score_b)

Kraken_bam_mapped_target_deam = expand([f"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{{score_b}}/N_score_{score_n}/deam/kraken/{indexlibid}.kraken_spc"
         for indexlibid, probesets in INDEXLIBID.items()
         for probeset in probesets],
         score_b=score_b)

Kraken_byread_bam_mapped_target_deam = expand([f"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{{score_b}}/N_score_{score_n}/deam/kraken/{indexlibid}.byread"
         for indexlibid, probesets in INDEXLIBID.items()
         for probeset in probesets],
         score_b=score_b)

Kraken_split_bam = expand([f"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{{score_b}}/N_score_{score_n}/split_kraken/{indexlibid}.uniq.L35MQ25_MD{{score_b}}_N{score_n}_K{kraken_group_name}.bam"
         for indexlibid, probesets in INDEXLIBID.items()
         for probeset in probesets],
         score_b=score_b)

Kraken_split_deam_bam = expand([f"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{{score_b}}/N_score_{score_n}/deam/split_kraken/{indexlibid}.uniq.L35MQ25_MD{{score_b}}_N{score_n}.deam_K{kraken_group_name}.bam"
         for indexlibid, probesets in INDEXLIBID.items()
         for probeset in probesets],
         score_b=score_b)
 

summary_table= expand([f"{project}_summary/{indexlibid}/{probeset}/Mam_div_score_{{score_b}}/N_score_{score_n}/{indexlibid}.pipeline_summary.txt"
        for indexlibid, probesets in INDEXLIBID.items()
        for probeset in probesets],
        score_b=score_b)

faunal_contam = expand([f"{project}_summary/{indexlibid}/{probeset}/Mam_div_score_{{score_b}}/N_score_{score_n}/{indexlibid}_cov_MD.pdf"
        for indexlibid, probesets in INDEXLIBID.items()
        for probeset in probesets],
        score_b=score_b)

kraken_plot = expand([f"{project}_summary/{indexlibid}/{probeset}/Mam_div_score_{{score_b}}/N_score_{score_n}/{indexlibid}_kraken_order_byburden.pdf"
        for indexlibid, probesets in INDEXLIBID.items()
        for probeset in probesets],
        score_b=score_b)


rule run_pipeline:
    input:
        bam_files=Final_bam,
        bam_filterON=lambda wildcards: bam_filterON if config.get("filter") == "HIGH" else [],
        fa_bam_mapped_target=fa_bam_mapped_target,
        fa_bam_mapped_target_deam=fa_bam_mapped_target_deam,
        Kraken_bam_mapped_target=Kraken_bam_mapped_target,
        Kraken_bam_mapped_target_deam=Kraken_bam_mapped_target_deam,
        Kraken_byread_bam_mapped_target_deam=Kraken_byread_bam_mapped_target_deam,
        Kraken_byread_bam_mapped_target=Kraken_byread_bam_mapped_target,
        Kraken_split_deam_bam=Kraken_split_deam_bam,
        Kraken_split_bam=Kraken_split_bam,
        summary_table=summary_table,
        cov_plot=faunal_contam,
        kraken_plot=kraken_plot
    run:
        print('Hello! The pipeline is running')
        pass


##############################################
#processing
##############################################


rule uniq_q25_l35:
    input:
        bam="{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.bam"
    output:
        bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/{indexlibid}.uniq.L35MQ25.bam",
        bai="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/{indexlibid}.uniq.L35MQ25.bam.bai",
        summary="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/summary_stats.uniq.L35MQ25.txt"
    shell: """
    ancient_dna_cpp_tools/analyzeBAM -out_folder $(dirname {output.bam}) -min_len 35 -min_map_qual 25 -remove_dups {input.bam}
           """
 
#If you wants to filter on a specific burden score
rule filter_control_sites_b_score:
    input:
        control_sites = get_bed,
        control= get_control,
    output:
        sites_filtered="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/sites_MD_{score_b}_N_{score_n}.filtered.txt",
    params:
       b_score = "{score_b}",
        #n_score=score_n,
    threads: 1
   # conda: "envs/processing.yaml"
    shell: """
        Rscript scripts_for_SediQuest/filter_SNPs.R {params.b_score} {input.control_sites} {input.control} {output.sites_filtered} 
    """

#To filter sites accordingly to you b score configuration choice
def get_sites_filtered(wildcards):
    if wildcards.score_b == "ALL":
        return get_control(wildcards) 
    else:
        return f"{wildcards.project}/mappedbams/{wildcards.indexlibid}/{wildcards.probeset}/rmdupL35MQ25/target/Mam_div_score_{wildcards.score_b}/N_score_{wildcards.score_n}/sites_MD_{wildcards.score_b}_N_{wildcards.score_n}.filtered.txt"


rule filter_bam_by_control_sites:
    input:
        sites_filtered=get_sites_filtered,
        bam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/{indexlibid}.uniq.L35MQ25.bam",
    output:
        sites_bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.bam",
        bai="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.bai",
    threads: 1
    #conda: "envs/processing.yaml"
    shell: """
        bedtools intersect -a {input.bam} -b {input.sites_filtered} > {output.sites_bam} 
        samtools index {output.sites_bam} {output.bai}
    """

rule deam_filter:
    input: bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.bam"
    output:
            bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.deam.bam",  
    shell: """
    ancient_dna_cpp_tools/filterBAM -p5 0,1,2 -p3 0,-1,-2 -suffix deam -out_folder $(dirname {output.bam}) {input.bam} 
    """



rule read_id_filter:
    input:
        burden_filter = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.read_summary.txt.gz",
        kraken_info = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/kraken/{indexlibid}.byread",
        burden_all = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_ALL/N_score_ALL/{indexlibid}.uniq.L35MQ25_MDALL_NALL.read_summary.txt.gz",
        burden_filter_deam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.deam.read_summary.txt.gz",
        kraken_info_deam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/kraken/{indexlibid}.byread",
        burden_all_deam =  "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_ALL/N_score_ALL/deam/{indexlibid}.uniq.L35MQ25_MDALL_NALL.deam.read_summary.txt.gz"
    output: 
        read_to_filter = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/combined_filter/reads_to_filter_MD{score_b}_N{score_n}.txt",
        read_to_filter_deam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/combined_filter/deam/reads_to_filter_MD{score_b}_N{score_n}.txt"
    params:
        score_b=score_b
    shell: """
        Rscript scripts_for_SediQuest/high_filter.R {params.score_b} {input.burden_filter} {input.kraken_info} {input.burden_all} {input.burden_filter_deam} {input.kraken_info_deam} {input.burden_all_deam} {output.read_to_filter} {output.read_to_filter_deam} 
        """

rule extract_reads:
    input:
        read_to_filter = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/combined_filter/reads_to_filter_MD{score_b}_N{score_n}.txt",
        bam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/{indexlibid}.uniq.L35MQ25.bam",
        read_to_filter_deam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/combined_filter/deam/reads_to_filter_MD{score_b}_N{score_n}.txt",
    output:
        bam_filter = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/combined_filter/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}_K{kraken_group_name}_filterON.bam",
        bam_filter_deam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/combined_filter/deam/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}_K{kraken_group_name}.deam_filterON.bam"
    shell: """
        samtools view -N {input.read_to_filter} -b {input.bam} -o {output.bam_filter}
        samtools view -N {input.read_to_filter_deam} -b {input.bam} -o {output.bam_filter_deam}
        """

##############################################
#kraken step 1 create a fasta file
##############################################


rule bam_to_fasta_1:
    input:
        bam_mapped_target="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.bam",
    output:
        fa_bam_mapped_target="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/kraken/{indexlibid}.fa.gz",
    threads: 1
    shell: """
         bam2fastx -a -A -Q {input.bam_mapped_target} | scripts_for_SediQuest/lenfilter.pl 35 | gzip -c > {output.fa_bam_mapped_target}
       """


rule bam_to_fasta_2:
    input:
        bam_mapped_target_deam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.deam.bam",
    output:
        fa_bam_mapped_target_deam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/kraken/{indexlibid}.fa.gz",
    threads: 1
    shell: """
         bam2fastx -a -A -Q {input.bam_mapped_target_deam} | scripts_for_SediQuest/lenfilter.pl 35 | gzip -c > {output.fa_bam_mapped_target_deam}
       """




##############################################
#kraken step 2 run kraken
##############################################


rule generic_kraken:
    input: 
        fa="{indexlibid}.fa.gz"
    output: 
        kraken="{indexlibid}.kraken"
    log:
        "{indexlibid}.kraken.log"
    threads: 1
    conda: "envs/kraken.yaml"
    shell: """
    kraken --threads {threads} --db {config[kraken_db]} \
    --output {output.kraken} {input.fa} 2>> {log}
    """

rule generic_kraken_summary:
    input: kraken="{indexlibid}.kraken"
    output: phylo="{indexlibid}.kraken_phylo"
    threads: 1
   # conda: "envs/kraken.yaml"
    shell: """
    python3 scripts_for_SediQuest/kraken_report.py --db {config[kraken_db]} {input.kraken} > {output.phylo}
    """

rule generic_kraken_spc_summary:
    input: spc_krak="{indexlibid}.kraken_phylo",
        #spc_groups="/mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/data/major_species_groups.txt"
    output: spc_summary="{indexlibid}.kraken_spc"
  #  conda: "envs/kraken.yaml"
    shell: """
    scripts_for_SediQuest/spc_kraken_summary.sh {input} {output}
    """


rule generic_kraken_translate:
    input: kraken="{indexlibid}.kraken"
    output: translate="{indexlibid}.translate"
    threads: 1
    shell: """
    python3 scripts_for_SediQuest/kraken_report.py --db {config[kraken_db]} {input.kraken} --translate {output.translate} > /dev/null

    """

rule generic_kraken_byread:
    input: translate="{indexlibid}.translate"
    output: byread="{indexlibid}.byread"
    threads: 1
    shell: """

     echo 'parsing translate file..'
     python3 scripts_for_SediQuest/parse_kraken_translate_file.py --translate-file {input.translate} > {output.byread}

    """


##############################################
#split kraken after deam
##############################################

rule generic_kraken_extract_after_deam:
    input: bam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.deam.bam",
           kraken = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/kraken/{indexlibid}.kraken",
           translate = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/kraken/{indexlibid}.translate"
    output: bam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/split_kraken/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.deam_K{kraken_group_name}.bam"
    threads: 1
    shell: """
    ofile=$(dirname {output.bam})/$(basename {output.bam} .bam)

    python3 scripts_for_SediQuest/kraken_report.py --db {config[kraken_db]} {input.kraken} --extractFile {input.bam} \
       --clades {config[kraken_group]} --extract-out-base $ofile   > /dev/null

    """


##############################################
#split kraken before deam
##############################################

rule generic_kraken_extract_before_deam:
    input: bam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.bam",
           kraken = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/kraken/{indexlibid}.kraken",
           translate = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/kraken/{indexlibid}.translate"
    output: bam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/split_kraken/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}_K{kraken_group_name}.bam"
    threads: 1
    shell: """
    ofile=$(dirname {output.bam})/$(basename {output.bam} .bam)

    python3 scripts_for_SediQuest/kraken_report.py --db {config[kraken_db]} {input.kraken} --extractFile {input.bam} \
       --clades {config[kraken_group]} --extract-out-base $ofile  > /dev/null

    """




##############################################
#summary
##############################################

#create coverage file for on target reads
rule cov_bed_files:
    input: "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.bam"
    output:
            cov="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.cov",
            count="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.count"
    shell:"""
            samtools mpileup -B {input} > {output.cov}
            awk '$4 != 0' {output.cov} | wc -l > {output.count}
            """


#get coverage also for deaminated reads
rule cov_bed_files_deam:
    input: "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.deam.bam"
    output:
            cov="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.deam.cov",
            count="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.deam.count"
    shell:"""
            samtools mpileup -B {input} > {output.cov}
            awk '$4 != 0' {output.cov} | wc -l > {output.count}
            """


#give an estimate of the contamination
rule contamination_estimaTe:
    input: 
        bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.bam"
    output:
        bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.authentict"
    shell: """
            samtools view {input.bam} | AuthentiCT deam2cont -o {output.bam} -s 10000 -
            """



#new deam summaries files
rule deam_stats_table:
    input: bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.bam",
    output: summary_new="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam_stats/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.summary_damage.txt"
    shell: """
            scripts_for_SediQuest/quick_substitutions.pl {input.bam} > {output.summary_new}
            """


#create the summary table
rule pipeline_summary:
    input:
        split_bam="{project}/split/{indexlibid}/{probeset}/{indexlibid}.bam",
        map_bam="{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.bam",
        target_bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.bam",
        rmdup_bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/{indexlibid}.uniq.L35MQ25.bam",
        deam_bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.deam.bam",
        split_kraken="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/split_kraken/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}_KPrimates.bam",
        split_kraken_deam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/split_kraken/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.deam_KPrimates.bam",
        summary_unique="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/summary_stats.uniq.L35MQ25.txt",
        deam_stats="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam_stats/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.summary_damage.txt",
        count="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.count",
        count_deam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.deam.count",
        contam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.authentict"
    output:
        summary_annotated="{project}_summary/{indexlibid}/{probeset}/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.pipeline_summary.txt",
    shell:
        """
        bash scripts_for_SediQuest/summary_table.sh {input.map_bam} {input.rmdup_bam} {input.target_bam} {input.deam_bam} {input.split_kraken} {input.split_kraken_deam} {input.deam_stats} {input.summary_unique} {input.contam} {input.count} {input.count_deam} {output.summary_annotated} {input.split_bam}  {wildcards.indexlibid} {wildcards.score_n} {wildcards.score_b} {wildcards.probeset} 
        """


#create a file indicating every info for each read
rule bam_read_summary_target:
     input: bam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.bam",
            control=get_control
     output: 
        summary ="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.read_summary.txt.gz",
        tmp_bed = temp("{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.tmp.bed")
     shell: """
            awk '{{print $2, $1, $3}}' OFS='\t' {input.control} > {output.tmp_bed}


            ## bam_basic_stats_pysam2.py is much faster - it goes through the bam read by read instead of looking for every position in the control file. should produce identical output to bam_basic_stats_pysam.py
     	    time python scripts_for_SediQuest/bam_basic_stats_pysam2.py \
	    	 --control {output.tmp_bed} \
		    --bam {input.bam} \
		    --tags lib --tags-fill {wildcards.indexlibid} \
		    --control-header scripts_for_SediQuest/probes_CONTROL_HEADER.txt \
		    | gzip -c > {output.summary}
            """

rule bam_read_summary_deam:
     input: bam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.deam.bam",
            control=get_control
     output: 
        summary ="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.deam.read_summary.txt.gz",
        tmp_bed = temp("{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}.tmp.bed")
     shell: """
            awk '{{print $2, $1, $3}}' OFS='\t' {input.control} > {output.tmp_bed}

            ## bam_basic_stats_pysam2.py is much faster - it goes through the bam read by read instead of looking for every position in the control file. should produce identical output to bam_basic_stats_pysam.py
     	    time python scripts_for_SediQuest/bam_basic_stats_pysam2.py \
	    	 --control {output.tmp_bed} \
		    --bam {input.bam} \
		    --tags lib --tags-fill {wildcards.indexlibid} \
		    --control-header scripts_for_SediQuest/probes_CONTROL_HEADER.txt  \
		    | gzip -c > {output.summary}
            """


#create plots kraken by burden score
rule plot_kraken_by_burden:
    input:
        burden=get_bed,
        kraken_target="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/kraken/{indexlibid}.byread",
        kraken_info_target="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.read_summary.txt.gz",
        kraken_deam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/kraken/{indexlibid}.byread",
        kraken_info_deam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.deam.read_summary.txt.gz",
    output:
        order="{project}_summary/{indexlibid}/{probeset}/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_kraken_order_byburden.pdf",
        fam_spe="{project}_summary/{indexlibid}/{probeset}/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_kraken_fam_spe_byburden.pdf"
    params:
        output_dir="{project}_summary/{indexlibid}/{probeset}/Mam_div_score_{score_b}/N_score_{score_n}/",
        n_score=score_n
    run:
        shell(f"""
			Rscript scripts_for_SediQuest/kraken_plot_by_burden.R {input.burden} {input.kraken_info_target} {input.kraken_target} {input.kraken_info_deam}  {input.kraken_deam} {params.output_dir}  {params.n_score}
		""")



#create plot combining cov and snps
rule plot_cov_snps:
    input:
        cov="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.uniq.L35MQ25_MD{score_b}_N{score_n}.cov",
        burden=get_bed,
    output:
        "{project}_summary/{indexlibid}/{probeset}/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_cov_MD.pdf"
    params:
        output_dir="{project}_summary/{indexlibid}/{probeset}/Mam_div_score_{score_b}/N_score_{score_n}/"
    run:
        shell(f"""
			Rscript scripts_for_SediQuest/combined_snps_cov.R {input.cov} {input.burden} {params.output_dir}  
		""")
