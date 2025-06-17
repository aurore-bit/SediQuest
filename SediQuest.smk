################################################################################
# 
#
# Aurore GALTIER, 25/11/2024
#
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

#define group to extract from kraken output
kraken_group = config["kraken_group"]

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

def get_control_info(wildcards):
    return probeset_df.loc[probeset_df['probeset'] == wildcards.probeset, "path_to_control_info"]



##############################################
#processing
##############################################

Final_bam = expand([f"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{{score_b}}/N_score_{score_n}/deam/{indexlibid}_uniqL35MQ25_MD{{score_b}}_N{score_n}_deam.bam"
         for indexlibid, probesets in INDEXLIBID.items()
         for probeset in probesets],
         score_b=score_b)

rule process_all:
    input:
        bam_files = Final_bam
    run:
        # Print the chosen indexlibid and probeset
        print('hey! Running target filtering, deamination, quality filtering, and duplicate removal :)')
        pass



rule uniq_q25_l35:
    input:
        bam="{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.bam"
    output:
        bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/{indexlibid}_uniqL35MQ25.bam",
        bai="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/{indexlibid}_uniqL35MQ25.bai"
    params:
        bam="../{indexlibid}.bam",
        unique_bam = "{indexlibid}.uniq.L35MQ25.bam",
        unique_bai= "{indexlibid}.uniq.L35MQ25.bam.bai",
        new_bam="{indexlibid}_uniqL35MQ25.bam",
        new_bai="{indexlibid}_uniqL35MQ25.bai",
        exhaust = "{indexlibid}.bam.exhaustion"
        #unique_sum="summary_stats_L35MQ25.txt"
    shell: """
    cd $(dirname {output.bam})
    /home/mmeyer/perlscripts/solexa/analysis/analyzeBAM.pl -nof -minlength 35 -qual 25 {params.bam} > {params.exhaust}
    mv {params.unique_bam} {params.new_bam}
    mv {params.unique_bai} {params.new_bai}
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
        Rscript scripts/filter_SNPs.R {params.b_score} {input.control_sites} {input.control} {output.sites_filtered} 
    """

rule filter_bam_by_control_sites:
    input:
      #  sites_filtered="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/sites_MD_{score_b}_N_{score_n}.filtered.txt",
        sites_filtered=get_control,
        bam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/{indexlibid}_uniqL35MQ25.bam",
    output:
        sites_bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{b_score}_N{score_n}.bam",
        bai="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{b_score}_N{score_n}.bai",
    threads: 1
    #conda: "envs/processing.yaml"
    shell: """
        bedtools intersect -a {input.bam} -b {input.sites_filtered} > {output.sites_bam} 
        samtools index {output.sites_bam} {output.bai}
    """


rule deam_filter:
    input: bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{b_score}_N{score_n}.bam"
    output:
            bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}_uniqL35MQ25_MD{b_score}_N{score_n}_deam.bam",
           # bai="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}_uniqL35MQ25_MD{b_score}_N{score_n}_deam.bai",
            stats="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}_uniqL35MQ25_MD{b_score}_N{score_n}_deam.deam53x3.stats"
    #conda: "envs/processing.yaml"
    params:
        unique_bam="{indexlibid}_uniqL35MQ25_MD{b_score}_N{score_n}.deam53x3.bam",
        new_bam="{indexlibid}_uniqL35MQ25_MD{b_score}_N{score_n}_deam.bam"
    shell: """
    /home/mmeyer/perlscripts/solexa/analysis/filterBAM.pl  -p5 0,1,2 -p3 0,-1,-2 -suffix deam53x3 {input.bam} &> {output.stats}
    mv {params.unique_bam}  {params.new_bam}
    mv {params.new_bam} $(dirname {output.bam})
    """



##############################################
#kraken step 1 create a fasta file
##############################################

fa_bam = [f"{project}/mappedbams/{indexlibid}/{probeset}/kraken/{indexlibid}.fa.gz"
       for indexlibid, probesets in INDEXLIBID.items()
       for probeset in probesets]

fa_split = [f"{project}/split/{indexlibid}/{probeset}/kraken/{indexlibid}.fa.gz"
         for indexlibid, probesets in INDEXLIBID.items()
         for probeset in probesets]

fa_bam_mapped_target = expand([f"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{{score_b}}/N_score_{score_n}/kraken/{indexlibid}.fa.gz"
                     for indexlibid, probesets in INDEXLIBID.items()
                     for probeset in probesets],
                     score_b=score_b)


fa_bam_mapped_target_deam = expand([f"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{{score_b}}/N_score_{score_n}/deam/kraken/{indexlibid}.fa.gz"
                         for indexlibid, probesets in INDEXLIBID.items()
                         for probeset in probesets],
                         score_b=score_b)

fa_bam_mapped = [f"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/kraken/{indexlibid}.fa.gz"
              for indexlibid, probesets in INDEXLIBID.items()
              for probeset in probesets]

wildcard_constraints:
    score_n="[^/]+"

rule kraken_step_1:
    input: fa_bam, fa_split, fa_bam_mapped_target, fa_bam_mapped_target_deam, fa_bam_mapped

rule bam_to_fasta_1:
    input:
        bam="{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.bam",
    output:
        fa_bam="{project}/mappedbams/{indexlibid}/{probeset}/kraken/{indexlibid}.fa.gz",
    threads: 1
    shell: """
        bam2fastx -a -A -Q {input.bam} | /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/lenfilter.pl 35 | gzip -c > {output.fa_bam}
    """


rule bam_to_fasta_2:
    input:
        split="{project}/split/{indexlibid}/{probeset}/{indexlibid}.bam",
    output:
        fa_split="{project}/split/{indexlibid}/{probeset}/kraken/{indexlibid}.fa.gz",
    threads: 1
    shell: """
         bam2fastx -a -A -Q {input.split} | /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/lenfilter.pl 35 | gzip -c > {output.fa_split}
       """



rule bam_to_fasta_3:
    input:
        bam_mapped_target="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.bam",
    output:
        fa_bam_mapped_target="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/kraken/{indexlibid}.fa.gz",
    threads: 1
    shell: """
         bam2fastx -a -A -Q {input.bam_mapped_target} | /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/lenfilter.pl 35 | gzip -c > {output.fa_bam_mapped_target}
       """

rule bam_to_fasta_4:
    input:
        bam_mapped="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/{indexlibid}_uniqL35MQ25.bam",
    output:
        fa_bam_mapped="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/kraken/{indexlibid}.fa.gz",
    threads: 1
    shell: """
         bam2fastx -a -A -Q {input.bam_mapped} | /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/lenfilter.pl 35 | gzip -c > {output.fa_bam_mapped}
       """


rule bam_to_fasta_5:
    input:
        bam_mapped_target_deam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam.bam",
    output:
        fa_bam_mapped_target_deam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/kraken/{indexlibid}.fa.gz",
    threads: 1
    shell: """
         bam2fastx -a -A -Q {input.bam_mapped_target_deam} | /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/lenfilter.pl 35 | gzip -c > {output.fa_bam_mapped_target_deam}
       """




##############################################
#kraken step 2 run kraken
##############################################
Kraken_split = [f"{project}/split/{indexlibid}/{probeset}/kraken/{indexlibid}.kraken_spc"
         for indexlibid, probesets in INDEXLIBID.items()
         for probeset in probesets]
         

Kraken_byread_split= [f"{project}/split/{indexlibid}/{probeset}/kraken/{indexlibid}.byread"
         for indexlibid, probesets in INDEXLIBID.items()
         for probeset in probesets]

Kraken_bam = [f"{project}/mappedbams/{indexlibid}/{probeset}/kraken/{indexlibid}.kraken_spc"
         for indexlibid, probesets in INDEXLIBID.items()
         for probeset in probesets]

Kraken_byread_bam = [f"{project}/mappedbams/{indexlibid}/{probeset}/kraken/{indexlibid}.byread"
         for indexlibid, probesets in INDEXLIBID.items()
         for probeset in probesets]

Kraken_bam_mapped= [f"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/kraken/{indexlibid}.kraken_spc"
         for indexlibid, probesets in INDEXLIBID.items()
         for probeset in probesets]

Kraken_byread_bam_mapped= [f"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/kraken/{indexlibid}.byread"
         for indexlibid, probesets in INDEXLIBID.items()
         for probeset in probesets]


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

Kraken_split_bam = expand([f"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{{score_b}}/N_score_{score_n}/split_kraken/{indexlibid}_uniqL35MQ25_MD{{score_b}}_N{score_n}_K{kraken_group}.bam"
         for indexlibid, probesets in INDEXLIBID.items()
         for probeset in probesets],
         score_b=score_b)

Kraken_split_deam_bam = expand([f"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{{score_b}}/N_score_{score_n}/deam/split_kraken/{indexlibid}_uniqL35MQ25_MD{{score_b}}_N{score_n}_deam_K{kraken_group}.bam"
         for indexlibid, probesets in INDEXLIBID.items()
         for probeset in probesets],
         score_b=score_b)
 
Kraken_split_capture_eff = [f"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/split_kraken/{indexlibid}_uniqL35MQ25_K{kraken_group}.bam"
         for indexlibid, probesets in INDEXLIBID.items()
         for probeset in probesets]

rule kraken_step_2:
    input: Kraken_split, Kraken_byread_split, Kraken_bam,Kraken_byread_bam, Kraken_bam_mapped,Kraken_byread_bam_mapped, Kraken_bam_mapped_target, Kraken_bam_mapped_target_deam, Kraken_byread_bam_mapped_target_deam, Kraken_byread_bam_mapped_target, Kraken_split, Kraken_split_deam_bam,Kraken_split_bam


rule generic_kraken:
    input: fa="{indexlibid}.fa.gz"
    output: kraken="{indexlibid}.kraken"
    threads: 1
    conda: "envs/kraken.yaml"
    shell: """

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


    echo 'making report..'
    ## translate file has info on each read (but not in an easy-to-digest way)
    ## phylo file has a summary for each level in the taxonomy
    python3 /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/kraken_report.py --db /mnt/ramdisk/refseqReleaseKraken {input.kraken} > {output.phylo}
    """

rule generic_kraken_spc_summary:
    input: spc_krak="{indexlibid}.kraken_phylo",
        #spc_groups="/mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/data/major_species_groups.txt"
    output: spc_summary="{indexlibid}.kraken_spc"
  #  conda: "envs/kraken.yaml"
    shell: """
    
    v=$((grep '  Afrotheria$' {input.spc_krak} || echo 0 0) | awk '{{print $2}}')
    echo "$v" {input.spc_krak}
    
    w=$((grep '  Primates$' {input.spc_krak} || echo 0 0) | awk '{{print $2}}')
    echo "$v: $w" {input.spc_krak}
    
    x=$((grep '  Laurasiatheria$' {input.spc_krak} || echo 0 0) | awk '{{print $2}}')
    echo "$v: $w: $x" {input.spc_krak}
    
    y=$((grep '  Glires$' {input.spc_krak} || echo 0 0) | awk '{{print $2}}')
    echo "$v: $w: $x: $y" {input.spc_krak}

    m=$((grep '  Mammalia$' {input.spc_krak} ||echo 0 0) | awk '{{print $2-'$v'-'$w'-'$x'-'$y'}}')
    echo " $v: $w: $x: $y: $m : " {input.spc_krak}

    b=$((grep '  Bacteria$' {input.spc_krak} || echo 0 0) | awk '{{print $2}}')
    echo "$v: $w: $x: $y : $m : $b" {input.spc_krak}

    s=$((grep '  Sauropsida$' {input.spc_krak} || echo 0 0) | awk '{{print $2}}')
    echo "$v: $w: $x: $y : $m : $b : $s :" {input.spc_krak}

    r=$(grep '	root$' {input.spc_krak} | awk '{{print $2-'$m'-'$b'-'$s'-'$v'-'$w'-'$x'-'$y'}}')
    echo "$v: $w: $x: $y : $m : $b : $s : $r :" {input.spc_krak}

    u=$(grep '	unclassified$' {input.spc_krak} | awk '{{print $2}}')
    echo "$v: $w: $x: $y : $m : $b : $s : $r : $u : " {input.spc_krak}

    echo "Mammalia $m" >> {output.spc_summary}
    echo "Bacteria $b" >> {output.spc_summary}
    echo "Sauropsida $s" >> {output.spc_summary}
    echo "root $r" >> {output.spc_summary}
    echo "unclassified $u" >> {output.spc_summary}
    echo "Afrotheria $v" >> {output.spc_summary}
    echo "Primates $w" >> {output.spc_summary}
    echo "Laurasiatheria $x" >> {output.spc_summary}
    echo "Glires $y" >> {output.spc_summary}

"""


rule generic_kraken_translate:
    input: kraken="{indexlibid}.kraken"
    output: translate="{indexlibid}.translate"
    threads: 1
    shell: """

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
#split kraken after deam
##############################################

rule generic_kraken_extract_after_deam:
    input: bam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam.bam",
           kraken = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/kraken/{indexlibid}.kraken",
           translate = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/kraken/{indexlibid}.translate"
    output: bam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/split_kraken/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam_K{kraken_group}.bam"
    wildcard_constraints:
       # kcat = "|\.swp_third|\.swp_original|\.swp_masked",
        kraken_group = "Primates"
    threads: 1
    shell: """
   

    echo 'finding clade'
    clade=$(awk '$2 == "'{wildcards.kraken_group}'" {{print $1}}' /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/data/clade_taxa_map_alphanum.txt)
    echo clade $clade
    ofile=$(dirname {output.bam})/$(basename {output.bam} .bam)
    echo ofile $ofile
    
    ## save the header first (doesn't kraken_report just overwrite this?)
    # samtools view -H -b {input.bam} > {output.bam}

    echo 'splitting bam..'
    ## translate file has info on each read (but not in an easy-to-digest way)
    ## have to write to the correct file - this doesn't take an ofile, which is maddening...
    python3 /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/kraken_report.py --db /mnt/ramdisk/refseqReleaseKraken {input.kraken} --extractFile {input.bam} \
       --clades $clade --extract-out-base $ofile  > /dev/null

    """







##############################################
#summary
##############################################

summary_damage= expand([f"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{{score_b}}/N_score_{score_n}/deam_stats/{indexlibid}_uniqL35MQ25_MD{{score_b}}_N{score_n}.summary_damage.txt"
        for indexlibid, probesets in INDEXLIBID.items()
        for probeset in probesets],
        score_b=score_b)

summary_table= expand([f"{project}_summary/{indexlibid}/{probeset}/Mam_div_score_{{score_b}}/N_score_{score_n}/{indexlibid}.pipeline_summary.txt"
        for indexlibid, probesets in INDEXLIBID.items()
        for probeset in probesets],
        score_b=score_b)

faunal_contam = expand([f"{project}_summary/{indexlibid}/{probeset}/Mam_div_score_{{score_b}}/N_score_{score_n}/{indexlibid}_SNPs_number_cov_combined.pdf"
        for indexlibid, probesets in INDEXLIBID.items()
        for probeset in probesets],
        score_b=score_b)


rule summaries:
    input:
        summary_damage=summary_damage,
        summary_table=summary_table,
        cov_plot=faunal_contam
    run:
        print('hey! Creating summaries')
        pass

#create coverage file for on target reads
rule cov_bed_files:
    input: "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.bam"
    output:
            cov="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.cov"
    shell:"""
            samtools mpileup -B {input} > {output.cov}
            """
""


#give average coverage deaminated not primates filtered 
rule cov_average_deam_all:
    input: 
            bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam.bam",
           # bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam_clean.bam",
            bed = "/mnt/archgen/Reference_Genomes/Human/hs37d5/SNPCapBEDs/1240K.pos.list_hs37d5.0based.bed",
            ref=lambda wildcards: get_ref(wildcards).replace("bwa-0.4.9", "whole_genome.fa")
    output:
            depth="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam.depth",
            count="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam.count"
    shell:"""
            /mnt/expressions/yaniv/Software/samtools-1.17/samtools mpileup -R -B -q25 -Q30 -l {input.bed} -f {input.ref} {input.bam} > {output.depth}
            awk '$4 ==1 {{count++}} END {{print count}}' {output.depth} > {output.count}
            """

#give average coverage deaminated primates filtered 
rule cov_average_deam_primates:
    input: 
            bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/split_kraken/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam_KPrimates.bam",
            bed = "/mnt/archgen/Reference_Genomes/Human/hs37d5/SNPCapBEDs/1240K.pos.list_hs37d5.0based.bed",
            ref=lambda wildcards: get_ref(wildcards).replace("bwa-0.4.9", "whole_genome.fa")

    output:
            depth="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/split_kraken/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam_KPrimates.depth",
            count="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/split_kraken/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam_KPrimates.count"
    shell:"""
            /mnt/expressions/yaniv/Software/samtools-1.17/samtools mpileup -R -B -q25 -Q30 -l {input.bed} -f {input.ref} {input.bam} > {output.depth}
            awk '$4 ==1 {{count++}} END {{print count}}' {output.depth} > {output.count}
            """


#give average coverage all primates filtered 
rule cov_average_primates:
    input: 
            bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/split_kraken/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_KPrimates.bam",
            bed = "/mnt/archgen/Reference_Genomes/Human/hs37d5/SNPCapBEDs/1240K.pos.list_hs37d5.0based.bed",
            ref=lambda wildcards: get_ref(wildcards).replace("bwa-0.4.9", "whole_genome.fa")

    output:
            depth="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/split_kraken/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_KPrimates.depth",
            count="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/split_kraken/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_KPrimates.count"
    shell:"""
            /mnt/expressions/yaniv/Software/samtools-1.17/samtools mpileup -R -B -q25 -Q30 -l {input.bed} -f {input.ref} {input.bam} > {output.depth}
            awk '$4 ==1 {{count++}} END {{print count}}' {output.depth} > {output.count}
            """



#give count coverage all not primates filtered 
rule cov_count_all:
    input: 
            bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.bam",
            bed="/mnt/archgen/Reference_Genomes/Human/hs37d5/SNPCapBEDs/1240K.pos.list_hs37d5.0based.bed",
            ref=lambda wildcards: get_ref(wildcards).replace("bwa-0.4.9", "whole_genome.fa")

    output:
            depth="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.depth",
            count="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.count"
    shell:"""
            /mnt/expressions/yaniv/Software/samtools-1.17/samtools mpileup -R -B -q25 -Q30 -l {input.bed} -f {input.ref} {input.bam} > {output.depth}
            awk '$4 ==1 {{count++}} END {{print count}}' {output.depth} > {output.count}
            """


#give count coverage for HO panel primates deaminated
rule cov_count_HO_primates_deam:
    input: 
            bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/split_kraken/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam_KPrimates.bam",
            bed="/mnt/expressions/Aurore/merlin_malta_project/v62.0_HO_public.snp.txt.bed",
            ref=lambda wildcards: get_ref(wildcards).replace("bwa-0.4.9", "whole_genome.fa")

    output:
            depth="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/split_kraken/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam_KPrimates.depth_ho",
            count="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/split_kraken/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam_KPrimates.count_ho"
    shell:"""
            /mnt/expressions/yaniv/Software/samtools-1.17/samtools mpileup -R -B -q25 -Q30 -l {input.bed} -f {input.ref} {input.bam} > {output.depth}
            awk '$4 ==1 {{count++}} END {{print count}}' {output.depth} > {output.count}
            """

#give count coverage for HO panel deaminated
rule cov_count_HO_deam:
    input: 
            bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam.bam",
            bed="/mnt/expressions/Aurore/merlin_malta_project/v62.0_HO_public.snp.txt.bed",
            ref=lambda wildcards: get_ref(wildcards).replace("bwa-0.4.9", "whole_genome.fa")

    output:
            depth="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam.depth_ho",
            count="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam.count_ho"
    shell:"""
            /mnt/expressions/yaniv/Software/samtools-1.17/samtools mpileup -R -B -q25 -Q30 -l {input.bed} -f {input.ref} {input.bam} > {output.depth}
            awk '$4 ==1 {{count++}} END {{print count}}' {output.depth} > {output.count}
            """


#give count coverage for HO panel primates
rule cov_count_HO_primates:
    input: 
            bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/split_kraken/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_KPrimates.bam",
            bed="/mnt/expressions/Aurore/merlin_malta_project/v62.0_HO_public.snp.txt.bed",
            ref=lambda wildcards: get_ref(wildcards).replace("bwa-0.4.9", "whole_genome.fa")
    output:
            depth="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/split_kraken/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_KPrimates.depth_ho",
            count="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/split_kraken/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_KPrimates.count_ho"
    shell:"""
            /mnt/expressions/yaniv/Software/samtools-1.17/samtools mpileup -R -B -q25 -Q30 -l {input.bed} -f {input.ref} {input.bam} > {output.depth}
            awk '$4 ==1 {{count++}} END {{print count}}' {output.depth} > {output.count}
            """

#give count coverage for HO panel 
rule cov_count_HO:
    input: 
            bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.bam",
            bed="/mnt/expressions/Aurore/merlin_malta_project/v62.0_HO_public.snp.txt.bed",
            ref=lambda wildcards: get_ref(wildcards).replace("bwa-0.4.9", "whole_genome.fa")

    output:
            depth="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.depth_ho",
            count="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.count_ho"
    shell:"""
            /mnt/expressions/yaniv/Software/samtools-1.17/samtools mpileup -R -B -q25 -Q30 -l {input.bed} -f {input.ref} {input.bam} > {output.depth}
            awk '$4 ==1 {{count++}} END {{print count}}' {output.depth} > {output.count}
            """



#give an estimate of the contamination
rule contamination_estimaTe:
    input: 
        bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.bam"
    output:
        bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.authentict"
    shell: """
            samtools view {input.bam} | AuthentiCT deam2cont -o {output.bam} -s 10000 -
            """

#just count number of reads
rule pipeline_summary:
    input:
        map_bam="{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.bam",
        target_bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.bam",
        rmdup_bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/{indexlibid}_uniqL35MQ25.bam",
        deam_bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam.bam",
        split_kraken="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/split_kraken/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam_KPrimates.bam",
        exhaustion="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/{indexlibid}.bam.exhaustion",
        deam_stats="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam_stats/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.summary_damage.txt",
        average_all_deam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam.count",
        average_primates_deam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/split_kraken/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam_KPrimates.count",
        average_all="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.count",
        average_primates="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/split_kraken/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_KPrimates.count",
        average_all_deam_ho="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam.count_ho",
        average_primates_deam_ho="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/split_kraken/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam_KPrimates.count_ho",
        average_primates_ho="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/split_kraken/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_KPrimates.count_ho",
        average_all_ho="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.count_ho",
        contam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.authentict"
    output:
        summary_annotated="{project}_summary/{indexlibid}/{probeset}/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}.pipeline_summary.txt",
    shell:
        """
        splitbam="{wildcards.project}/split/{wildcards.indexlibid}/{wildcards.probeset}/{wildcards.indexlibid}.bam"
        
        if [ -e $splitbam ] ; then
          sc=$(samtools view -c $splitbam)
        else
          sc="NA"
        fi

        # Get read counts for each input BAM file
        map_bam_reads=$(/home/bioinf/usr/bin/samtools view -c {input.map_bam})
        rmdup_bam_reads=$(/home/bioinf/usr/bin/samtools view -c {input.rmdup_bam})
        target_bam_reads=$(/home/bioinf/usr/bin/samtools view -c {input.target_bam})
        deam_bam_reads=$(/home/bioinf/usr/bin/samtools view -c {input.deam_bam})
        primates=$(/home/bioinf/usr/bin/samtools view -c {input.split_kraken})

        # Get deam stats (5'CT_95CI, 3'CT_95CI, cond5'CT_95CI, cond3'CT_95CI) from deam_stats file
        deam_stats_file="{input.deam_stats}"
        deam_5CT_95CI=$(awk 'NR==2 {{print $4}}' $deam_stats_file)
        deam_3CT_95CI=$(awk 'NR==2 {{print $5}}' $deam_stats_file)
        deam_cond5CT_95CI=$(awk 'NR==2 {{print $10}}' $deam_stats_file)
        deam_cond3CT_95CI=$(awk 'NR==2 {{print $11}}' $deam_stats_file)

        # Get %unique and %exhausted from exhaustion file
        exhaustion_file="{input.exhaustion}"
        unique_percentage=$(awk 'NR==2 {{print $8}}' $exhaustion_file)
        exhausted_percentage=$(awk 'NR==2 {{print $9}}' $exhaustion_file)

        #get contamination estimate
       contamination_file="{input.contam}"
       contamination=$(awk 'NR==8 {{print $2}}' $contamination_file)
        err_estimate=$(awk 'NR==8 {{print $3}}' $contamination_file)

        tmp={output.summary_annotated}

        # Write header to the summary file
        echo IndexLibID N_score MD_score probeset split mapped rmdup target deam primates 5'CT_95CI 3'CT_95CI cond5'CT_95CI cond3'CT_95CI %unique %exhausted Count_Krall Count_KrPrimates Count_KrPrimates_deam Count_KrALL_deam Count_Krall_ho Count_KrPrimates_ho Count_KrPrimates_deam_ho Count_KrALL_deam_ho Contamination Contamination_err_estimate | tr ' ' '\t' > $tmp

        # Append the values for each field to the summary file
        echo {wildcards.indexlibid} {wildcards.score_n} {wildcards.score_b} {wildcards.probeset} \
        $sc \
        $map_bam_reads \
        $rmdup_bam_reads \
        $target_bam_reads \
        $deam_bam_reads \
        $primates \
        "$deam_5CT_95CI" \
        "$deam_3CT_95CI" \
        "$deam_cond5CT_95CI" \
        "$deam_cond3CT_95CI" \
        "$unique_percentage" \
        "$exhausted_percentage" \
        $(cat {input.average_all}) \
        $(cat {input.average_primates}) \
        $(cat {input.average_primates_deam}) \
        $(cat {input.average_all_deam}) \
        $(cat {input.average_all_ho}) \
        $(cat {input.average_primates_ho}) \
        $(cat {input.average_primates_deam_ho}) \
        $(cat {input.average_all_deam_ho}) \
        "$contamination" \
        "$err_estimate" \
        | tr ' ' '\t' >> $tmp
        """


#create a file indicating every info for each read
rule bam_read_summary_target:
     input: bam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.bam",
            control=get_control_info
     output: summary ="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.read_summary.txt.gz"
     shell: """
            ## bam_basic_stats_pysam2.py is much faster - it goes through the bam read by read instead of looking for every position in the control file. should produce identical output to bam_basic_stats_pysam.py
     	    time /home/benjamin_vernot/miniconda3/bin/python /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/bam_basic_stats_pysam2.py \
	    	 --control {input.control} \
		    --bam {input.bam} \
		    --tags lib --tags-fill {wildcards.indexlibid} \
		    --control-header /mnt/expressions/benjamin_vernot/soil_capture_2017/site_categories_for_capture/soil_probe_designs/probes_CONTROL_HEADER.txt \
		    | gzip -c > {output.summary}
            """

rule bam_read_summary_deam:
     input: bam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam.bam",
            control=get_control_info
     output: summary ="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam.read_summary.txt.gz"
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
rule deam_stats:
    input: bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}_deam.bam",
    output: stats="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam_stats/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.L30-1000_MQ0.3p_dinucleotide_refbase_composition.txt ",
    params: bam="../{indexlibid}_uniqL35MQ25_MD{score_b}_N{n_score}.bam"
    shell: """
            cd $(dirname {output.plots})
            /home/mmeyer/perlscripts/solexa/analysis/substitution_patterns.pl {params.bam} 
            """

#new deam summaries files
rule deam_stats_table:
    input: bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.bam",
    output: summary_new="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/deam_stats/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.summary_damage.txt"
   # params:  summary_new="{indexlibid}_uniqL35MQ25_MD{score_b}_N{n_score}.summary_damage.txt"
    shell: """
            /home/mmeyer/perlscripts/solexa/analysis/quick_substitutions.pl {input.bam} > {output.summary_new}
            """


#create plots kraken by burden score
rule plot_kraken_by_burden:
    input:
        burden=get_bed,
        kraken_target="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/{indexlibid}_n_{n_score}.byread",
        kraken_info_target="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/{indexlibid}_n_{n_score}.read_summary.txt.gz",
        kraken_deam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/{indexlibid}_n_{n_score}.byread",
        kraken_info_deam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/{indexlibid}_n_{n_score}.read_summary.txt.gz",
    output:
        order="{project}_summary/{indexlibid}/{probeset}/{score_b}_filter/{indexlibid}_n_{n_score}_kraken_order_byburden.pdf",
        fam_spe="{project}_summary/{indexlibid}/{probeset}/{score_b}_filter/{indexlibid}_n_{n_score}_kraken_fam_spe_byburden.pdf"
    params:
        output_dir="{project}_summary/{indexlibid}/{probeset}/{score_b}_filter/",
        n_score=score_n
    run:
        shell(f"""
			Rscript scripts/kraken_plot_by_burden.R {input.burden} {input.kraken_info_target} {input.kraken_target} {input.kraken_info_deam}  {input.kraken_deam} {params.output_dir}  {params.n_score}
		""")



#create plot combining cov and snps
rule plot_cov_snps:
    input:
        cov="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_uniqL35MQ25_MD{score_b}_N{score_n}.bam.cov",
        burden=get_bed,
    output:
        "{project}_summary/{indexlibid}/{probeset}/Mam_div_score_{score_b}/N_score_{score_n}/{indexlibid}_SNPs_number_cov_combined.pdf"
    params:
        output_dir="{project}_summary/{indexlibid}/{probeset}/Mam_div_score_{score_b}/N_score_{score_n}/"
    run:
        shell(f"""
			Rscript scripts/combined_snps_cov.R {input.cov} {input.burden} {params.output_dir}  
		""")