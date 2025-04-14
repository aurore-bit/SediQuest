################################################################################
# sedquest: nuclear capture sediment pipeline
#
# Aurore GALTIER, 14/04/25
#
#FOR NOW YOU NEED TO:
#-have a config folder with config file, probeset information and a sample table
#-copy split bam (raw file) in /output_v0/split/{indexlibid}/{probeset}/
#-copy mapped bam in /output_v0/mappedbams/{indexlibid}/{probeset}/

#WHAT IT DOES
#-check_mapped rule to check if the reference genome was indeed mapped to the good reference (third allele reference genome)
#-process_all rule to run unique_qual_length target and deam filtering (in this order)
#-kraken rule to run all the kraken steps
#-summaries to create summaries plot and table

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

#define n_score
score_n = config["score_n"]

#define group to extract from kraken output
kraken_group = config["kraken_group"]

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
# This rule is to sort your bam and check 
#that you mapped it to the right reference
##############################################
sorted_bams = expand('{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.sorted.bam',
    indexlibid=INDEXLIBID,
    probeset=indexlibid_df.loc[INDEXLIBID, "probeset"].unique(),
    project=project)

rule check_mapped:
    input:
        bam_files = sorted_bams,
    run:
        print("Hey, sorting and checking bam file")
        pass


rule samtools:
    input:
        mapped = "{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.bam"
    output:
        sorted = "{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.sorted.bam"
    params:
        ref = get_ref
    threads: 8
    shell:
        """
        python3.7 scripts/check_mapped_ref.py {input.mapped} {params.ref}
        /home/visagie/.local/bin/samtools sort --threads {threads} -@30 -o {output.sorted} {input.mapped}
        """


##############################################
#processing rule to extract reads with quality, 
#length filtering, duplicate removal
# on site target filtering
#and C->T substitution filtering 
##############################################

rule process_all:
    input: 
        expand([ "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/{indexlibid}_n_{n_score}.bam",
                "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/{indexlibid}_n_{n_score}.bai"
                ],
            indexlibid=INDEXLIBID,
            probeset=indexlibid_df.loc[INDEXLIBID, "probeset"],
            project=project,
            score_b=score_b,
            n_score=score_n)
    run:
        print('hey! Running target filtering, deamination, quality filtering and duplicate removal :)')
        pass


rule uniq_q25_l35:
    input:
        bam="{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.bam"
    output:
        bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/{indexlibid}.bam",
        bai="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/{indexlibid}.bai",
        exhaust = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/{indexlibid}.bam.exhaustion"
    params:
        unique_bam = "{indexlibid}.uniq.L35MQ25.bam",
        unique_bai= "{indexlibid}.uniq.L35MQ25.bam.bai",
        new_bam="{indexlibid}.bam",
        new_bai="{indexlibid}.bai",
        #unique_sum="summary_stats_L35MQ25.txt"
    shell: """
    /home/mmeyer/perlscripts/solexa/analysis/analyzeBAM.pl -nof -minlength 35 -qual 25 {input.bam} > {output.exhaust}
    mv {params.unique_bam} {params.new_bam}
    mv {params.unique_bai} {params.new_bai}
    mv {params.new_bam}  $(dirname {output.bam})
    mv {params.new_bai}  $(dirname {output.bam})
           """

#this is needed if you want to filter for burden score!!!!! 
rule filter_control_sites_b_score:
    input:
        control_sites = get_bed,
        control= get_control,
    output:
         sites_filtered = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/sites_b_{score_b}_n_{n_score}.filtered.txt",
    params:
       b_score = "{score_b}",
       n_score = "{n_score}",
    threads: 1
    conda: "envs/processing.yaml"
    shell: """
        Rscript scripts/filter_SNPs.R {params.b_score} {input.control_sites} {input.control} {output.sites_filtered} {params.n_score}
    """

rule filter_bam_by_control_sites:
    input:
        sites_filtered="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/sites_b_{score_b}_n_{n_score}.filtered.txt",
        #this is if you don't want to filter for burden score
	#sites_filtered=get_control,
        bam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/{indexlibid}.bam",
    output:
        sites_bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/{indexlibid}_n_{n_score}.bam",
        bai="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/{indexlibid}_n_{n_score}.bai",
    threads: 1
    #conda: "envs/processing.yaml"
    shell: """
        bedtools intersect -a {input.bam} -b {input.sites_filtered} > {output.sites_bam} 
        samtools index {output.sites_bam} {output.bai}
    """


rule deam_filter:
    input: bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/{indexlibid}_n_{n_score}.bam"
    output:
            bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/{indexlibid}_n_{n_score}.bam",
            bai="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/{indexlibid}_n_{n_score}.bai",
            stats="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/{indexlibid}_n_{n_score}.deam53x3.stats"
    #conda: "envs/processing.yaml"
    params:
        unique_bam="{indexlibid}_n_{n_score}.deam53x3.bam",
        new_bam="{indexlibid}_n_{n_score}.bam"
    shell: """
    /home/mmeyer/perlscripts/solexa/analysis/filterBAM.pl  -p5 0,1,2 -p3 0,-1,-2 -suffix deam53x3 {input.bam} &> {output.stats}
    mv {params.unique_bam}  {params.new_bam}
    mv {params.new_bam} $(dirname {output.bam})
    samtools index {output.bam} {output.bai}
    """



##############################################
#kraken rule to produce all kraken outputs
##############################################

rule kraken:
    input:
        expand(["{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/{indexlibid}_n_{n_score}.kraken_spc",
            "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/{indexlibid}_n_{n_score}.byread",
            "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/{indexlibid}.kraken_spc",
            "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/{indexlibid}.byread",
            "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/{indexlibid}_n_{n_score}.kraken_spc",
            "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/{indexlibid}_n_{n_score}.byread",
            "{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.kraken_spc",
            "{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.byread",
            "{project}/split/{indexlibid}/{probeset}/{indexlibid}.kraken_spc",
            "{project}/split/{indexlibid}/{probeset}/{indexlibid}.byread",
            "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/split_kraken/{indexlibid}_n_{n_score}.k_Primates.bam",
            "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/split_kraken/deam/{indexlibid}_n_{n_score}.k_{kraken_group}.bam"],
            indexlibid=INDEXLIBID,
            probeset=indexlibid_df.loc[INDEXLIBID, "probeset"],
            project=project,
            score_b=score_b,
            n_score=score_n,
            kraken_group = kraken_group)
    run:
        print('hey! Running Kraken')
        pass


rule bam_to_fasta:
    input:
        bam="{indexlibid}_n_{n_score}.bam"
    output:
        fa="{indexlibid}_n_{n_score}.fa.gz"
    threads: 1
    shell: """
     bam2fastx -a -A -Q {input.bam} | /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/lenfilter.pl 35 | gzip -c > {output.fa}
     """

rule generic_kraken:
    input: fa="{indexlibid}_n_{n_score}.fa.gz"
    output: kraken="{indexlibid}_n_{n_score}.kraken"
    threads: 1
    shell: """

    echo THIS FAILS IF THE KRAKEN DB IS NOT COPIED - fix it

     nproc=$(time ~frederic_romagne/kraken/install/kraken --threads {threads} --db /mnt/ramdisk/refseqReleaseKraken \
     --output {output.kraken} {input.fa} 2>&1 | grep 'processed' | cut -f1 -d' ')
     echo "count seqs {input.fa}"
     nseq=$(gunzip -c {input.fa} | grep -c -e '>' -e '@') || echo "no seqs found $nseq"

     if [ ! $nseq -eq $nproc ] ; then echo KRAKEN RUN FAILED; stop; fi
    """

rule generic_kraken_summary:
    input: kraken="{indexlibid}_n_{n_score}.kraken"
    output: phylo="{indexlibid}_n_{n_score}.kraken_phylo"
    threads: 1
    shell: """

    echo 'making report..'
    python3 /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/kraken_report.py --db /mnt/ramdisk/refseqReleaseKraken {input.kraken} > {output.phylo}
    """

rule generic_kraken_spc_summary:
    input: spc_krak="{indexlibid}_n_{n_score}.kraken_phylo",
    output: spc_summary="{indexlibid}_n_{n_score}.kraken_spc"
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
    input: kraken="{indexlibid}_n_{n_score}.kraken"
    output: translate="{indexlibid}_n_{n_score}.translate"
    threads: 1
    shell: """

    echo 'making report..'
    python3 /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/kraken_report.py --db /mnt/ramdisk/refseqReleaseKraken {input.kraken} --translate {output.translate} > /dev/null

    """

rule generic_kraken_byread:
    input: translate="{indexlibid}_n_{n_score}.translate"
    output: byread="{indexlibid}_n_{n_score}.byread"
    threads: 1
    shell: """

     echo 'parsing translate file..'
     python3 /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/parse_kraken_translate_file.py --translate-file {input.translate} > {output.byread}

    """


##############################################
#extract kraken primate reads after deam
#filtering
##############################################

rule generic_kraken_extract_after_deam:
    input: bam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/{indexlibid}_n_{n_score}.bam",
           kraken = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/{indexlibid}_n_{n_score}.kraken",
           translate = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/{indexlibid}_n_{n_score}.translate"
    output: bam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/split_kraken/{indexlibid}_n_{n_score}.k_{kraken_group}.bam"
    wildcard_constraints:
    threads: 1
    shell: """
    echo 'finding clade'
    clade=$(awk '$2 == "'{wildcards.kraken_group}'" {{print $1}}' /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/data/clade_taxa_map_alphanum.txt)
    echo clade $clade
    ofile=$(dirname {output.bam})/$(basename {output.bam} .bam)
    echo ofile $ofile
    echo 'splitting bam..'
    python3 /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/kraken_report.py --db /mnt/ramdisk/refseqReleaseKraken {input.kraken} --extractFile {input.bam} \
       --clades $clade --extract-out-base $ofile  > /dev/null

    """

##############################################
#extract kraken primate reads before deam                                                                                            
#filtering
##############################################


rule generic_kraken_extract_before_deam:
    input: bam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/{indexlibid}_n_{n_score}.bam",
           kraken = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/{indexlibid}_n_{n_score}.kraken",
           translate = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/{indexlibid}_n_{n_score}.translate"
    output: bam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/split_kraken/{indexlibid}_n_{n_score}.k_{kraken_group}.bam"
    wildcard_constraints:
    threads: 1
    shell: """
    echo 'finding clade'
    clade=$(awk '$2 == "'{wildcards.kraken_group}'" {{print $1}}' /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/data/clade_taxa_map_alphanum.txt)
    echo clade $clade
    ofile=$(dirname {output.bam})/$(basename {output.bam} .bam)
    echo ofile $ofile
    echo 'splitting bam..'
    python3 /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/kraken_report.py --db /mnt/ramdisk/refseqReleaseKraken {input.kraken} --extractFile {input.bam} \
       --clades $clade --extract-out-base $ofile  > /dev/null

    """

rule deam_filter_after_kraken:
    input: bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/split_kraken/{indexlibid}_n_{n_score}.k_{kraken_group}.bam"
    output:
            bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/split_kraken/deam/{indexlibid}_n_{n_score}.k_{kraken_group}.bam",
            bai="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/split_kraken/deam/{indexlibid}_n_{n_score}.k_{kraken_group}.bai",
            stats="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/split_kraken/deam/{indexlibid}_n_{n_score}.k_{kraken_group}.deam53x3.stats"
    params:
        unique_bam="{indexlibid}_n_{n_score}.k_{kraken_group}.deam53x3.bam",
        new_bam="{indexlibid}_n_{n_score}.k_{kraken_group}.bam"
    shell: """
    /home/mmeyer/perlscripts/solexa/analysis/filterBAM.pl  -p5 0,1,2 -p3 0,-1,-2 -suffix deam53x3 {input.bam} &> {output.stats}
    mv {params.unique_bam}  {params.new_bam}
    mv {params.new_bam} $(dirname {output.bam})
    samtools index {output.bam} {output.bai}
    """


##############################################
#summary rules to produce tables and plots
##############################################
rule summaries:
    input:
        expand([#"{project}_summary/{indexlibid}/{probeset}/{score_b}_filter/{indexlibid}_n_{n_score}_coverage_plot.pdf",
        "{project}_summary/{indexlibid}/{probeset}/{score_b}_filter/{indexlibid}_n_{n_score}.k_{kraken_group}.new.pipeline_summary.txt",
       # "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam_stats/{indexlibid}_n_{n_score}.summary_damage.txt",
       # "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/{indexlibid}_n_{n_score}.read_summary.txt.gz",
        #"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/{indexlibid}_n_{n_score}.read_summary.txt.gz",
       # "{project}_summary/{indexlibid}/{probeset}/{score_b}_filter/{indexlibid}_n_{n_score}_kraken_order_byburden.pdf",
        #"{project}_summary/{indexlibid}/{probeset}/{score_b}_filter/{indexlibid}_n_{n_score}_kraken_fam_spe_byburden.pdf",
        #"{project}_summary/{indexlibid}/{probeset}/{score_b}_filter/{indexlibid}_n_{n_score}_SNPs_number_cov_combined.pdf",
       # "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/split_kraken/deam_stats/{indexlibid}_n_{n_score}.k_{kraken_group}.summary_damage.txt"
       "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam_stats/{indexlibid}_n_{n_score}_new.summary_damage.txt"],
        indexlibid=INDEXLIBID,
        probeset=indexlibid_df.loc[INDEXLIBID, "probeset"],
        project=project,
        score_b=score_b,
        n_score=score_n,
        kraken_group = kraken_group)
    run:
        print('hey! Creating summaries')
        pass

#create coverage file
rule cov_bed_files:
    input: "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/{indexlibid}_n_{n_score}.bam"
    output:
            cov="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/{indexlibid}_n_{n_score}.cov"
    shell:"""
            samtools mpileup -B {input} > {output.cov}
            """

#just count number of reads
rule pipeline_summary:
    input:
        map_bam="{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.bam",
        target_bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/{indexlibid}_n_{n_score}.bam",
        rmdup_bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/{indexlibid}.bam",
       # rmdup_stats="summary_stats_L35MQ25.txt",
        deam_bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/{indexlibid}_n_{n_score}.bam",
        split_kraken="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/split_kraken/{indexlibid}_n_{n_score}.k_{kraken_group}.bam",
        deam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam_stats/{indexlibid}_n_{n_score}_new.summary_damage.txt"
    output:
        summary_annotated="{project}_summary/{indexlibid}/{probeset}/{score_b}_filter/{indexlibid}_n_{n_score}.k_{kraken_group}.new.pipeline_summary.txt"
    shell:
        """
        splitbam="{wildcards.project}/split/{wildcards.indexlibid}/{wildcards.probeset}/{wildcards.indexlibid}.bam"
        #ls $splitbam

        if [ -e $splitbam ] ; then
          sc=$(samtools view -c $splitbam)
        else
          sc="NA"
        fi

        # Get read counts for each input BAM file
        map_bam_reads=$(samtools view -c {input.map_bam})
        rmdup_bam_reads=$(samtools view -c {input.rmdup_bam})
        target_bam_reads=$(samtools view -c {input.target_bam})
        deam_bam_reads=$(samtools view -c {input.deam_bam})
        primates=$(samtools view -c {input.split_kraken})

        # Get deam stats (5'CT, 3'CT, 5'CT_95CI, 3'CT_95CI, cond5'CT_95CI, cond3'CT_95CI) from deam_stats file
        deam_stats_file="{input.deam}"
        deam_5CT=$(awk 'NR==2 {{print $2}}' $deam_stats_file)
        deam_3CT=$(awk 'NR==2 {{print $3}}' $deam_stats_file)
        deam_5CT_95CI=$(awk 'NR==2 {{print $4}}' $deam_stats_file)
        deam_3CT_95CI=$(awk 'NR==2 {{print $5}}' $deam_stats_file)
        deam_cond5CT_95CI=$(awk 'NR==2 {{print $10}}' $deam_stats_file)
        deam_cond3CT_95CI=$(awk 'NR==2 {{print $11}}' $deam_stats_file)

        tmp={output.summary_annotated}

        echo IndexLibID probeset split mapped rmdup target deam primates  | tr ' ' '\t' > $tmp

        echo {wildcards.indexlibid} {wildcards.n_score} {wildcards.score_b} {wildcards.probeset} \
        $sc \
        $map_bam_reads \
        $rmdup_bam_reads \
        $target_bam_reads \
        $deam_bam_reads \
        $primates \
        "$deam_5CT" \
        "$deam_3CT" \
        "$deam_5CT_95CI" \
        "$deam_3CT_95CI" \
        "$deam_cond5CT_95CI" \
        "$deam_cond3CT_95CI" \
        | tr ' ' '\t' >> $tmp
        """

#create a file indicating every info for each read
rule bam_read_summary_target:
     input: bam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/{indexlibid}_n_{n_score}.bam",
            control=get_control_info
     output: summary ="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/{indexlibid}_n_{n_score}.read_summary.txt.gz"
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
     input: bam = "{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/{indexlibid}_n_{n_score}.bam",
            control=get_control_info
     output: summary ="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/{indexlibid}_n_{n_score}.read_summary.txt.gz"
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
    input: bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/{indexlibid}_n_{n_score}.bam",
    output: stats="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam_stats/{indexlibid}_n_{n_score}.L30-1000_MQ0.3p_dinucleotide_refbase_composition.txt ",
    params: bam="../{indexlibid}_n_{n_score}.bam",
            summary="{indexlibid}_n_{n_score}.summary_damage.txt"
    shell: """
            cd $(dirname {output.plots})
            /home/mmeyer/perlscripts/solexa/analysis/substitution_patterns.pl {params.bam} 
            """

#deam summaries files
#rule deam_stats_split:
 #   input: bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/split_kraken/deam/{indexlibid}_n_{n_score}.k_{kraken_group}.bam",
  #  output: plots="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/split_kraken/deam_stats/{indexlibid}_n_{n_score}.k_{kraken_group}.L30-1000_MQ0.plots.pdf",
   #         summary="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/split_kraken/deam_stats/{indexlibid}_n_{n_score}.k_{kraken_group}.summary_damage.txt"
    #params: bam="../{indexlibid}_n_{n_score}.k_{kraken_group}.bam",
     #       summary="{indexlibid}_n_{n_score}.k_{kraken_group}.summary_damage.txt"
   # shell: """
    #        cd $(dirname {output.plots})
     #       /home/mmeyer/perlscripts/solexa/analysis/substitution_patterns.pl {params.bam} 
      #      /home/mmeyer/perlscripts/solexa/analysis/summarize_CT_frequencies.pl -screen  > {params.summary}
       #     """


#new deam summaries files
rule new_deam_stats_split:
    input: bam="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/{indexlibid}_n_{n_score}.bam",
    output: summary_new="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam_stats/{indexlibid}_n_{n_score}_new.summary_damage.txt"
    params:  summary_new="{indexlibid}_n_{n_score}_new.summary_damage.txt"
    shell: """
            /home/mmeyer/perlscripts/solexa/analysis/quick_substitutions.pl {input.bam} > {output.summary_new}
            """


#rule deam_summaries:
 #   input: plots="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/deam_stats/n_{n_score}/{indexlibid}_n_{n_score}.L30-1000_MQ0.plots.pdf",
  #  output: summary="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/deam_stats/n_{n_score}/{indexlibid}_n_{n_score}.summary_damage.txt"
   # params:"{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam/deam_stats/n_{n_score}/"
    #shell: """    
     #       cd {params}
      #      /home/mmeyer/perlscripts/solexa/analysis/summarize_CT_frequencies.pl -screen 
       #    """

#rule annotate_deam_summary:
 #   input: damage="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam_stats/{indexlibid}_n_{n_score}.summary_damage.txt"
  #  output: damage_annotated="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/deam_stats/{indexlibid}_n_{n_score}.anno_damage.txt",
   # shell:
    #    """
     #   echo -e "IndexLibID\tProbeset\t5'CT\t3'CT\t5'#refC\t3'#refC\t5'CT_95CI\t3'CT_95CI" > {output.damage_annotated}
      #  sed 1d {input.damage} | tr '-' '\t' | cut -f 3- | awk  'BEGIN {{OFS="\t"}} {{print "{wildcards.indexlibid}", "{wildcards.probeset}", $0}}' > {output.damage_annotated}
       # """

#create plots kraken, coverage, SNPs, %primates
rule plot_per_lib:
    input:
        cov="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/{indexlibid}_n_{n_score}.cov",
        kraken="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/{indexlibid}_n_{n_score}.byread",
        primate="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/{indexlibid}_n_{n_score}.read_summary.txt.gz",
        burden=get_bed,
    output:
        "{project}_summary/{indexlibid}/{probeset}/{score_b}_filter/{indexlibid}_n_{n_score}_coverage_plot.pdf"
    params:
        output_dir="{project}_summary/{indexlibid}/{probeset}/{score_b}_filter/",
         n_score=score_n
    run:
        shell(f"""
			Rscript scripts/primates_coverage_kraken_SNPs_burden_plot.R {score_b} {input.burden} {input.cov} {input.kraken} {input.primate} {params.output_dir} {params.n_score}
		""")


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
        cov="{project}/mappedbams/{indexlibid}/{probeset}/rmdupL35MQ25/target/{score_b}_filter/{indexlibid}_n_{n_score}.cov",
        burden=get_bed,
    output:
        "{project}_summary/{indexlibid}/{probeset}/{score_b}_filter/{indexlibid}_n_{n_score}_SNPs_number_cov_combined.pdf"
    params:
        output_dir="{project}_summary/{indexlibid}/{probeset}/{score_b}_filter/",
        n_score=score_n
    run:
        shell(f"""
			Rscript scripts/combined_snps_cov.R {input.cov} {input.burden} {params.output_dir}  {params.n_score}
		""")

#create deam and summary table for all libraries analyzed in the pipeline
#rule summary_deam_all:
 #   input:
  #      dir_deam= "{project}/",
   #     dir_summary="{project_summary}/"
    #output:
     #   summary_all="{project_summary}/all_summaries.txt",
      #  deam_all="{project_summary}/all_deam.txt"
   # params:
    #    b
     #   n_score=score_n
   # run:
    #    shell(f"""
	#		Rscript scripts/save_all_summaries_pipeline.R {input.dir_deam} {input.dir_summary} {output.summary_all} {output.deam_all}
	#	""")


