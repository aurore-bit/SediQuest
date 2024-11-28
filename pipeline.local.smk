################################################################################
# Draft nuclear sediment pipeline
#
# Aurore GALTIER, 25/11/2024
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

# Define summary folder
summary = config["summary"]

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

##############################################
# map_all Rule
##############################################

mapped_bams = expand('{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.bam',
                     indexlibid=INDEXLIBID,	
		             probeset=indexlibid_df.loc[INDEXLIBID, "probeset"],
                     project=project)
  
rule map_all:
    input:
        bam_files = mapped_bams  # Use the list of files here
    run:
        print("Hey, mapping is done")
        pass

rule map:
    input:
        bam = get_bam,  
    output:
        mapped = "{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.bam",
        sorted = "{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.sorted.bam"
    threads: 8
    params: ref = get_ref
    shell:
        """
        echo 'Hey! Mapping BAM files...'
        bwa bam2bam -t {threads} -g {params.ref} -n 0.01 -o 2 -l 16500 --only-aligned {input.bam}  | samtools sort -@30 -o {output.mapped} 
        /home/visagie/.local/bin/samtools sort --threads {threads} -o {output.sorted} {output.mapped}
        """

##############################################
#processing
##############################################

rule processing:
    input: 
        expand('{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/deam/{indexlibid}.uniq.L35MQ25.deam53x3.bam',
               indexlibid=INDEXLIBID,
               probeset=indexlibid_df.loc[INDEXLIBID, "probeset"],
	       project=project)
    run:
        print('hey! Running target filtering, deamination, quality filtering and duplicate removal :)')
        pass


rule filter_bam_by_control_sites:
    input:
        control_sites = get_bed,
        bam = "{project}/mappedbams/{indexlibid}/{probeset}/{indexlibid}.bam"
    output:
        sites_bam="{project}/mappedbams/{indexlibid}/{probeset}/target/{indexlibid}.bam",
    threads: 1
    #conda: "envs/processing.yaml"
    shell: """
        bedtools intersect -a {input.bam} -b {input.control_sites} > {output.sites_bam} 
    """

rule uniq_q25_l35:
    input:
        bam = "{project}/mappedbams/{indexlibid}/{probeset}/target/{indexlibid}.bam",
    output:
        bam="{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/{indexlibid}.uniq.L35MQ25.bam",
    shell:
            """
            bam-rmdup -o {output.bam} -q 25 -l 35 {input.bam}
            """

rule deam_filter:
    input: bam="{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/{indexlibid}.uniq.L35MQ25.bam"
    output: bam="{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/deam/{indexlibid}.uniq.L35MQ25.deam53x3.bam",
        stats="{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/deam/{indexlibid}.uniq.L35MQ25.deam53x3.stats"
    #conda: "envs/processing.yaml"
    shell: """
    /home/mmeyer/perlscripts/solexa/analysis/filterBAM.pl  -p5 0,1,2 -p3 0,-1,-2 -suffix deam53x3 {input.bam} &> {output.stats}
    """


##############################################
#kraken
##############################################
def get_bam_files():
    bam_files = []
    for indexlibid in INDEXLIBID:
        for probeset in indexlibid_df.loc[indexlibid, "probeset"]:
            paths = [
                f"{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/deam/*.bam",
                f"{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/*.bam",
                f"{project}/mappedbams/{indexlibid}/{probeset}/target/*.bam",
                f"{project}/mappedbams/{indexlibid}/{probeset}/*.bam"
            ]
            for path in paths:
                bam_files.extend(glob.glob(path))
    return bam_files

bam_files_list = get_bam_files()

base_names = [os.path.basename(bam).replace('.bam', '') for bam in bam_files_list]

rule kraken:
    input:
        expand("{base}.kraken_spc", base=base_names)
    run:
        print('hey! Running Kraken')
        pass

rule bam_to_fasta:
    input:
        split=get_bam,
        bam="{base}/{base_name}.bam"
    output:
        bam="{base}/{base_name}.bam.fa.gz",
        split="{base}/kraken_split/{base_name}.bam.fa.gz"
    threads: 1
   # conda: "envs/kraken.yaml"
    shell: """
     bam2fastx -a -A -Q {input.bam} | /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/lenfilter.pl 35 | gzip -c > {output.bam},
     bam2fastx -a -A -Q {input.split} | /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/lenfilter.pl 35 | gzip -c > {output.split}
    """

rule generic_kraken:
    input: fa="{base}.fa.gz"
    output: kraken="{base}.kraken",
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
    input: kraken="{base}.kraken"
    output: phylo="{base}.kraken_phylo"
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

#rule generic_kraken_translate:
 #   input: kraken="{base}.kraken"
  #  output: translate="{base}.translate"
   # threads: 1
   # conda: "envs/kraken.yaml"
   # shell: """

    #if [ $(hostname) != "bionc13" ] ; then 
     #    echo KRAKEN MUST BE RUN ON BIONC13. SLEEPING 30 SECONDS.
      #   sleep 30
       #  stop
   # fi

    #echo 'making report..'
    ## translate file has info on each read (but not in an easy-to-digest way)
    ## phylo file has a summary for each level in the taxonomy
   # python3 /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/kraken_report.py --db /mnt/ramdisk/refseqReleaseKraken {input.kraken} --translate {output.translate} > /dev/null

    #"""

#rule generic_kraken_byread:
 #   input: translate="{base}.translate"
  #  output: byread="{base}.byread"
   # threads: 1
   # conda: "envs/kraken.yaml"
    #shell: """

     #echo 'parsing translate file..'
    # python3 /mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/bin/parse_kraken_translate_file.py --translate-file {input.translate} > {output.byread}

   # """

### summarize the major species groups in the kraken phylogeny file

rule generic_kraken_spc_summary:
    input: spc_krak="{base}.kraken_phylo",
        spc_groups="/mnt/expressions/benjamin_vernot/soil_capture_2017/process_sequencing/data/major_species_groups.txt"
    output: spc_summary="{base}.kraken_spc"
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


##############################################
#summary
##############################################
rule summaries:
    input:
        expand("{summary}/{indexlibid}/{indexlibid}.kraken_mapped.pdf",
    summary=summary,
    indexlibid=INDEXLIBID),
        "{project}/mappedbams/{indexlibid}/{probeset}/{seqrun}/split_exact/map_{ref}/target/rmdupL35MQ25/{indexlibid}.pipeline_summary.txt",
        "{project}/mappedbams/{indexlibid}/{probeset}/{seqrun}/split_exact/map_{ref}/target/rmdupL35MQ25/deam/deam_stats/{indexlibid}.uniq.L35MQ25.deam53x3.anno_damage.txt",
    run:
        print('hey! Creating summaries')
        pass



rule plot:
    input: "{project}/mappedbams/{indexlibid}/{probeset}/target/rmdupL35MQ25/deam/{indexlibid}.uniq.L35MQ25.deam53x3.bam"
    output: plot="{summary}/{wildcards.indexlibid}/{wildcards.probeset}/{wildcards.seqrun}/coverage_plot.png"
    run:
        shell(f"""
			cd {summary_dir}
			mkdir -p {wildcards.indexlibid}
			samtools mpileup -B {project}/mappedbams/{wildcards.indexlibid}/{wildcards.probeset}/{wildcards.seqrun}/target/rmdupL35MQ25/{wildcards.indexlibid}.uniq.L35MQ25.bam > {project}/mappedbams/{wildcards.indexlibid}/{wildcards.probeset}/{wildcards.seqrun}/target/rmdupL35MQ25/{wildcards.indexlibid}.uniq.L35MQ25.bam.cov
			Rscript scripts/%primates_coverage_kraken_SNPs_burden_plot.R
		""")

rule pipeline_summary:
    input: map_bam="{project}/mappedbams/{indexlibid}/{probeset}/{seqrun}/{indexlibid}.bam",
        target_bam="{project}/mappedbams/{indexlibid}/{probeset}/{seqrun}/target/{indexlibid}.bam",
        rmdup_bam="{project}/mappedbams/{indexlibid}/{probeset}/{seqrun}/target/rmdupL35MQ25/{indexlibid}.uniq.L35MQ25.bam",
        rmdup_stats="{project}/mappedbams/{indexlibid}/{probeset}/{seqrun}/target/rmdupL35MQ25/{indexlibid}.summary_stats_L35MQ25.txt"
    output: summary_annotated="{project}/mappedbams/{indexlibid}/{probeset}/{seqrun}/target/rmdupL35MQ25/{indexlibid}.pipeline_summary.txt"
    shell:
        """

        splitbam="{wildcards.project}/splitbams/{wildcards.seqrun}/{wildcards.indexlibid}.bam"

        if [ -e $splitbam ] ; then
          sc=$(samtools view -c $splitbam)
        else
          sc="NA"
        fi

        tmp={output.summary_annotated}

        echo IndexLibID probeset ref {params.tags} split mapped \
        $(head -n1 {input.rmdup_stats} | cut -f2- | sed 's/&//g' | sed 's/%/pct_/g' | sed 's/raw/target/g') \
        primate | tr ' ' '\t' > $tmp

        echo {wildcards.indexlibid} {wildcards.probeset} {wildcards.ref} {params.tag_vals} $sc \
        $(samtools view -c {input.map_bam}) \
        $(tail -n1 {input.rmdup_stats} | cut -f2-) \
        $(samtools view -c {input.primate_bam}) \
        | tr ' ' '\t' >> $tmp

        """


#deam summaries files
rule deam_summaries:
    input: bam="{project}/mappedbams/{indexlibid}/{probeset}/{seqrun}/target/rmdupL35MQ25/deam/{indexlibid}.uniq.L35MQ25.deam53x3.bam"
    output: plots="{project}/mappedbams/{indexlibid}/{probeset}/{seqrun}/target/rmdupL35MQ25/deam/deam_stats/{indexlibid}.uniq.L35MQ25.deam53x3.L30-1000_MQ0.plots.pdf",
        summary="{project}/mappedbams/{indexlibid}/{probeset}/{seqrun}/target/rmdupL35MQ25/deam/deam_stats/{indexlibid}.uniq.L35MQ25.deam53x3.summary_damage.txt"
    shell: """
           /home/mmeyer/perlscripts/solexa/analysis/substitution_patterns.pl {input.bam} > {output.plots}
	   /home/mmeyer/perlscripts/solexa/analysis/summarize_CT_frequencies.pl -screen > {output.summary]
           """

rule annotate_deam_summary:
    input: damage="{project}/mappedbams/{indexlibid}/{probeset}/{seqrun}/target/rmdupL35MQ25/deam/deam_stats/{indexlibid}.uniq.L35MQ25.deam53x3.summary_damage.txt"
    output: damage_annotated="{project}/mappedbams/{indexlibid}/{probeset}/{seqrun}/target/rmdupL35MQ25/deam/deam_stats/{indexlibid}.uniq.L35MQ25.deam53x3.anno_damage.txt",
    shell:
        """
        sed 1d {input.damage} | tr '-' '\t' | cut -f 3- | awk  'BEGIN {{OFS="\t"}} {{print "{wildcards.indexlibid}", "{wildcards.probeset}", "{wildcards.seqrun}", "{wildcards.ref}", $0}}' > {output.damage_annotated}
        """


