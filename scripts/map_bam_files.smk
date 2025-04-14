import pandas as pd
import os

# Read the CSV files using pandas
indexlibid_df = pd.read_csv("config/nuclear_data_subset_2025.csv", sep=";")
probeset_df = pd.read_csv("config/probeset.csv", sep="\t")

INDEXLIBID_TO_PROBESETS = indexlibid_df.groupby('library')['probeset'].apply(list).to_dict()
PROBESET_TO_REF = probeset_df.set_index('probeset')['path_to_ref'].to_dict()

INDEXLIBID = indexlibid_df['library'].unique().tolist()

def get_ref(wildcards):
    return PROBESET_TO_REF[wildcards.probeset]

rule all:
    input:
        [f"/mnt/expressions/Aurore/sediment_pipeline_v0/overview_output/split/{indexlibid}/{probeset}/{indexlibid}.sorted.bam"
         for indexlibid, probesets in INDEXLIBID_TO_PROBESETS.items()
         for probeset in probesets]
    run:
        print('hey! Mapping bam files check mappr finger crossed lolilol :)')


rule sort_bam:
    input:
        bam="/mnt/expressions/Aurore/sediment_pipeline_v0/overview_output/split/{indexlibid}/{probeset}/{indexlibid}.bam"
    output:
        sorted="/mnt/expressions/Aurore/sediment_pipeline_v0/overview_output/split/{indexlibid}/{probeset}/{indexlibid}.sorted.bam"
    shell:
        "/home/visagie/.local/bin/samtools sort -n {input.bam} -o {output.sorted}"

mapped_bams = [f"overview_output/mappedbams/{indexlibid}/{probeset}/{indexlibid}.bam"
	for indexlibid, probesets in INDEXLIBID_TO_PROBESETS.items()
         for probeset in probesets]


rule map_all:
    input: [f + '.submitted' for f in mapped_bams]
    run:
        print('\n\n')
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print('hey! Submitting mapped bam jobs (if any exist). Run rule "check_mapped" to confirm they are done, or just try running the whole snakemake file.')
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print('\n\n')
        pass



rule submit_map_bam:
    input:
        bam = "overview_output/split/{indexlibid}/{probeset}/{indexlibid}.bam",
	sorted="overview_output/split/{indexlibid}/{probeset}/{indexlibid}.sorted.bam"
#        bam = _get_and_print_splitbam_loc
    output:
        bam = "overview_output/mappedbams/{indexlibid}/{probeset}/{indexlibid}.bam.submitted",
#        out="overview_output/mappedbams/{indexlibid}/{probeset}/{indexlibid}.bam"   
    params:
        job_name="{indexlibid}_{probeset}_mapping",
        ref=get_ref
    run:
        print('hey mapping bam')
        print(input.bam)
        print(output.bam)
        obam = re.sub(r"\.submitted$", "", output.bam)
        # print(obam)
        cmd = f"/mnt/solexa/bin/mappr-cli -z {params.job_name} -a -g {params.ref} -f {obam} {input.bam} -n 0.01 -o 2 -l 16500 --only-aligned"
#        cmd = "/mnt/solexa/bin/mappr-cli -z %s -a -g %s -f %s %s -n 0.01 -o 2 -l 16500 --only-aligned" % (job_name, genome, obam, input.bam)
        print('command:', cmd)
        shell(cmd)
        # touch the bam.submitted flag file
        Path(output.bam).touch()
        pass



rule check_mapped:
#    input: [f + '.submitted' for f in mapped_bams]
    run:
        print('\n\n')
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print('\nhey! Checking submitted mapped bam jobs:')
        not_finished = 0
        not_submitted = 0
        not_fixed = 0
        for i,f in enumerate(mapped_bams):
            my_sub = Path(f + '.submitted')

            #indexlibid, _, _, _, _ = f.split('/')
            my_splitbam = f"overview_output/split/{indexlibid}/{probeset}/{indexlibid}.bam"
            #my_splitbam = '%s/splitbams/%s/%s/%s/split_exact/%s.bam' % (project, seqrun, proc, lane, indexlibid)
            
            if not my_sub.is_file():
                #print(f'Bam file {i} has not been submitted to mappr: {f}.submitted does not exist!')
                print('Bam file %d has not been submitted to mappr: %s.submitted does not exist!' % (i,f))
                not_submitted += 1
                continue
            
            my_bam = Path(f)
            if my_bam.is_file():
                print(os.path.getmtime(my_sub), my_sub)
                print(os.path.getmtime(my_bam), my_bam)
                print(os.path.getmtime(my_splitbam), my_splitbam)
                if not os.path.exists(my_splitbam):
                    print('mapped bam file has no splitbam:', my_sub)
                    print(' - splitbam:', my_splitbam)
                    # print('la', f"{seqrun}/split_exact/{lane}/*/{indexlibid}.bam")
                elif os.path.getmtime(my_bam) < os.path.getmtime(my_splitbam) or os.path.getmtime(my_sub) < os.path.getmtime(my_splitbam):
                    if not fix_times:
                        print('splitbam is modified after mapped file(s): rerun with fix_times=T to fix.')
                        print(' : ', my_sub)
                        print(' : ', my_splitbam)
                        not_fixed += 1
                    else:
                        mtime = min(os.path.getmtime(my_bam), os.path.getmtime(my_sub))
                        print('fixing times:')
                        print(mtime, my_sub)
                        print(os.path.getmtime(my_splitbam), my_splitbam)
                        os.utime(my_splitbam, (mtime-20, mtime-20))
                        # exit(-1)
                        pass
                    pass
                continue
            #print(f'Mapped bam file {f} does not exist!')
            print('Mapped bam file %s does not exist!' % f)
            print(' - splitbam:', my_splitbam)
            not_finished += 1
            pass

        if not_fixed > 0:
            #print(f'\n{not_fixed} mapped bams have modify times after their split bams. Run with fix_times=T to fix! [by changing mod time of split bam]')
            print('\n%d mapped bams have modify times after their split bams. Run with fix_times=T to fix! [by changing mod time of split bam]' % not_fixed)
            pass
        
        if not_submitted > 0:
            #print(f'\n{not_submitted} mappr jobs have not been submitted yet - run map_all rule first!')
            print('\n%d mappr jobs have not been submitted yet - run map_all rule first!' % not_submitted)
            pass
        
        if not_finished > 0:
            #print(f'\n{not_finished} mappr jobs have not completed.')
            print('\n%d mappr jobs have not completed.' % not_finished)
            pass
        
        if not_submitted == 0 and not_finished == 0:
            #print(f'\nAll mappr jobs submitted and finished!')
            print('\nAll mappr jobs submitted and finished!')
            pass
        
        print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        print('\n\n')

        pass




