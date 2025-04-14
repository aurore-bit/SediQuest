
import sys, gzip
import argparse
import random
#AA187
args.sites= "/mnt/expressions/benjamin_vernot/soil_capture_2017/site_categories_for_capture/soil_probe_designs/probes_archaic_hets_n3_b6_CONTROL_BV11-BV11.txt.gz.bam_sorted.bed"
args.fasta = "/mnt/solexa/Genomes/hg19_evan/whole_genome.fa"
args.fasta_ofile ="/mnt/expressions/Aurore/sediment_pipeline_v0/kraken_database/AA187_whole_genome.fa"
args.replace="nonref"


if args.sites_format == 'control':
    args.sites = gzip.open(args.sites, 'rt')
elif args.sites_format in ('bed', 'bed3'):
    args.sites = open(args.sites, 'r')
    pass

args.fasta = open(args.fasta, 'r')
args.fasta_ofile = open(args.fasta_ofile, 'w')

print('reading sites', args.sites, file=sys.stderr)
sites = {}
for i,line in enumerate(args.sites):

    ## control format file has a header, skip it
    if i == 0 and args.sites_format == 'control':
        continue

    line = line.rstrip().split()

    if args.sites_format == 'control':
        chrom, pos, ref, a1, a2 = line[1:6]
    elif args.sites_format == 'bed':
        chrom, pos0, pos, ref, a2 = line[0:5]
        a1 = ref
    elif args.sites_format == 'bed3':
        chrom, pos0, pos, ref, a2, a3 = line[0:6]
        a1 = ref
        pass

    pos = int(pos)
    # print chrom, pos
    if not chrom in sites: sites[chrom] = {}
    if args.sites_format == 'bed3':
        sites[chrom][pos] = (ref,a1,a2,a3)
    else:
        sites[chrom][pos] = (ref,a1,a2)
        pass
    # break
    pass

def modify_fasta(myseq, mychr, header):
    fail = False
    myseq = list(myseq)
    bases = ('A', 'C', 'G', 'T')

    if mychr not in sites:
        print('fasta chromosome', mychr, 'not in list of sites - outputing unmodified', file=sys.stderr)
        return myseq

    site_count = 0
    ref_mismatch_count = 0
    for s in sorted(sites[mychr]):

        site_count += 1

        if s > len(myseq):
            print('FASTA too short for site, skipping rest:', mychr, s, '>', len(myseq), file=sys.stderr)
            # fail = True
            sys.exit(-1)
            break

        r = sites[mychr][s][0]  # ref from sites file
        a1 = sites[mychr][s][1] # allele 1 from sites file
        a2 = sites[mychr][s][2] # allele 2 from sites file
        r_fa = myseq[s-1]       # ref from fasta

        if args.replace == 'third' and args.sites_format == 'bed3':
            choices = [b for b in bases if b != a1 and b != a2]
            a3 = sites[mychr][s][3] # allele 3 from sites file
            r_replace = a3
            if a3 not in choices:
                print( 'HEY! trying to put in the third allele specified in the sites file, but it does not match what I expected.', mychr, s )
                print( 'Expected non-a1/a2 alleles:', choices)
                print( 'a3:', a3 )
                pass

            if a3 not in bases:
                print( 'HEY! trying to put in the third allele specified in the sites file, but it IS NOT A PROPER NUCLEOTIDE.', mychr, s )
                print( 'Expected non-a1/a2 alleles:', choices)
                print( 'a3:', a3 )
                print( 'SKIPPING THIS SITE' )
                continue

        elif args.replace == 'third':
            choices = [b for b in bases if b != a1 and b != a2]
            r_replace = random.sample(choices,1)[0]
            pass
        elif args.replace == 'nonref' and a1 == r_fa:
            # if a2 == ref, then this won't actually modify the ref (but I think this is the correct behavior
            r_replace = a2
        elif args.replace == 'nonref' and a2 == r_fa:
            r_replace = a1
        elif args.replace == 'nonref':
            print( "HEY! trying to replace ref w/ nonref allele, but both alleles are nonref.  leaving as is!", mychr, s, r, a1, a2, r_fa, file=sys.stderr)
            print(s//80*80, s//80*80+80, file=sys.stderr)
            print(''.join(myseq[s//80*80:s//80*80+80]), file=sys.stderr)
            r_replace = r_fa
            pass

        ## THIS DOESN'T SEEM TO DO WHAT I SAY IT DOES?
        if r != r_fa:
            ref_mismatch_count += 1
            print( "HEY! ref doesn't match fasta, replacing with a1 instead:", mychr, s,
                   'fasta_ref=%s; sites ref=%s, a1=%s, a2=%s [%d / %d sites on chr]' % (r_fa, r, a1, a2, ref_mismatch_count, site_count), file=sys.stderr)
            # print(s//80*80, s//80*80+80, file=sys.stderr)
            # print(''.join(myseq[s//80*80:s//80*80+80]), file=sys.stderr)
            # sys.ex
            fail = True
            pass

        ## chromosome, position, ref allele [from control file], a1, a2, the replacement allele
        print( mychr, s, r, a1, a2, r_replace, file=args.mod_ofile)
        myseq[s-1] = r_replace
        pass

    # if fail:
    #     print('At least one site with reference missmatch, exiting.', file=sys.stderr)
    #     sys.exit(-1)
    #     pass

    return myseq
