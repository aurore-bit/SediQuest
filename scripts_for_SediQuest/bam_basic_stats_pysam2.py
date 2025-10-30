import sys, argparse, csv, gzip, pysam
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--bam", help="hg19 mapped sam file")
parser.add_argument("--control", help="control file with list of sites of interest")
parser.add_argument("--control-lim", help="only read lim lines of control file", type=int, default=None)
parser.add_argument("--control-header", help="header for control file")
parser.add_argument("--tags", help="extra fields to add", nargs="+", default=[])
parser.add_argument("--tags-fill", help="content of extra fields to add", nargs="+", default=[])
parser.add_argument("--debug", action='store_true')
parser.add_argument("--ngs-bams", action='store_true', help='Different processing for bams aligned with NextGenMap')
args = parser.parse_args()


print('tags', args.tags, file=sys.stderr)


args.control_header = open(args.control_header, 'rt').readline().rstrip().split('\t')

if args.control.endswith('.gz'):
    args.control = csv.DictReader(gzip.open(args.control, 'rt'), delimiter='\t', fieldnames=args.control_header)
else:
    args.control = csv.DictReader(open(args.control, 'rt'), delimiter='\t', fieldnames=args.control_header)
    pass

samfile = pysam.AlignmentFile(args.bam, "rb")


def get_edit_dist(hread, dist_method, q_pos):
    hread_pairs = hread.get_aligned_pairs(with_seq=True)

    if args.debug_dist:
        print(hread_pairs)
        print(hread.seq + ' <- read')
        print(''.join((p[2] if p[2] != None else '-') for p in hread_pairs if p[0] != None) + ' <- ref')
        print(''.join('.' if p[0] != q_pos else hread.seq[p[0]] for p in hread_pairs if p[0] != None) + ' <- target site')
        pass
    
    ## always count read gaps
    read_gaps = sum(1 for p in hread_pairs if p[0] == None)

    ## remove read gaps
    hread_pairs = [p for p in hread_pairs if p[0] != None]

    if dist_method == 1:
        ## count every difference (gap + snp, even at ends and at target site)
        snps = sum(1 for p in hread_pairs if p[2] != hread.seq[p[0]])
    elif dist_method == 2:
        ## count every difference EXCEPT at target site (gap + snp, even at ends)
        snps = sum(1 for p in hread_pairs if p[0] != q_pos and p[2] != hread.seq[p[0]])
    elif dist_method == 3:
        ## count differences EXCEPT at target site and within 3bp of end of read (gap + snp)
        snps = sum(1 for p in hread_pairs[3:-3] if p[0] != q_pos and p[2] != hread.seq[p[0]])
        pass

    if args.debug_dist:
        print('found %d edits' % (snps + read_gaps))
        pass
    
    return snps + read_gaps



def get_end_deam(hread):
    

    # list of tuples with (query position, ref position, ref base)
    # query is None if there's a del in the read (i.e., we skip bases in the ref while staying put on the read)
    # ref fields are None if there's a del in the ref (i.e., we skip bases in the read while staying put on the ref)
    hread_pairs = hread.get_aligned_pairs(with_seq=True)

    # remove ref sites with no corresponding read site
    hread_pairs = [p for p in hread_pairs if p[0] != None]

    seqlen = len(hread_pairs)
    if seqlen < 6: return ('tooshort', '.')

    if seqlen != len(hread.query_sequence):
        print('query sequence length error when looking at deamination??')
        print(len(hread.query_sequence), hread.query_sequence)
        print(seqlen, hread_pairs)
        return ('lenerror', 'lenerror')

    deam_stats = ['.'] * 6
    for didx,sidx in enumerate((0,1,2,seqlen-3,seqlen-2,seqlen-1)):
        read_base = hread.query_sequence[sidx]
        ref_base = hread_pairs[sidx][2]
        ## sometimes there is no ref base
        if ref_base == None: continue
        if not hread.is_reverse and ref_base.upper() == 'C': deam_stats[didx] = 'CG'
        if hread.is_reverse and ref_base.upper() == 'G': deam_stats[5-didx] = 'CG'
        if read_base != ref_base:
            if args.debug: print( hread.query_sequence)
            if args.debug: print( ' ' * sidx + '*' + '     site?', 'ENDMISMATCH', sidx, read_base, hread_pairs[sidx], 'rev' if hread.is_reverse else 'forward')
            if not hread.is_reverse and ref_base.upper() == 'C' and read_base == 'T':
                deam_stats[didx] = 'deam'
                pass
            if hread.is_reverse and ref_base.upper() == 'G' and read_base == 'A':
                deam_stats[5-didx] = 'deam'
                pass
            # if (hread.is_reverse and ref_base.upper() == 'G' and read_base == 'A') or \
            #    (not hread.is_reverse and ref_base.upper() == 'C' and read_base == 'T'):
            #     if args.debug: print( ' ' * sidx + '*' + '     site?', 'DEAM-ENDMISMATCH', sidx, read_base, hread_pairs[sidx], hread.is_reverse)
            #     if sidx <= 2: d5 = 'deam'
            #     if sidx >= 3: d3 = 'deam'
            #     pass
            pass
        pass
    
    ## return (d5, d3, 'T' if d5 == 'deam' or d3 == 'deam' else 'F')
    return (deam_stats, 'deam' in deam_stats, 'CG' in deam_stats or 'deam' in deam_stats)

    # d5, d3 = ('.', '.')
    # d5cg, d3cg = ('F', 'F')
    # for sidx in (0,1,2,seqlen-3,seqlen-2,seqlen-1):
    #     read_base = hread.query_sequence[sidx]
    #     ref_base = hread_pairs[sidx][2]
    #     ## sometimes there is no ref base
    #     if ref_base == None: continue
    #     if read_base != ref_base:
    #         if args.debug: print( hread.query_sequence)
    #         if args.debug: print( ' ' * sidx + '*' + '     site?', 'ENDMISMATCH', sidx, read_base, hread_pairs[sidx], hread.is_reverse)
    #         if not hread.is_reverse and ref_base.upper() == 'C':
    #             if sidx <= 2: d5cg = 'T'
    #             if sidx >= 3: d3cg = 'T'
    #             if read_base == 'T':
    #                 if sidx <= 2: d5 = 'deam'
    #                 if sidx >= 3: d3 = 'deam'
    #                 pass
    #             pass
    #         if hread.is_reverse and ref_base.upper() == 'G':
    #             if sidx <= 2: d3cg = 'T'
    #             if sidx >= 3: d5cg = 'T'
    #             if read_base == 'A':
    #                 if sidx <= 2: d3 = 'deam'
    #                 if sidx >= 3: d5 = 'deam'
    #                 pass
    #             pass
    #         # if (hread.is_reverse and ref_base.upper() == 'G' and read_base == 'A') or \
    #         #    (not hread.is_reverse and ref_base.upper() == 'C' and read_base == 'T'):
    #         #     if args.debug: print( ' ' * sidx + '*' + '     site?', 'DEAM-ENDMISMATCH', sidx, read_base, hread_pairs[sidx], hread.is_reverse)
    #         #     if sidx <= 2: d5 = 'deam'
    #         #     if sidx >= 3: d3 = 'deam'
    #         #     pass
    #         pass
    #     pass
    
    # ## return (d5, d3, 'T' if d5 == 'deam' or d3 == 'deam' else 'F')
    # return (d5, d3, d5cg, d3cg, 'T' if d5 == 'deam' or d3 == 'deam' else 'F')



first_line = True




####################################
## go through every line in the control file

control = defaultdict(lambda : defaultdict(list))

print('reading control file..', file=sys.stderr)
for i,row in enumerate(args.control):
    control[row['chrom']][int(row['pos'])] += [row]
    if not args.control_lim == None and i > args.control_lim: break
    # print(row['pos'])
    # if int(row['pos']) > 100000000: break
    pass

print('processing bam file..', file=sys.stderr)

## added nofilter on April 28 2023, to fix issue where our NextSeq runs were getting 512 / 516 failed flags
for pileupcolumn in samfile.pileup(min_base_quality = 0, stepper = 'nofilter'):

    #print(pileupcolumn.pos, pileupcolumn.reference_pos, pileupcolumn.reference_name)
    #print(control[pileupcolumn.reference_name].keys())

    #print('hey', pileupcolumn.reference_pos, pileupcolumn.reference_pos+1 in control[pileupcolumn.reference_name], pileupcolumn.reference_pos in control[pileupcolumn.reference_name])
    
    if pileupcolumn.reference_pos+1 not in control[pileupcolumn.reference_name]: continue
    # if pileupcolumn.reference_pos > 100000000: break

    # if pileupcolumn.reference_pos >= 6615103 and pileupcolumn.reference_pos <= 6615110:
    #     for pileupread in pileupcolumn.pileups:
    #         print(pileupread)
    #         ## skip this read if it has a deletion at this site, or skips the ref (insertion?)
    #         print('isdel/isrefskip', pileupread.is_del, pileupread.is_refskip)
    #         if args.debug: print(pileupread.query_position)
    #         if args.debug: print(pileupread.alignment.query_sequence)
    #         if args.debug: print(' '*pileupread.query_position + '*')
    #         if args.debug: print(pileupread.alignment.get_reference_sequence())
    #         if args.debug: print ('\tbase in read %s = %s' %
    #                (pileupread.alignment.query_name,
    #                 pileupread.alignment.query_sequence[pileupread.query_position]))
    #         pass
    #     pass

    for row in control[pileupcolumn.reference_name][pileupcolumn.reference_pos+1]:

        ####################################
        ## get the pileup column at this position - shoooould only be one, but it still gives back an iterator
        pos = int(row['pos'])
        
        if args.debug: print('\n~~~~~~~~~~~~~~~~~~~~~', row['chrom'], row['pos'], row['ref'], row['a1'], row['a2'])
        if args.debug: print ("coverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))

        ####################################
        ## sample all reads from the bam file that overlap this position

        for pileupread in pileupcolumn.pileups:
            
            ## skip this read if it has a deletion at this site, or skips the ref (insertion?)
            # print('isdel/isrefskip', pileupread.is_del, pileupread.is_refskip)
            if pileupread.is_del or pileupread.is_refskip: continue

            read = pileupread.alignment.query_name
            if args.debug: print(pileupread.query_position)
            if args.debug: print(pileupread.alignment.query_sequence, '<- read')
            if args.debug: print(' '*pileupread.query_position + '*')
            if args.debug: print(pileupread.alignment.get_reference_sequence(), '<- ref')
            if args.debug: print ('\tbase in read %s = %s' %
                   (pileupread.alignment.query_name,
                    pileupread.alignment.query_sequence[pileupread.query_position]))

            # if args.debug and pileupread.alignment.query_sequence[pileupread.query_position] != pileupread.alignment.get_reference_sequence()[pileupread.query_position]:
            #    print('QUERY SITE MISSMATCH - only works if no idel!')

            ## the base in the read at the query position
            read_base = pileupread.alignment.query_sequence[pileupread.query_position]
            read_len = len(pileupread.alignment.query_sequence)
            read_start = pileupread.alignment.reference_start
            read_end = pileupread.alignment.reference_end

            ## check for deamination
            #deam5,deam3,d5cg,d3cg,deam53 = get_end_deam(pileupread.alignment)
            deam_stats,deam53,cg53 = get_end_deam(pileupread.alignment)
            deam_tags = ['deam5.%d' % i for i in (1,2,3)] + ['deam3.%d' % i for i in (3,2,1)]

            ## print something if we don't match either expected state, but still process the read
            if read_base != row['a1'] and read_base != row['a2'] and args.debug: print("BAD READ - DOESN'T MATCH EXPECTED STATES")

            if len(pileupread.alignment.get_reference_sequence()) != len(pileupread.alignment.query_sequence) and \
               args.debug: print('LENGTH MISMATCH REF QUERY')
            if len(pileupread.alignment.get_aligned_pairs(with_seq=True)) != len(pileupread.alignment.query_sequence) and \
               args.debug: print('LENGTH MISMATCH QUERY/PAIRS', pileupread.alignment.get_aligned_pairs(with_seq=True))
            
            if pileupread.alignment.query_sequence[pileupread.query_position] != read_base:
                print("ERROR query position doesn't match read base?", pileupread.alignment.query_sequence[pileupread.query_position],  read_base)
                pass

            read_masked = list(pileupread.alignment.query_sequence)
            read_masked[pileupread.query_position] = 'N'
            read_masked = ''.join(read_masked)

            base_qual = pileupread.alignment.query_qualities[pileupread.query_position]
            
            # print('PAIRS', pileupread.alignment.get_aligned_pairs(with_seq=True))
            indel_bases = sum(None in column for column in pileupread.alignment.get_aligned_pairs(with_seq=True))
            
            if first_line:
                print('read_id', 'read_base', 'read_len', 'read_start', 'read_end', 'strand', *deam_tags, 'deam53', 'cg53',
                      *args.control_header, *args.tags, 'read_seq', 'read_seq.x', 'base_qual', 'read_target_pos', 'read_has_indel', 'indel_bases', sep='\t')
                first_line = False
                pass
            print(read, read_base, read_len, read_start, read_end, 'B' if pileupread.alignment.is_reverse else 'T', *deam_stats, deam53, cg53,
                  *[row[x] for x in args.control_header], *args.tags_fill,
                  # print the read sequence, plus the read sequence with an N replacing the target base
                  pileupread.alignment.query_sequence, read_masked,
                  base_qual, pileupread.query_position,
                  'T' if indel_bases > 0 else 'F', indel_bases,
                  sep='\t')
            pass
            
        pass
    
    pass


# for line in args.sam:
#     line = line.rstrip().split('\t')
#     read = line[0]

#     print(line)
#     for spc in args.species:
#         if read in sdict[spc]: print(spc, sdict[spc][read])
#         pass

#     hum_dist = get_min_dist(read, ('homo_sapiens',))
#     ape_dist = get_min_dist(read, apes)
#     monkey_dist = get_min_dist(read, monkeys)
#     quadraped_dist = get_min_dist(read, quadrapeds)
#     carnivore_dist = get_min_dist(read, carnivorus)
#     hopper_dist = get_min_dist(read, hoppers)

#     spc_groups = ['hum', 'ape', 'monkey', 'glires', 'carn', 'ung']
#     spc_group_dists = (hum_dist, ape_dist, monkey_dist, hopper_dist, carnivore_dist, quadraped_dist)
#     min_spc_groups = [spc_groups[i] for i in range(6) if spc_group_dists[i] == min(spc_group_dists)]
#     if min(spc_group_dists) == args.global_min:
#         min_spc_groups_label = 'none'
#     else:
#         min_spc_groups_label = '.'.join(sorted(min_spc_groups))
        
#     print('REPORT', read, hum_dist, ape_dist, monkey_dist, hopper_dist, carnivore_dist, quadraped_dist, ' '.join(min_spc_groups), min_spc_groups_label, sep='\t')
#     pass

# for read in sdict['homo_sapiens']:

#     print()
#     print(read)
#     for spc in args.species:
#         if read in sdict[spc]: print(spc, sdict[spc][read])

        
