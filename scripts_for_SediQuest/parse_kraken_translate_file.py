import sys, argparse, csv, gzip, pysam, random

parser = argparse.ArgumentParser()
parser.add_argument("--translate-file", help="output file of fred's translate tool", required=True)
parser.add_argument("--debug", action='store_true')
args = parser.parse_args()




for line in open(args.translate_file, 'rt'):

    if args.debug: print(line)
    
    line = line.strip().split(maxsplit=2)
    if len(line) == 3:
        (read_id, root, lineage) = line
    else:
        continue

    lineage = lineage.split(';')
    if args.debug: print(lineage)

    if args.debug: print('TABLE', '\t'.join(lineage))

    if 'Primates' in lineage: read_cat = 'Primates'
    elif 'Rodentia' in lineage: read_cat = 'Rodentia'
    elif 'Glires' in lineage: read_cat = 'Glires'
    elif 'Caniformia' in lineage: read_cat = 'Caniformia'
    elif 'Feliformia' in lineage: read_cat = 'Feliformia'
    elif 'Carnivora' in lineage: read_cat = 'Carnivora'
    elif 'Ruminantia' in lineage: read_cat = 'Ruminantia'
    elif 'Equidae' in lineage: read_cat = 'Equidae'
    elif 'Afrotheria' in lineage: read_cat = 'Afrotheria'
    elif 'Suina' in lineage: read_cat = 'Suina'
    elif 'Cetacea' in lineage: read_cat = 'Cetacea'
    elif 'Chiroptera' in lineage: read_cat = 'Chiroptera'
    # elif '' in lineage: read_cat = ''
    # elif '' in lineage: read_cat = ''
    # elif '' in lineage: read_cat = ''
    else: read_cat = 'Other'

    keepcharacters = ('.','_')
    final_lin = "".join(c if c.isalnum() or c in keepcharacters else '_' for c in lineage[-1]).rstrip()
            
    print(read_id, read_cat, final_lin, sep='\t')

    # print(read_id)
