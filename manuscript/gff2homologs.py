#!/usr/bin/env python3
# Copyright 2017 Thaddeus D. Seher

# USAGE
#  python3 gff2homologs.py C_albicans_SC5314_A22_current_features.gff > C_albicans_SC5314_A22_current_homologs.txt

import sys
import urllib.parse

new = []
with open(sys.argv[1], 'r') as flo:
    for line in flo:
        if not line.startswith('#'):
            sline = line.rstrip().split('\t')
            if (sline[2] == 'gene'):
                fields = urllib.parse.unquote(sline[8]).split(';')
                gene_tag = None
                for f in fields:
                    if f.startswith('Gene='):
                        gene_tag = f.split('=')[1]
                        break
                else:
                    for f in fields:
                        if f.startswith('Note='):
                            gene_tag = f.split('=')[1].split(' ')[0][1:-1]
                            break
                
                id_tag = None
                for f in fields:
                    if f.startswith('ID='):
                        id_tag = f.split('=')[1]
                        break
                
                if ((len(new) != 0) and (id_tag[:-2] == new[-1][1][:-2])):
                    new[-1].append(id_tag)
                    if (new[-1][0].startswith('orf') and not gene_tag.startswith('orf')):
                        new[-1][0] = gene_tag
                else:
                    new.append([gene_tag, id_tag])

for line in new:
    print("\t".join(line))
