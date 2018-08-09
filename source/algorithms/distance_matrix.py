#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/distance_matrix.py

# Import non-standard packages
import regex

def build_score(match=1, transition=-1, transversion=-1.2, gap_open=-4, gap_extend=-2):
    characters = ['-', 'A', 'C', 'G', 'T']
    scoring = {}
    for c1 in characters:
        for c2 in characters:
            if ((c1 == '-') or (c2 == '-')):
                scoring[(c1, c2, 'open')] = gap_open
                scoring[(c1, c2, 'extend')] = gap_extend
            elif (c1 == c2):
                scoring[(c1, c2)] = match
            elif (((c1 in ['A', 'G']) and (c2 in ['A', 'G'])) or ((c1 in ['C', 'T']) and (c2 in ['C', 'T']))):
                scoring[(c1, c2)] = transition
            else:
                scoring[(c1, c2)] = transversion
    return scoring

# Scoring matrix of log-likelihood ratios with weights for
# transitions/transversions, gap penalties, etc
# SCORES = {
#     # Matches
#     ('A', 'A'): 1,
#     ('C', 'C'): 1,
#     ('G', 'G'): 1,
#     ('T', 'T'): 1,
#     # Transitions (Purines: A<->G, Pyrimidines: C<->T)
#     ('A', 'G'): -1,
#     ('G', 'A'): -1,
#     ('C', 'T'): -1,
#     ('T', 'C'): -1,
#     # Transversions
#     ('A', 'T'): -1.2, # Weak
#     ('T', 'A'): -1.2, # Weak
#     ('C', 'G'): -1.2, # Strong
#     ('G', 'C'): -1.2, # Strong
#     ('A', 'C'): -1.2, # Amino
#     ('C', 'A'): -1.2, # Amino
#     ('G', 'T'): -1.2, # Keto
#     ('T', 'G'): -1.2, # Keto
#     # Insertions
#     ('-', 'A', 'open'): -4,
#     ('-', 'C', 'open'): -4,
#     ('-', 'G', 'open'): -4,
#     ('-', 'T', 'open'): -4,
#     ('-', 'A', 'extend'): -2,
#     ('-', 'C', 'extend'): -2,
#     ('-', 'G', 'extend'): -2,
#     ('-', 'T', 'extend'): -2,
#     # Deletions
#     ('A', '-', 'open'): -4,
#     ('C', '-', 'open'): -4,
#     ('G', '-', 'open'): -4,
#     ('T', '-', 'open'): -4,
#     ('A', '-', 'extend'): -2,
#     ('C', '-', 'extend'): -2,
#     ('G', '-', 'extend'): -2,
#     ('T', '-', 'extend'): -2,
#     # Gaps should never be aligned to each other
#     ('-', '-', 'open'): -4,
#     ('-', '-', 'extend'): -2,
# }
    
def build_iupac_score(single_scoring):
    """
    Extrapolate input scoring matrix to handle IUPAC ambiguity codes.
    Input argument should be the scoring dict, which can be generated
    with the 'build_score()' function.
    """
    #single_scoring = build_score()
    iupac_scoring = {}
    
    iupac = {
        '-': ['-'],
        'A': ['A'],
        'C': ['C'],
        'G': ['G'],
        'T': ['T'],
        'R': ['A', 'G'],
        'Y': ['C', 'T'],
        'M': ['A', 'C'],
        'K': ['G', 'T'],
        'W': ['A', 'T'],
        'S': ['C', 'G'],
        'B': ['C', 'G', 'T'],
        'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'],
        'V': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'T'],
    }
    
    for c1 in iupac.keys():
        for c2 in iupac.keys():
            scores = []
            go_scores = []
            ge_scores = []
            for i1 in iupac[c1]:
                for i2 in iupac[c2]:
                    if ((i1 == '-') or (i2 == '-')):
                        go_scores.append(single_scoring[(i1, i2, 'open')])
                        ge_scores.append(single_scoring[(i1, i2, 'extend')])
                    else:
                        scores.append(single_scoring[(i1, i2)])
            
            if (len(scores) > 0):
                iupac_scoring[(c1, c2)] = sum(scores)/len(scores)
            
            go_concat = scores + go_scores
            if (len(go_scores) > 0):
                iupac_scoring[(c1, c2, 'open')] = sum(go_concat)/len(go_concat)
            
            ge_concat = scores + ge_scores
            if (len(ge_scores) > 0):
                iupac_scoring[(c1, c2, 'extend')] = sum(ge_concat)/len(ge_concat)
    return iupac_scoring

def load_scores(file_path, sep="\t"):
    """Load the alignment scores"""
    scoring = {}
    with open(file_path, 'r') as flo:
        for line in flo:
            line = line.rstrip()
            if (len(line) > 0):
                m = regex.match(r'^\s*#', line)
                if not m:
                    sline = line.split(sep)
                    if (sline[2] == ''):
                        scoring[(sline[0], sline[1])] = float(sline[3])
                    else:
                        scoring[(sline[0], sline[1], sline[2])] = float(sline[3])
    return scoring
