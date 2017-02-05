#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# utils.py

# Import standard packages
import string

def lcs(string1, string2):
    """Find the longest common substring between two strings"""
    import difflib
    
    matcher = difflib.SequenceMatcher(None, string1, string2, True)
    match = matcher.find_longest_match(0, len(string1), 0, len(string2))
    # Match(a=0, b=15, size=9)
    return match
    # print(string1[match.a: match.a + match.size])
    # print(string2[match.b: match.b + match.size])

def rc(seq, kind="dna"):
    """Returns the reverse-complement of a string containing DNA or RNA characters"""
    if (kind == "dna"):
        complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB') # exclude ws, WS
    elif (kind == "rna"):
        complements = string.maketrans('acgurymkbdhvACGTRYMKBDHV', 'ugcayrkmvhdbTGCAYRKMVHDB') # exclude ws, WS
    else:
        raise ValueError("'" + str(kind) + "' is an invalid argument for rc()")
    return seq.translate(complements)[::-1]

def read_fasta_file(args):
    """Load contig sequences from file into dict()
    Primary sequence headers must be unique
    """
    contigs = {}
    with open(args.fasta, 'r') as flo:
        name = None
        for line in flo:
            line = line.rstrip()
            if line.startswith('>'):
                name = regex.split(r'\s+', line[1:], 1)[0]
                contigs[name] = ''
            else:
                # Handle malformatted FASTA
                if ((name == None) or (name == "")):
                    raise ValueError('FASTA file malformatted')
                else:
                    contigs[name] += line
    
    return contigs