#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/utils.py

# Import standard packages
import sys
import gzip
import regex

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
        complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB') # exclude ws, WS
    elif (kind == "rna"):
        complements = str.maketrans('acgturymkbdhvACGTURYMKBDHV', 'ugcaayrkmvhdbTGCAAYRKMVHDB') # exclude ws, WS
    else:
        raise ValueError("'" + str(kind) + "' is an invalid argument for rc()")
    return seq.translate(complements)[::-1]

def read_fasta_file(filename):
    """Load contig sequences from file into dict()
    Primary sequence headers must be unique
    Can open a *.gz compressed file
    """
    
    if filename.endswith('.gz'):
        flo = gzip.open(filename, 'rt')
    else:
        flo = open(filename, 'r')
    
    contigs = {}
    #with open(filename, 'r') as flo:
    with flo:
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
    
    print('FASTA file parsed: {!r}'.format(filename), file=sys.stderr)
    return contigs

def disambiguate_iupac(iupac_sequence):
    """converts a string containing IUPAC nucleotide sequence to a list
    of non-iupac sequences.
    """
    iupac = {
        'a': ['a'],
        'c': ['c'],
        'g': ['g'],
        't': ['t'],
        'r': ['a', 'g'],
        'y': ['c', 't'],
        'm': ['a', 'c'],
        'k': ['g', 't'],
        'w': ['a', 't'],
        's': ['c', 'g'],
        'b': ['c', 'g', 't'],
        'd': ['a', 'g', 't'],
        'h': ['a', 'c', 't'],
        'v': ['a', 'c', 'g'],
        'n': ['a', 'c', 'g', 't'],
        
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
    sequences = ['']
    for char in iupac_sequence:
        # If the current character is an ambiguous base,
        # then copy all existing sequences
        for i in range(len(sequences)):
            for j in range(len(iupac[char]) - 1):
                sequences.append(sequences[i])
        # Sort the sequences... not very efficient
        sequences = sorted(sequences)
        # Add to each sequence a new, disambiguated nt
        for i in range(len(sequences)):
            sequences[i] += iupac[char][i % len(iupac[char])]
    
    # Throw an error if there is a non-cannonical nucleotide
    assert all(map(lambda seq: all(map(lambda x: x in 'AaCcGgTtUu', seq)), sequences)) == True
    
    return sequences

def build_regex(iupac_sequence, case_sensitive=False, max_substitutions=0, max_insertions=0, max_deletions=0, max_errors=0):
    """Build a regular expression for the nucleotide search, taking IUPAC
    ambiguities into account.
    """
    # Format the fuzzy matching restrictions
    fuzzy = []
    if (max_substitutions > 0):
        fuzzy.append('s<={}'.format(max_substitutions))
    if (max_insertions > 0):
        fuzzy.append('i<={}'.format(max_insertions))
    if (max_deletions > 0):
        fuzzy.append('d<={}'.format(max_deletions))
    if (max_errors > 0):
        fuzzy.append('e<={}'.format(max_errors))
    if (len(fuzzy) > 0):
        fuzzy = '{' + ','.join(fuzzy) + '}'
    else:
        fuzzy = ''
    
    # Choose the regex flags
    myflags = regex.ENHANCEMATCH | regex.IGNORECASE
    if case_sensitive:
        myflags = regex.ENHANCEMATCH
    #myflags = regex.IGNORECASE
    #if case_sensitive:
    #    myflags = 0
    
    # Convert the IUPAC sequence to an equivalent regex
    iupac = {
        'a': 'a',
        'c': 'c',
        'g': 'g',
        't': 't',
        'r': '[ag]',
        'y': '[ct]',
        'm': '[ac]',
        'k': '[gt]',
        'w': '[at]',
        's': '[cg]',
        'b': '[cgt]',
        'd': '[agt]',
        'h': '[act]',
        'v': '[acg]',
        'n': '[acgt]',
        
        'A': 'A',
        'C': 'C',
        'G': 'G',
        'T': 'T',
        'R': '[AG]',
        'Y': '[CT]',
        'M': '[AC]',
        'K': '[GT]',
        'W': '[AT]',
        'S': '[CG]',
        'B': '[CGT]',
        'D': '[AGT]',
        'H': '[ACT]',
        'V': '[ACG]',
        'N': '[ACGT]',
    }
    sequence = ''.join(map(lambda x: iupac[x], iupac_sequence))
    pattern = '(' + sequence + ')' + fuzzy
    print('Compiled regex: {!r}'.format(pattern), file=sys.stderr)
    compiled_regex = regex.compile(pattern, flags=myflags)
    return compiled_regex

def find_target_matches(compiled_regex, contigs, overlap=False):
    '''Finds all instances of the compiled_regex within the contigs(dict),
    and returns list of matches as 0-indexed
    [contig, start, stop, orientation, oriented-sequence]
    '''
    
    matches = []
    
    for contig in contigs:
        # '+' orientation
        ref = contigs[contig]
        ref_len = len(ref)
        for m in compiled_regex.finditer(ref, overlapped=overlap):
            s = m.start()
            e = m.end()
            matches.append([contig, s, e, '+', ref[s:e]])
        # '-' orientation
        rc_ref = rc(contigs[contig])
        rc_ref_len = len(rc_ref)
        for m in compiled_regex.finditer(rc_ref, overlapped=overlap):
            s = m.start()
            e = m.end()
            #print("match: ", s, e, m.groups(), rc_ref[s:e], ref[ref_len-e:ref_len-s], file=sys.stderr)
            assert m.groups()[0] == rc_ref[s:e]
            assert m.groups()[0] == rc(ref[ref_len-e:ref_len-s]) # m.groups()[0] should be what was matched...
            assert ref_len == rc_ref_len
            matches.append([contig, rc_ref_len-s, rc_ref_len-e, '-', rc_ref[s:e]])
    print("Matching finished", file=sys.stderr)
    
    return matches

