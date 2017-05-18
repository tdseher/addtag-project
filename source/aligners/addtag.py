#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/aligners/addtag.py

# List general Python imports
import sys
import os

# import non-standard package
import regex

# import AddTag-specific packages
from .. import utils
from .. import nucleotides

def find_target_matches(compiled_regex, contigs, overlap=False):
    '''
    Finds all instances of the compiled_regex within the contigs(dict),
    and returns list of matches as 0-indexed
    
    overlap = Include exhaustive search for overlapping sites.
              May increase computation time.
    returns [contig, start, stop, orientation, oriented-sequence]
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

def align(query_fasta, subject_fasta, folder=os.getcwd()):
    """
    Aligns query FASTA to subject FASTA using regular expressions
    """
    name = os.path.splitext(os.path.basename(query_fasta))[0]
    sam_file = os.path.join(folder, name + '.sam')
    
    targets = utils.load_fasta_file(query_fasta)
    contigs = utils.load_fasta_file(subject_fasta)
    #target = 'TCCGGTACAKTGAKTTGTAC'
    with open(sam_file, 'w') as flo:
        for header in targets:
            regex = nucleotides.build_regex(targets[header], max_errors=5)
            matches = find_target_matches(regex, contigs, overlap=True)
            for m in matches:
                print(m, file=flo)
                # This should print a SAM file
                #for seq in nucleotides.disambiguate_iupac(m[4]):
                #    print(seq, len(seq), hsuzhang.hsuzhang_score(target, seq))

def test():
    """Code to test the classes and functions in 'source/nucleotides.py'"""
    
    print("=== align ===")
    
    

if (__name__ == '__main__'):
    test()