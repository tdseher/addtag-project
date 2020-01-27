#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/aligners/inbuilt.py

# List general Python imports
import sys
import os
import logging
logger = logging.getLogger(__name__)

# import non-standard package
import regex

# import AddTag-specific packages
sys.path.append( os.path.join( os.path.dirname(__file__), os.path.pardir ) )
import utils
import nucleotides
import cigarstrings
from .aligner import PairwiseAligner, Record

class Inbuilt(PairwiseAligner):
    logger = logger.getChild(__qualname__)

    def __init__(self):
        super().__init__(
            name="inbuilt",
            authors=['Seher, Thaddeus D.'],
            title='',
            journal='',
            issuing='',
            year=2016,
            doi='',
            input='fasta',
            output='sam',
            truncated=False,
            #classification='pairwise'
        )

    def index(self, fasta, output_filename, output_folder, threads, *args, **kwargs):
        pass

    def align(self, query, subject, output_filename, output_folder, threads, *args, **kwargs):
        pass

    def load_record(self, flo, *args, **kwargs):
        pass

def find_target_matches(compiled_regex, contigs, overlap=False):
    '''
    Finds all instances of the compiled_regex within the contigs(dict),
    and returns list of matches as 0-indexed
    
    overlap = Include exhaustive search for overlapping sites.
              May increase computation time.
    returns [contig, start, stop, orientation, oriented-sequence]
    '''
    
    matches = []
    
    # target AAACCC
    #                               vvvvvv                    19:25
    #    contig  GGGGGGGGGGGGAAAAAAAAAACCCGGGTTTTTTTTTTTTTT
    #                       vvvvvv                            11:17 -to-> 42-17:42-11 -to-> 25:31
    # rc(contig) AAAAAAAAAAAAAACCCGGGTTTTTTTTTTCCCCCCCCCCCC
    
    for contig in contigs:
        # '+' orientation
        ref = contigs[contig]
        ref_len = len(ref)
        for m in compiled_regex.finditer(ref, overlapped=overlap):
            #s = m.start()
            #e = m.end()
            #matches.append([contig, s, e, '+', ref[s:e]])
            
            cigar = cigarstrings.match2cigar(m, specific=True)
            matches.append([contig, m.start(), m.end(), '+', m.group(0), cigar])
        
        # '-' orientation
        rc_ref = nucleotides.rc(ref)
        for m in compiled_regex.finditer(rc_ref, overlapped=overlap):
            #s = m.start()
            #e = m.end()
            ##print("match: ", s, e, m.groups(), rc_ref[s:e], ref[ref_len-e:ref_len-s], file=sys.stderr)
            #assert m.groups()[0] == rc_ref[s:e]
            #assert m.groups()[0] == rc(ref[ref_len-e:ref_len-s]) # m.groups()[0] should be what was matched...
            #assert ref_len == rc_ref_len
            #matches.append([contig, rc_ref_len-s, rc_ref_len-e, '-', rc_ref[s:e]])
            
            # Make the CIGAR string, and reverse it
            cigar = cigarstrings.match2cigar(m, specific=True)
            cigar = cigarstrings.reverse_cigar(cigar)
            
            matches.append([contig, ref_len-m.end(), ref_len-m.start(), '-', m.group(0), cigar])
            
    logger.info("AddTag matching finished")
    
    return matches

def align(query_fasta, subject_fasta, folder=os.getcwd()):
    """
    Aligns query FASTA to subject FASTA using regular expressions
    """
    name = os.path.splitext(os.path.basename(query_fasta))[0]
    sam_file = os.path.join(folder, name + '.sam')
    
    targets = utils.old_load_fasta_file(query_fasta)
    contigs = utils.old_load_fasta_file(subject_fasta)
    #target = 'TCCGGTACAKTGAKTTGTAC'
    with open(sam_file, 'w') as flo:
        for target_name, target_seq in targets.items():
            regex = nucleotides.build_regex(target_seq, max_errors=5)
            matches = find_target_matches(regex, contigs, overlap=True)
            for m in matches:
                #print(m, file=flo)
                # This should print a SAM file
                
                flags = 0
                if (m[3] == '-'):
                    flags |= 16
                
                mapq = 1
                cigar = m[5]
                qqual = 'I'*len(target_seq) # 'I' is chr(ord('!')+40)
                print("\t".join(map(str, [target_name, flags, m[0], m[1]+1, mapq, cigar, '*', 0, 0, target_seq, qqual])), file=flo)
                
                # I don't think the query needs to be disambiguated
                #for seq in nucleotides.disambiguate_iupac(m[4]): # Disambiguate the subject... (Don't need this)
                #    print("\t".join(map(str, [target_name, flag, m[0], m[1]+1, mapq, cigar, '*', 0, 0, seq, len(seq), hsuzhang.hsuzhang_score(target, seq)])), file=flo)

def test():
    """Code to test the classes and functions in 'source/nucleotides.py'"""
    
    print("=== align ===")
    align(r"D:\VirtualBox share\addtag-project\ADE2g\excision-query.fasta", r"D:\VirtualBox share\temp-genome.fasta", folder=r"D:\VirtualBox share")
    

if (__name__ == '__main__'):
    test()
