#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/scores.py

# List general Python imports
import sys

# import non-standard package
import regex

# import AddTag-specific packages
from . import nucleotides

class Score():
    """
    General class that should be subclassed when adding a new score
    """
    def __init__(self):
        """
        Create a new instance of the Score class
        """
        self.name = 'name'
        self.citation = 'AddTag 2017'
        self.include_in_off_target_score_calculation = True
        
        # Only include this score in 'off-target-score' calculation when this
        # score is between these values
        self.threshold = (1, 100)
        
    
    def calculate(self, sequence):
        """
        Algorithm to calculate score based on input sequence
        """
        return 0.0

def linear_score(seq1, seq2):
    """Scores lower if substitutions near 3' end of the sequence
    Should be gRNA only, with no PAM
      Scores >94 are good.
    Insertions and deletions are strongly penalized (as they are not aligned)
    Does not take ambiguities into account (yet)
    """
    shorter_seq_len, longer_seq_len = sorted([len(seq1), len(seq2)])
    x = list(range(longer_seq_len))
    x_sum = sum(x)
    y = list(map(lambda i: i/x_sum, x))
    score = 1
    for i in range(longer_seq_len):
        if (i < shorter_seq_len):
            if (seq1[i] != seq2[i]):
                score -= y[i]
        else:
            score -= y[i]
    
    return score*100

def gc_score(seq):
    """Caclulates the %GC content for seq. IUPAC ambiguities okay."""
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
    
    gc = 0.0
    for i in seq:
        k = iupac[i]
        for j in ['G', 'g', 'C', 'c']:
            if j in k:
                gc += 1.0/len(k)
    
    return 100*gc/len(seq)

# Off-target score
#  Hsu specificity score
#   The higher the specificity score, the lower are off-target effects in the genome.
#   The specificity score ranges from 0-100 and measures the uniqueness
#   of a guide in the genome. See Hsu et al. Nat Biotech 2013.
#   (http://dx.doi.org/10.1038/nbt.2647) We recommend values >50, where possible.
def off_target_score(off_target_scores, on_target_scores=(100,)):
    """Calculate the off-target score using the Zhong lab MIT Guide Score algorithm
    
    'off_target_scores' should be a list/tuple, with each element ranging from 0-100.
    'on_target_scores' should be a list/tuple, with each element ranging from 0-100.
    Returns off target score ranging from 0-100.
    
    Returns the naive proportion of gRNAs that will bind to a single target
    out of all possible targets.
    
    Higher values indicate greater specificity to intended target. Higher
    scores are better. Lower values mean that the candidate gRNA may bind
    nonspecifically. Guides with scores over 50 are good enough to be candidates
    (about 50% of the gRNA binding events will be to the intended target).
    The off-target score tells you the inverse probability of Cas9 off-target
    binding. A higher score means the sequence has less chance to bind to
    off-target sequences in the rest of the genome.
    
    The off-target score is from Hsu et al. (2013) doi:10.1038/nbt.2647
    and measures how specific the guide is to the target location. Scores should
    be used to rank guides relative to each other.
    
    This "MIT Guide Score" is defined on http://crispr.mit.edu/about
    Also called "Efficiency score", I think (or Eff Score) (see CRISPOR).
    """
    return  100.0*sum(on_target_scores)/(sum(on_target_scores)+sum(off_target_scores))

def r_score(seq1, seq2, length):
    """
    Input should not include PAM sequence.
    
    nnnnnnnnnnnnnnnnnnnnPAM query
                    iiii--- no mismatch within 4 nt of PAM
                iiiiiiii--- no mismatch within 8 nt of PAM
            iiiiiiiiiiii--- no mismatch within 12 nt of PAM
        iiiiiiiiiiiiiiii--- no mismatch within 16 nt of PAM
    """
    if (nucleotides.ridentities(seq1, seq2) >= length):
        return 1.0
    else:
        return 0.0

#def score_histogram(sequence, algorithms):
#    """Code that scores a gRNA sequence
#    Returns scores"""
#    # two types of off-target scores
#    #  CFD off-target score
#    #  MIT off-target score
#    
#    # Histogram of off-targets:
#    #  For each number of mismatches, the number of off-targets is indicated.
#    #  Example:
#    #   1-3-20-50-60    This means 1 off-target with 0 mismatches, 3 off-targets with 1 mismatch, 20 off-targets with 2 mismatches, etc.
#    #   0-2-5-10-20     These are the off-targets that have no mismatches in the 12 bp adjacent to the PAM. These are the most likely off-targets.
#    #   
#    #   Off-targets are considered if they are flanked by one of the motifs NGG, NAG or NGA.


def test():
    """Code to test the classes and functions in 'source/scores.py'"""
    
    print("=== linear_score ===")
    seq1 = 'GATGGATCGACCGAGATTTA'
    seq2 = 'GATCGATCGACCGAGATGTA'
    print(seq1, seq1, linear_score(seq1, seq1))
    print(seq1, seq2, linear_score(seq1, seq2))
    
    print("=== gc_score ===")
    
    print("=== off_target_score ===")
    

if (__name__ == '__main__'):
    test()
