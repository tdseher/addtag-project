#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/stemmer.py

# Import standard packages
import sys
import math

class Stemmer(PairedSequenceAlgorithm):
    def __init__(self):
        super().__init__("StemmerOffTarget", "Stemmer, et al", 2015,
            citation="Stemmer, et al. CCTop: An Intuitive, Flexible and Reliable CRISPR/Cas9 Target Prediction Tool. PLoS One 10(4), e0124633 (2015).",
            off_target=True,
            prefilter=False,
            postfilter=False,
            minimum=1.0,
            maximum=100.0,
            default=100.0
        )
    
    def calculate(self, intended, potential, *args, **kwargs):
        on_sequence, on_target, on_pam, on_upstream, on_downstream = intended
        off_sequence, off_target, off_pam, off_upstream, off_downstream = potential
        
        return self.score(on_target, off_target, on_pam)
    
    #def calculate(self, seq1, seq2, pam, *args, **kwargs):
    def score(self, seq1, seq2, pam, max_length=20):
        """
        The Stemmer algorithm computes a naive likelihood based on the number
        of mismatches and the positions of the mismatches within the off-target
        alignment
        """
        # Given that seq1 and seq2 are aligned
        mismatch_positions = [1, 4, 7] # positions of all mismatches between the two sequences, counted from the 5' end (starting with 1?)
        score = sum(1.2**x for x in mismatch_positions) # the more mismatches, and the closer the mismatch is to the PAM, the higher the score
        
        # if there is a mismatch at every position in a 20 nt sequence, the score is this
        # >>> sum(1.2**x for x in range(1, 21))
        # 224.02559954684835
        
        # If there is a mismatch at the 20th nt, the score is this:
        # >>> sum(1.2**x for x in [20])
        # 38.33759992447472
        
        
        # inverting...
        score = sum(1.2**(len(seq1)-x) for x in mismatch_positions)
        
        return 0.0
    
        # calculations 4/11/2019
        # Converting to a scale between 0 and 1
        # The Stemmer paper only considers mismatches,
        # However, we expand it so we can consider indels too.
        mismatch_positions = [1, 8, 15]
        denom = sum(1.2**n for n in range(1, len(seq1)+1))
        score = (denom - sum(1.2**n for n in mismatch_positions))/denom
        return score
    
    def off_target_score(self, off_targets):
        """
        We assume that the distance between off-target sites and exons are
        negligible.
        """
        off_targets = [score1, score2]
        
        # >>> import math
        # >>> math.log10(0)
        # Traceback (most recent call last):
        #   File "<pyshell#15>", line 1, in <module>
        #     math.log10(0)
        # ValueError: math domain error
        # >>> math.log10(1)
        # 0.0
        # >>> math.log10(10)
        # 1.0
        # >>> math.log10(100)
        # 2.0
        # >>> math.log10(100000)
        # 5.0
        
        dist = 1 # distance in nt between genomic position of potential off-target and the exon closest to it.
        number_off_targets = len(off_targets)
        
        if (number_off_targets > 0):
            score = sum(((math.log10(dist) + x)/number_off_targets) - number_off_targets for x in off_targets)
        else:
            score = 100.0
        
        #score = 0
        #for ot_score in off_targets:
        #    score += ((math.log10(dist) + ot_score)/number_off_targets) - number_off_targets
        return score
        
        # sum(1.2**(20-x) for x in mismatch_positions)
        # sum(((math.log10(dist) + x)/len(off_targets)) - len(off_targets) for x in off_targets)

def test():
    """Code to test the functions and classes"""
    
    a = ('',   'CGATGGCTAGGATCGATTGA', 'TGG', '', '')
    b = ('',   'RYMKWSACGTbDHVNACGTA', 'TGG', '', '')
    c = ('',     'ATGSCTCGGATCGATTGA', 'AGG', '', '')
    d = ('',   'GCGATGCGCAGCTAGGCCGG', 'CGG', '', '')
    e = ('',   'CGAAGGCTCGGACCGATTGA', 'GGG', '', '')
    f = ('',   'CGCTGGCTAGGATCGATTGA', 'AGG', '', '')
    g = ('',   'AAAATTAACTATAGGTAAAG', 'TGG', '', '')
    h = ('',   'AACATCAACTCTAGCTAACG', 'CGG', '', '')
    i = ('',   'AACATCAACTCTACCTAACG', 'CGG', 'CCGA', 'AACA')
    j = ('',   'GTTAGCGGTATGTATATGTG', 'TGG', 'GGGA', 'CTCA')
    k = ('', 'CTCAACATGGTATGTATATGTG', 'TGG', 'TCGA', 'TTCA')
    
    print("=== Doench2014 ===")
    C = Doench2014()
    print(C.calculate(a)) # 49.18420140306892
    print(C.calculate(b)) # 16.668179672077244
    print(C.calculate(c)) # 32.30335056113513
    print(C.calculate(d)) # 43.173287661422385
    print(C.calculate(e)) # 62.564461061593754
    print(C.calculate(f)) # 22.59673843434246
    
    print("=== Doench2016 ===")
    C = Doench2016()
    print(C.calculate(a, a)) # 100.0
    print(C.calculate(a, b)) # 0.0
    print(C.calculate(a, c)) # 0.0
    print(C.calculate(a, d)) # 0.008843903254636097
    print(C.calculate(a, e)) # 21.482277090941643
    print(C.calculate(a, f)) # 42.8571429
    print(C.calculate(j, k)) # 12.26204765820043

if (__name__ == "__main__"):
    test()
