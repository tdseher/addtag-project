#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/proximalg.py

# The idea for this filter/score comes from the CHOPCHOPv2 implementation.
# This score will return a 0 if there is no G in the SPACER immediately adjacent to the PAM
# If the nt in SPACER immediately adjacent to PAM is a G, it will return 1.
# If the nt is ambiguous, then it will return a number between 0 and 1. E.G. N will return 0.25

# Import non-standard packages
if (__name__ == "__main__"):
    from algorithm import SingleSequenceAlgorithm
else:
    from .algorithm import SingleSequenceAlgorithm

class ProximalG(SingleSequenceAlgorithm):
    def __init__(self):
        super().__init__(
            name="ProximalG",
            authors=['Seher, Thaddeus D.'],
            title='',
            journal='',
            issuing='',
            year=2020,
            doi='',
            off_target=False,
            prefilter=True,
            postfilter=False,
            minimum=0.0,
            maximum=1.0,
            default=None,
            weight_str=None
        )
    
    def is_available(self):
        """
        Determines if the prerequisites for the Algorithm have been met.
        :return: True or False
        """
        return True
    
    def calculate(self, intended, *args, **kwargs):
        sequence, side, target, pam, upstream, downstream = intended
        
        return self.score(target, side)
    
    def score(self, seq, side):
        """
        Caclulates the if the single nt closest to PAM is a 'G' residue. IUPAC ambiguities okay.
        Treats IUPAC ambiguities as fractions.
        """
        
        iupac = {
            'a': 0.0,
            'c': 0.0,
            'g': 1.0,
            't': 0.0,
            'r': 0.5,
            'y': 0.0,
            'm': 0.0,
            'k': 0.5,
            'w': 0.0,
            's': 0.5,
            'b': 1.0/3,
            'd': 1.0/3,
            'h': 0.0,
            'v': 1.0/3,
            'n': 0.25,
            
            'A': 0.0,
            'C': 0.0,
            'G': 1.0,
            'T': 0.0,
            'R': 0.5,
            'Y': 0.0,
            'M': 0.0,
            'K': 0.5,
            'W': 0.0,
            'S': 0.5,
            'B': 1.0/3,
            'D': 1.0/3,
            'H': 0.0,
            'V': 1.0/3,
            'N': 0.25,
        }
        if (side == '>'):
            return iupac[seq[-1]]
        else: # (side == '<'):
            return iupac[seq[0]]

def test():
    """Code to test the functions and classes"""
    data = [
        ('', '>', 'GCTAGCGCTAACGCACTCAG', 'AGG', '', ''), # 1.00
        ('', '>', 'GCATATCGGCCGCCTATAAC', 'TGG', '', ''), # 0.00
        ('', '>', 'GGCTCGCAATAGAGAACCCN', 'GGG', '', ''), # 0.25
        ('', '<', 'GCTAGCGCTAACGCACTCAC', 'TTTC', '', ''), # 1.00
        ('', '<', 'CCATATCGGCCGCCTATAAG', 'TTTG', '', ''), # 0.00
        ('', '<', 'NGCTCGCAATAGAGAACCCA', 'TTTA', '', ''), # 0.25
    ]
    
    print("=== ProximalG ===")
    C = ProximalG()
    for d in data:
        print(C.calculate(d))

if (__name__ == "__main__"):
    test()
