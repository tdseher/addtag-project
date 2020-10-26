#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/homopolymer.py

# Import non-standard packages
if (__name__ == "__main__"):
    from algorithm import SingleSequenceAlgorithm
else:
    from .algorithm import SingleSequenceAlgorithm

class Homopolymer(SingleSequenceAlgorithm):
    def __init__(self):
        super().__init__(
            name='Homopolymer',
            authors=[
                'Hough, Soren H.',
                'Kancleris, Kris',
                'Brody, Leigh',
                'Humphryes-Kirilov, Neil',
                'Wolanski, Joseph',
                'Dunaway, Keith',
                'Ajetunmobi, Ayokunmi',
                'Dillard, Victor',
            ],
            title='Guide Picker is a comprehensive design tool for visualizing and selecting guides for CRISPR experiments',
            journal='BMC Bioinformatics',
            issuing='18(1):167',
            year=2017,
            doi='https://doi.org/10.1186/s12859-017-1581-4',
            off_target=False,
            prefilter=True,
            postfilter=False,
            minimum=0.0,
            maximum=4.5,
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
        
        return self.score(target)
    
    def score(self, seq):
        """
        Caclulates the longest homopolymer substring in seq. IUPAC ambiguities okay.
        Treats IUPAC ambiguities as fractions. So W=(1/2)T, and H=(1/3)T
        """
        # Source implementation:
        #   "A score of 100 means the guide RNA does not contain a consecutive 4+ nucleotide homopolymer (desirable) 
        #   and a score of zero means it does contain a 4+ nucleotide homopolymer (undesirable).
        iupac = {
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
        
        #reset = (('A', 0.0), ('C', 0.0), ('G', 0.0), ('T', 0.0))
        #max_count = dict(reset) # {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        #current_count = dict(reset) # {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        max_count = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        current_count = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}
        for s in seq:
            group = iupac[s.upper()]
            for nt in ['A', 'C', 'G', 'T']: # should also not forget about lower case characters
                if nt in group:
                    current_count[nt] += 1.0/len(group)
                else:
                    max_count[nt] = max(max_count[nt], current_count[nt])
                    current_count[nt] = 0.0
        
        for nt in ['A', 'C' , 'G', 'T']:
            max_count[nt] = max(max_count[nt], current_count[nt])
        
        longest = max(max_count.values())
        
        return longest

def test():
    """Code to test the functions and classes"""
    data = [
        ('', '>', 'CGATGGCTAGGATCGATTGA', 'TGG', '', ''), # 2.0
        ('', '>', 'RYMKWSACGTbDHVNACGTA', 'TGG', '', ''), # 2.25
        ('', '>',   'ATGSCTCGGATCGATTGA', 'AGG', '', ''), # 2.0
        ('', '>',             'ACATGTTG', 'AGG', '', ''), # 2.0
        ('', '>',           'ACATGTTNTG', 'AGG', '', ''), # 3.25
        ('', '>',                    'W', 'AGG', '', ''), # 0.5
        ('', '>',                    'h', 'AGG', '', ''), # 0.3333333333333333
        ('', '>',                    'N', 'AGG', '', ''), # 0.25
        ('', '>', 'AACTTTTTYCCATTntTTTT', 'GGG', '', ''), # 7.25
        ('', '<', 'NNNAACC', 'TTTG', '', ''), # 2.75
    ]
    
    print("=== Homopolymer ===")
    C = Homopolymer()
    for d in data:
        print(C.calculate(d))

if (__name__ == "__main__"):
    test()