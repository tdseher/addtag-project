#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/polyt.py

# Import non-standard packages
if (__name__ == "__main__"):
    from algorithm import SingleSequenceAlgorithm
else:
    from .algorithm import SingleSequenceAlgorithm

class PolyT(SingleSequenceAlgorithm):
    def __init__(self):
        super().__init__(
            name="PolyT",
            authors=['Seher, Thaddeus D.'],
            title='',
            journal='',
            issuing='',
            year=2017,
            doi='',
            #citation="AddTag",
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
        sequence, target, pam, upstream, downstream = intended
        
        return self.score(target)
    
    def score(self, seq):
        """
        Caclulates the longest Poly-T+ substring in seq. IUPAC ambiguities okay.
        Treats IUPAC ambiguities as fractions. So W=(1/2)T, and H=(1/3)T
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
        
        max_t_count = 0.0
        current_t_count = 0.0
        for i in seq:
            group = iupac[i]
            for j in ['T', 't']:
                if j in group:
                    current_t_count += 1.0/len(group)
                    break
            else:
                max_t_count = max(current_t_count, max_t_count)
                current_t_count = 0.0
        max_t_count = max(current_t_count, max_t_count)
        
        return max_t_count

def test():
    """Code to test the functions and classes"""
    
    a = ('', 'CGATGGCTAGGATCGATTGA', 'TGG', '', '') # 2.0
    b = ('', 'RYMKWSACGTbDHVNACGTA', 'TGG', '', '') # 1.9999999999999998
    c = ('',   'ATGSCTCGGATCGATTGA', 'AGG', '', '') # 2.0
    d = ('',             'ACATGTTG', 'AGG', '', '') # 2.0
    e = ('',           'ACATGTTNTG', 'AGG', '', '') # 3.25
    f = ('',                    'W', 'AGG', '', '') # 0.5
    g = ('',                    'h', 'AGG', '', '') # 0.3333333333333333
    h = ('',                    'N', 'AGG', '', '') # 0.25
    i = ('', 'AACTTTTTYCCATTNTTTTT', 'GGG', '', '') # 7.25
    
    print("=== PolyT ===")
    C = PolyT()
    print(C.calculate(a))
    print(C.calculate(b))
    print(C.calculate(c))
    print(C.calculate(d))
    print(C.calculate(e))
    print(C.calculate(f))
    print(C.calculate(g))
    print(C.calculate(h))
    print(C.calculate(i))

if (__name__ == "__main__"):
    test()