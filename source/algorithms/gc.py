#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/gc.py

if (__name__ == "__main__"):
    from algorithm import SingleSequenceAlgorithm
else:
    from .algorithm import SingleSequenceAlgorithm

class GC(SingleSequenceAlgorithm):
    def __init__(self):
        super().__init__("GC", "Seher", 2017,
            citation="AddTag",
            off_target=True,
            prefilter=True,
            postfilter=False,
            minimum=10.0,
            maximum=90.0,
            default=None
        )
    
    def calculate(self, intended, *args, **kwargs):
        sequence, target, pam, upstream, downstream = intended
        
        return self.score(target)
    
    def score(self, seq):
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

def test():
    a = ('', 'AAAATTAACTATAGGTAAAG', 'TGG', '', '')
    b = ('', 'AACATCAACTCTAGCTAACG', 'CGG', '', '')
    c = ('', 'AACATCACCTCTGGCTAACG', 'CGG', '', '')
    d = ('', 'ACCACCAACTCTAGCTGACG', 'CGG', '', '')
    e = ('',  'CCACCAACTCTAGCTGACG', 'CGG', '', '')
    
    print("=== GC ===")
    C = GC()
    print(C.calculate(a))
    print(C.calculate(b))
    print(C.calculate(c))
    print(C.calculate(d))
    print(C.calculate(e))

if (__name__ == "__main__"):
    test()