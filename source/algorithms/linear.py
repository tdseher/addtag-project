#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/linear.py

if (__name__ == "__main__"):
    from algorithm import PairedSequenceAlgorithm
else:
    from .algorithm import PairedSequenceAlgorithm

class Linear(PairedSequenceAlgorithm):
    def __init__(self):
        super().__init__("Linear", "Seher", 2017,
            citation="AddTag",
            off_target=True,
            prefilter=False,
            postfilter=False,
            minimum=80.0,
            maximum=100.0,
            default=100.0
        )
    
    def calculate(self, intended, potential, *args, **kwargs):
        on_sequence, on_target, on_pam, on_upstream, on_downstream = intended
        off_sequence, off_target, off_pam, off_upstream, off_downstream = potential
        
        return self.score(on_target, off_target)
    
    def score(self, seq1, seq2):
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

def test():
    a = ('', 'AAAATTAACTATAGGTAAAG', 'TGG', '', '')
    b = ('', 'AACATCAACTCTAGCTAACG', 'CGG', '', '')
    
    print("=== Linear ===")
    C = Linear()
    print(C.calculate(a, b))

if (__name__ == "__main__"):
    test()
