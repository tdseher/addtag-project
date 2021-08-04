#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/seher.py

if (__name__ == "__main__"):
    from algorithm import PairedSequenceAlgorithm
else:
    from .algorithm import PairedSequenceAlgorithm

class Linear(PairedSequenceAlgorithm):
    def __init__(self):
        super().__init__(
            name="Linear",
            authors=['Seher, Thaddeus D.', 'Nguyen, Namkha', 'Ramos, Diana', 'Bapat, Priyanka', 'Nobile, Clarissa J.', 'Sindi, Suzanne S.', 'Hernday, Aaron D.'],
            title='AddTag, a two-step approach with supporting software package that facilitates CRISPR/Cas-mediated precision genome editing',
            journal='G3 Genes|Genomes|Genetis',
            issuing='',
            year=2021,
            doi='https://doi.org/10.1093/g3journal/jkab216',
            #citation="AddTag",
            off_target=True,
            prefilter=False,
            postfilter=False,
            minimum=75.0,
            maximum=100.0,
            default=100.0,
            weight_str=None,
            rgn_list=('Cas9', 'Cas12a', 'Cas3')
        )
    
    def is_available(self):
        """
        Determines if the prerequisites for the Algorithm have been met.
        :return: True or False
        """
        return True
    
    def calculate(self, intended, potential, *args, **kwargs):
        on_sequence, on_side, on_target, on_pam, on_upstream, on_downstream = intended
        off_sequence, off_side, off_target, off_pam, off_upstream, off_downstream = potential
        
        if (on_side == '>'):
            # Cases for 5'-SPACER>PAM-3'
            return self.score(on_target, off_target)
        else:
            # Cases for 5'-PAM<SPACER-3'
            return self.score(on_target[::-1], off_target[::-1])
    
    def score(self, seq1, seq2):
        """Scores lower if substitutions near 3' end of the sequence
        Should be gRNA only, with no PAM
          Scores >94 are good.
        Insertions and deletions are strongly penalized (as they are not aligned)
        Does not take ambiguities into account (yet)
        """
        shorter_seq_len, longer_seq_len = sorted([len(seq1), len(seq2)])
        x = list(range(1, longer_seq_len+1))
        x_sum = sum(x)
        
        # Normalize scores so all elements of 'y' sum to 1.0
        y = [i/x_sum for i in x] # list(map(lambda i: i/x_sum, x))
        score = 1
        for i in range(longer_seq_len):
            if (i < shorter_seq_len):
                if (seq1[i] != seq2[i]):
                    score -= y[i]
            else:
                score -= y[i]
        
        return score*100

def test():
    Cas9_a = ('', '>', 'AAAATTAACTATAGGTAAAG', 'TGG', '', '')
    Cas9_b = ('', '>', 'AACATCAACTCTAGCTAACG', 'CGG', '', '')
    Cas12a_a = ('', '<', 'AAAATTAACTATAGGTAAAGACCA', 'TTTA', '', '')
    Cas12a_b = ('', '<', 'AACATCAACTCTAGCTAACGACCA', 'TTTG', '', '')
    
    print("=== Linear ===")
    C = Linear()
    print('Cas9', C.calculate(Cas9_a, Cas9_b))
    print('Cas12a', C.calculate(Cas12a_a, Cas12a_b))

if (__name__ == "__main__"):
    test()
