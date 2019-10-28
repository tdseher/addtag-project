#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/housden.py

# Import standard packages
import os

# Import non-standard packages
import regex

if (__name__ == "__main__"):
    from algorithm import SingleSequenceAlgorithm
else:
    from .algorithm import SingleSequenceAlgorithm

# Metric for "efficiency"
# Moreno-Mateos et al (2015)
#   Range: mostly 0-100. Linear regression model, trained on data from 1000
#   guides on >100 genes, from zebrafish 1-cell stage embryos injected with
#   mRNA.
# See Moreno-Mateos et al (2015):
#   http://www.nature.com/nmeth/journal/v12/n10/full/nmeth.3543.html
#   http://dx.doi.org/10.1038/nmeth.3543
#   http://crisprscan.org

class MorenoMateos(SingleSequenceAlgorithm):
    def __init__(self):
        super().__init__(
            name="Moreno-Mateos",
            authors=['Moreno-Mateos, Miguel A.', 'Vejnar, Charles E.', 'Beaudoin, Jean-Denis', 'Fernandez, Juan P.', 'Mis, Emily K.', 'Khokha, Mustafa K.', 'Giraldez, Antonio J.'],
            title='CRISPRscan: designing highly efficient sgRNAs for CRISPR-Cas9 targeting in vivo',
            journal='Nature Methods',
            issuing='12:982',
            year=2015,
            doi='https://doi.org/10.1038/nmeth.3543',
            #citation="Moreno-Mateos, et al. CRISPRscan: designing highly efficient sgRNAs for CRISPR-Cas9 targeting in vivo. Nature Methods 12, 982–988 (2015)",
            off_target=True,
            on_target=True,
            prefilter=False,
            postfilter=False,
            minimum=1.0,
            maximum=100.0,
            default=None
        )
    
    def calculate(self, intended, *args, **kwargs):
        sequence, target, pam, upstream, downstream = intended
        
        return self.score(target, pam, upstream, downstream)
    
    def score(self, seq, pam, upstream='', downstream=''):
        """
        CrisprScan: Moreno-Mateos et al (2015)
        """
        #assert(len(seq)==35)
        
        # The scoring of the 35bp input sequence:
        # positions
        #  1-6    6 nt   5' upstream sequence
        #  7-26  20 nt   gRNA
        #  27-29  3 nt   PAM
        #  30-35  6 nt   3' downstream sequence
        
        # We anchor the scoring algorithm at the PAM sequence
        # Inefficient code
        new = ['-'] * 35
        for i, nt in enumerate(pam):
            new[26+i] = nt

        for i, nt in enumerate(seq[::-1]):
            new[26-1-i] = nt

        for i, nt in enumerate(upstream[::-1]):
            ind = 26-len(seq)-1-i
            if (ind < 0):
                break
            new[ind] = nt

        for i, nt in enumerate(downstream):
            ind = 26+len(pam)+i
            if (ind >= len(new)):
                break
            new[ind] = nt
        new = ''.join(new)
        
        intercept = 0.183930943629
        score = intercept
        for pos, model_seq, weight in SCORES:
            if (model_seq == new[pos:pos+len(model_seq)]):
                score += weight
        
        return 100*score

def load_scores(file_path, sep="\t"):
    """
    Function to open and parse the tab-delimited 'morenomateos_scores.txt'
    file for Moreno-Mateos et al (2015), and return a dict
    """
    #weights = {}
    #for i in ['A', 'C', 'G', 'T']:
    #    weights[i] = [0]*35
    #    for j in ['A', 'C', 'G', 'T']:
    #        weights[i+j] = [0]*35
    #
    #with open(file_path, 'r') as flo:
    #    for line in flo:
    #        line = line.rstrip()
    #        if (len(line) > 0):
    #            m = regex.match(r'^\s*#', line)
    #            if not m:
    #                sline = line.split(sep)
    #                weights[sline[1]][int(sline[0])-1] = float(sline[2])
    #
    #return weights
    weights = []
    with open(file_path, 'r') as flo:
        for line in flo:
            line = line.rstrip()
            if (len(line) > 0):
                m = regex.match(r'^\s*#', line)
                if not m:
                    sline = line.split(sep)
                    weights.append((int(sline[0])-1, sline[1], float(sline[2])))
    return weights

# Load scores from the data files
try:
    SCORES = load_scores(os.path.join(os.path.dirname(__file__), 'morenomateos_scores.txt'))
except FileNotFoundError:
    raise Exception("Could not find file with scores")

def test():
    """Code to test the classes and functions in 'source/morenomateos/__init__.py'"""
    
    seqs = [
        'TCCTCTGGTGGCGCTGCTGGATGGACGGGACTGTA', # [77.3863488629] 77.36304776290001
        'TCCTCTNGTGGCGCTGCTGGATGGACGGGACTGTA', # [77.3863488629] 77.36304776290001
        'AACGCTGTACGCTAGCTACCGATYYGCGACGCAAT', # [36.0535749629] 37.916735162900004
        'AAGTGTCGACTCCCGCTCTCAAAGAGCGGAGCTCC', # [50.76752636289999] 50.767526362900014
    ]
    print("=== MorenoMateos ===")
    C = MorenoMateos()
    for s in seqs:
        #print(calcCrisprScanScores([s]), morenomateos_score(s[6:26], s[26:29]))
        sequence, target, pam, upstream, downstream = s, s[6:26], s[26:29], s[:6], s[29:]
        print(C.calculate((sequence, target, pam, upstream, downstream)))

if (__name__ == "__main__"):
    test()
