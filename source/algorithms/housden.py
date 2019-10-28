#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/housden.py

# Import standard packages
import os
from math import log10

# Import non-standard packages
import regex

if (__name__ == "__main__"):
    from algorithm import SingleSequenceAlgorithm
else:
    from .algorithm import SingleSequenceAlgorithm

class Housden(SingleSequenceAlgorithm):
    def __init__(self):
        super().__init__(
            name="Housden",
            authors=['Housden, Benjamin E.', 'Valvezan, Alexander J.', 'Kelley, Colleen', 'Sopko, Richelle', 'Hu, Yanhui', 'Roesel, Charles', 'Lin, Shuailiang', 'Buckner, Michael', 'Tao, Rong', 'Yilmazel, Bahar', 'Mohr, Stephanie E.', 'Manning, Brendan D.', 'Perrimon, Norbert'],
            title='Identification of potential drug targets for tuberous sclerosis complex by synthetic screens combining CRISPR-based knockouts with RNAi',
            journal='Science Signaling',
            issuing='8(393):rs9',
            year=2015,
            doi='https://doi.org/10.1126/scisignal.aab3729',
            #citation="Housden, et al. Identification of potential drug targets for tuberous sclerosis complex by synthetic screens combining CRISPR-based knockouts with RNAi. Science Signaling 8(393), rs9 (2015)",
            off_target=True,
            on_target=True,
            prefilter=False,
            postfilter=False,
            minimum=0.0,
            maximum=13.0,
            default=None
        )
    
    def calculate(self, intended, *args, **kwargs):
        sequence, target, pam, upstream, downstream = intended
        
        return self.score(target)
    
    def score(self, seq):
        """
        Calculate Housden et al (2015).
        Removed rounding to one decimal point.
        Based on code from the CRISPOR paper.
        
        Input should be <=20nt spacer, with 1-distal, and 20-proximal to PAM.
        
        The efficiency score was calculated on the basis of a probability matrix
        computed using the in vitro cell line data. It reflects a cumulative P-value
        for high efficiency of each nucleotide from position 1 to 20, with higher
        values representing higher efficiency.
        
        Return scores range are 0 if invalid, or from 1.14-12.81 if valid.
        """
        seq = seq[-20:].upper()
        m = regex.search('[^ATCGatcg]', seq)
        if (m != None):
            return 0.0
        else:
            #assert(len(seq)==20)
            
            score = 1.0
            for i in range(-len(seq), 0):
                score *= SCORES[seq[i]][i]
            
            score = -1*log10(score)
            
            return score

def load_scores(file_path, sep="\t"):
    """Function to open and parse the tab-delimited 'scores' file for
    Housden et al (2015), and returns a dict"""
    scores = {
        'A': [],
        'C': [],
        'G': [],
        'T': [],
    }
    with open(file_path, 'r') as flo:
        for line in flo:
            line = line.rstrip()
            if (len(line) > 0):
                m = regex.match(r'^\s*#', line)
                if not m:
                    sline = line.split(sep)
                    scores['A'].append(float(sline[1]))
                    scores['C'].append(float(sline[2]))
                    scores['G'].append(float(sline[3]))
                    scores['T'].append(float(sline[4]))
    return scores

# Load scores from the data file
try:
    # Housden matrix: an array of 4x20=80 floats. The first 20 are for A,
    # the next 20 for C, then G, then T. Listed in Figure S1B
    SCORES = load_scores(os.path.join(os.path.dirname(__file__), 'housden_scores.txt'))
except FileNotFoundError:
    raise Exception("Could not find file with Housden scores")

# Theoretical minimum value
score_min = 1.0
for i in range(20):
    score_min *= max([SCORES['A'][i], SCORES['C'][i], SCORES['G'][i], SCORES['T'][i]])
score_min = -1*log10(score_min) # 1.140626469071239

# Theoretical maximum value
score_max = 1.0
for i in range(20):
    score_max *= min([SCORES['A'][i], SCORES['C'][i], SCORES['G'][i], SCORES['T'][i]])
score_max = -1*log10(score_max) # 12.812924105117725

def test():
    a = ('', 'AAAATTAACTATAGGTAAAG', '', '', '')
    b = ('', 'AACATCAACTCTAGCTAACG', '', '', '')
    c = ('', 'ATCTGACCTCCCGGCTAATT', '', '', '')
    d = ('', 'GGATGGAAATCTGACCTCCCGGCTAATT', '', '', '')
    e = ('', 'TGACCTCCCGGCTAATT', '', '', '')
    f = ('', 'ATCTGACCTCCCGGCTTATT', '', '', '')
    g = ('', 'ATCAGGCATCCCGGCTAGGT', '', '', '')
    
    print("=== Housden ===")
    C = Housden()
    print(C.calculate(a)) # 5.069672946242666
    print(C.calculate(b)) # 4.475098675704691
    print(C.calculate(c)) # 6.932429898673437
    print(C.calculate(d)) # -0.0/6.932429898673437
    print(C.calculate(e)) #  0.0/6.299595848542124
    print(C.calculate(f)) # 7.4666080475144145
    print(C.calculate(g)) # 6.784819433460606

if (__name__ == "__main__"):
    test()
