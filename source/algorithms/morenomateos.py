#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/morenomateos/__init__.py

# Import standard packages
import sys
import os

# Import non-standard packages
import regex

# Metric for "efficiency"
# Moreno-Mateos et al (2015)
#   Range: mostly 0-100. Linear regression model, trained on data from 1000
#   guides on >100 genes, from zebrafish 1-cell stage embryos injected with
#   mRNA.
# See Moreno-Mateos et al (2015):
#   http://www.nature.com/nmeth/journal/v12/n10/full/nmeth.3543.html
#   http://dx.doi.org/10.1038/nmeth.3543
#   http://crisprscan.org

paramsCRISPRscan = [
# converted excel table of logistic regression weights with 1-based positions
('AA',18,-0.097377097),
('TT',18,-0.094424075),('TT',13,-0.08618771),('CT',26,-0.084264893),('GC',25,-0.073453609),
('T',21,-0.068730497),('TG',23,-0.066388075),('AG',23,-0.054338456),('G',30,-0.046315914),
('A',4,-0.042153521),('AG',34,-0.041935908),('GA',34,-0.037797707),('A',18,-0.033820432),
('C',25,-0.031648353),('C',31,-0.030715556),('G',1,-0.029693709),('C',16,-0.021638609),
('A',14,-0.018487229),('A',11,-0.018287292),('T',34,-0.017647692),('AA',10,-0.016905415),
('A',19,-0.015576499),('G',34,-0.014167123),('C',30,-0.013182733),('GA',31,-0.01227989),
('T',24,-0.011996172),('A',15,-0.010595296),('G',4,-0.005448869),('GG',9,-0.00157799),
('T',23,-0.001422243),('C',15,-0.000477727),('C',26,-0.000368973),('T',27,-0.000280845),
('A',31,0.00158975),('GT',18,0.002391744),('C',9,0.002449224),('GA',20,0.009740799),
('A',25,0.010506405),('A',12,0.011633235),('A',32,0.012435231),('T',22,0.013224035),
('C',20,0.015089514),('G',17,0.01549378),('G',18,0.016457816),('T',30,0.017263162),
('A',13,0.017628924),('G',19,0.017916844),('A',27,0.019126815),('G',11,0.020929039),
('TG',3,0.022949996),('GC',3,0.024681785),('G',14,0.025116714),('GG',10,0.026802158),
('G',12,0.027591138),('G',32,0.03071249),('A',22,0.031930909),('G',20,0.033957008),
('C',21,0.034262921),('TT',17,0.03492881),('T',13,0.035445171),('G',26,0.036146649),
('A',24,0.037466478),('C',22,0.03763162),('G',16,0.037970942),('GG',12,0.041883009),
('TG',18,0.045908991),('TG',31,0.048136812),('A',35,0.048596259),('G',15,0.051129717),
('C',24,0.052972314),('TG',15,0.053372822),('GT',11,0.053678436),('GC',9,0.054171402),
('CA',30,0.057759851),('GT',24,0.060952114),('G',13,0.061360905),('CA',24,0.06221937),
('AG',10,0.063717093),('G',10,0.067739182),('C',13,0.069495944),('GT',31,0.07342535),
('GG',13,0.074355848),('C',27,0.079933922),('G',27,0.085151052),('CC',21,0.088919601),
('CC',23,0.095072286),('G',22,0.10114438),('G',24,0.105488325),('GT',23,0.106718563),
('GG',25,0.111559441),('G',9,0.114600681)]

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


def morenomateos_score(seq, pam, upstream='', downstream=''):
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

def calcCrisprScanScores(seqs):
    """ input is a 35bp long sequence: 6bp 5', 20bp guide, 3 bp PAM and 6bp 3'
    >>> calcCrisprScanScores(["TCCTCTGGTGGCGCTGCTGGATGGACGGGACTGTA"])
    [77]
    >>> calcCrisprScanScores(["TCCTCTNGTGGCGCTGCTGGATGGACGGGACTGTA"])
    [77]
    """
    scores = []
    for seq in seqs:
        assert(len(seq)==35)
        intercept = 0.183930943629
        score = intercept
        for modelSeq, pos, weight in paramsCRISPRscan:
            subSeq = seq[pos-1:pos+len(modelSeq)-1]
            if subSeq==modelSeq:
                score += weight
        #scores.append(int(100*score))
        scores.append(100*score)
    return scores

def test():
    """Code to test the classes and functions in 'source/morenomateos/__init__.py'"""
    
    seqs = [
        'TCCTCTGGTGGCGCTGCTGGATGGACGGGACTGTA', # [77.3863488629] 77.36304776290001
        'TCCTCTNGTGGCGCTGCTGGATGGACGGGACTGTA', # [77.3863488629] 77.36304776290001
        'AACGCTGTACGCTAGCTACCGATYYGCGACGCAAT', # [36.0535749629] 37.916735162900004
        'AAGTGTCGACTCCCGCTCTCAAAGAGCGGAGCTCC', # [50.76752636289999] 50.767526362900014
    ]
    for s in seqs:
        print(calcCrisprScanScores([s]), morenomateos_score(s[6:26], s[26:29]))

# Load scores from the data files
try:
    SCORES = load_scores(os.path.join(os.path.dirname(__file__), 'morenomateos_scores.txt'))
except FileNotFoundError:
    raise Exception("Could not find file with scores")


if (__name__ == "__main__"):
    test()
