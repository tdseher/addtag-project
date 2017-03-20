#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/housden/__init__.py

# Import standard packages
import os
from math import log10

# Import non-standard packages
import regex

# Metric for "efficiency"
#  Housden
#   Range: ~ 1-10. Weight matrix model trained on data from Drosophila mRNA
#   injections. See Housden et al.:
#   http://stke.sciencemag.org/content/8/393/rs9.long

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

def housden_score(seq):
    """
    Calculate Housden et al (2015) score and round to one decimal point.
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
        
        return round(score, 1) # round to one decimal point

def calcHousden(seqs):
    """
    Calc housden score and round to one decimal point.
    Based on java file Crispr.java received from the authors.
    >>> calcHousden(["ATCTGACCTCCCGGCTAATT"])
    [6.9]
    """
    # Housden matrix (see function below):
    # an array of 4x20=80 floats. The first 20 are for A, the next 20 for T, then C, then G
    # imported from original file received from the authors: matrix_final.txt
    factors = [
        0.4979, 0.7959, 0.7553, 0.6569, 0.9481, 0.7147, 0.437, 0.6212, 0.9077, 1.0, 0.1957, 0.7959, 0.6212, 0.8912, 1.0, 0.5485, 0.9942, 0.5485, 0.455, 1.0, \
        0.6699, 0.5485, 0.275, 0.5972, 0.6212, 0.7555, 1.0, 0.5131, 0.8608, 0.7553, 0.6569, 0.3417, 1.0, 0.016, 0.9146, 0.7555, 0.2906, 0.4979, 0.5485, 0.5131, 
        0.4979, 0.6869, 0.8528, 0.7643, 0.5325, 0.3417, 0.3417, 0.7643, 0.6434, 0.0092, 0.9331, 0.5325, 0.7272, 0.9708, 0.2905, 0.7272, 0.2957, 0.7918, 0.6434, 0.5062, \
        0.7918, 0.4461, 0.4851, 0.4461, 0.3417, 0.6869, 0.2417, 0.5485, 0.0947, 0.9256, 0.5325, 0.8308, 0.1255, 0.7918, 0.2544, 0.4461, 0.4979, 0.6212, 0.7918, 0.4461
    ]
    
    scores = []
    for seq in seqs:
        seq = seq.upper()
        if "N" in seq: # cannot do Ns
            scores.append(-1.0)
            continue
        if len(seq) != 20:
            scores.append(-1.0)
            continue
        #assert(len(seq)==20)
        
        # Columns associated with each nucleotide
        nuclToIndex = {"A":0,"T":1,"C":2,"G":3}

        score = 1.0
        for i in range(0, 20):
            nuclIndex = nuclToIndex[seq[i]]
            idx = (nuclIndex*20)+i
            score *= factors[idx]
        score = -1*log10(score)
        score = float("%0.1f" % score) # round to one decimal point
        scores.append(score)
    return scores

def test():
    """Code to test the classes and functions in 'source/housden/__init__.py'"""
    s = 'ATCTGACCTCCCGGCTAATT' # should be 6.9
    print(calcHousden([s]), housden_score(s))
    
    s = 'GGATGGAAATCTGACCTCCCGGCTAATT' # should be -1.0/6.9
    print(calcHousden([s]), housden_score(s))
    
    s = 'TGACCTCCCGGCTAATT' # should be -1.0/6.3
    print(calcHousden([s]), housden_score(s))
    
    s = 'ATCTGACCTCCCGGCTTATT' # should be 7.5
    print(calcHousden([s]), housden_score(s))
    
    s = 'ATCAGGCATCCCGGCTAGGT' # should be 6.8
    print(calcHousden([s]), housden_score(s))
    

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

if (__name__ == '__main__'):
    test()
