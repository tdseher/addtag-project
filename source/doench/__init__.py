#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/doench/__init__.py

# Import standard packages
import os
import math
import fractions

# Import non-standard packages
import regex

# Import included AddTag-specific modules
from ..utils import rc

# This code only deals with mismatches
# It can be expanded to take insertions and deletions into account
# See Supplementary table 19

# Metric for "efficiency"
#  Doench
#   Range: 0-100. Linear regression model trained on 880 guides transfected
#   into human MOLM13/NB4/TF1 cells (three genes) and mouse cells
#   (six genes). Delivery: lentivirus. The Fusi score can be considered an
#   updated version this score, as their training data overlaps a lot.
#   See Doench et al.: http://www.nature.com/nbt/journal/v32/n12/full/nbt.3026.html

def load_params(file_path, sep="\t"):
    """Function to open and parse the tab-delimited 'params' file for
    Doench et al (2014), and return a list of tuples"""
    params = []
    with open(file_path, 'r') as flo:
        for line in flo:
            line = line.rstrip()
            if (len(line) > 0):
                m = regex.match(r'^\s*#', line)
                if not m:
                    sline = line.split(sep)
                    params.append((int(sline[0]), sline[1], float(sline[2])))
    return params

def on_target_score_2014(seq):
    """Function to calculate the sgRNA on-target efficacy score,
    as defined in Doench et al 2014"""
    # Code retrieved from:
    #  https://github.com/maximilianh/crisporWebsite/doenchScore.py
    # Calculates the sgRNA on-target efficacy score from the article
    # "Rational design of highly active sgRNAs for CRISPR-Cas9-mediated gene inactivation"
    # by J Doench et al. 2014
    # The authors' web tool is available at http://www.broadinstitute.org/rnai/public/analysis-tools/sgrna-design
    # Thanks to Cameron Mac Pherson at Pasteur Paris for fixing my original version. Maximilian Haeussler 2014
    
    params = load_params(os.path.join(os.path.dirname(__file__), 'doench_params.txt'))
    
    gcHigh = -0.1665878
    gcLow = -0.2026259
    
    score = 0.59763615 # also called  'intercept'
    
    guideSeq = seq[4:24]
    gcCount = guideSeq.count("G") + guideSeq.count("C")
    if gcCount <= 10:
        gcWeight = gcLow
    if gcCount > 10:
        gcWeight = gcHigh
    score += abs(10-gcCount)*gcWeight
    
    for pos, modelSeq, weight in params:
        subSeq = seq[pos:pos+len(modelSeq)]
        if subSeq==modelSeq:
            score += weight
    return 1.0/(1.0+math.exp(-score))

def test_on_target_score_2014():
    # Code for testing this function:
    test_sequences = {
        "TATAGCTGCGATCTGAGGTAGGGAGGGACC": 0.713089368437,
        "TCCGCACCTGTCACGGTCGGGGCTTGGCGC": 0.0189838463593,
    }
    for seq in test_sequences:
        print("Observed:", on_target_score_2014(seq))
        print("Expected:", test_sequences[seq])

def dump_mismatch_scores(file_path, sep="\t"):
    """Convert the pickled dictionary of mismatch scores, defined by
    Doench et al (2016), into a tab-delimited format
    """
    # This looks to be 5-dimensional data:
    #  1  2  3 4/5                 4  5
    # rG:dA, 5 0.3                 3/10
    # rC:dT,12 0.7142857140000001  5/ 7
    # rC:dT,18 0.538461538         7/13
    # rU:dG,14 0.28571428600000004 2/ 7
    # rC:dA,13 0.7                 7/10
    # rC:dC,18 0.133333333         2/15
    # ...etc...
    import pickle
    
    # Unpickle the input data file
    # by default, file_path should equal 'mismatch_score.pkl'
    mismatch_scores = pickle.load(open(file_path,'rb'))
    
    # Extract the data from each key:value pair
    lines = []
    for k in mismatch_scores:
        temp1 =  k.split(':')
        temp2 = temp1[1].split(',')
        r = temp1[0].replace('r', '') # string
        d = temp2[0].replace('d', '') # string
        pos = temp2[1] # already a string
        score = str(mismatch_scores[k]) # float
        approx = str(fractions.Fraction(mismatch_scores[k]).limit_denominator()) # fraction
        lines.append([r, d, pos, score, approx])
    
    output = sep.join(['# r', 'd', 'position', 'score', 'approximation']) + '\n'
    for sline in sorted(lines, key=lambda x: (x[0], x[1], int(x[2]))):
        output += sep.join(map(str, sline)) + '\n'
    
    return output

def dump_pam_scores(file_path, sep="\t"):
    """Convert the pickled dictionary of pam scores, defined by
    Doench et al (2016), into a tab-delimited format
    """
    import pickle
    # Unpickle the input data file
    # original file_path is 'pam_scores.pkl'
    pam_scores = pickle.load(open(file_path,'rb'))
    
    output = sep.join(['# ?', 'score', 'approximation']) + '\n'
    for k in pam_scores:
        output += sep.join(map(str, [k, pam_scores[k], fractions.Fraction(pam_scores[k]).limit_denominator()])) + '\n'
    
    return output

def load_mismatch_scores(file_path, sep="\t", approximation=False):
    """Load the mismatch scores defined by Doench et al (2016)"""
    # Data extracted into tab-delimited text as 'doench_mismatch_scores.txt'
    
    mismatch_scores = {}
    with open(file_path, 'r') as flo:
        for line in flo:
            line = line.rstrip()
            if (len(line) > 0):
                m = regex.match(r'^\s*#', line)
                if not m:
                    sline = line.split(sep)
                    k = 'r' + sline[0] + ':d' + sline[1] + ',' + sline[2]
                    if approximation:
                        mismatch_scores[k] = float(fractions.Fraction(sline[4]))
                    else:
                        mismatch_scores[k] = float(sline[3])
    return mismatch_scores

def load_pam_scores(file_path, sep="\t", approximation=False):
    """Load the PAM scors defined by Doench et al (2016)
    These represent the 'NGG' scores, excluding the 'N'
    """
    pam_scores = {}
    with open(file_path, 'r') as flo:
        for line in flo:
            line = line.rstrip()
            if (len(line) > 0):
                m = regex.match(r'^\s*#', line)
                if not m:
                    sline = line.split(sep)
                    if approximation:
                        pam_scores[sline[0]] = float(fractions.Fraction(sline[2]))
                    else:
                        pam_scores[sline[0]] = float(sline[1])
    return pam_scores

def on_target_score_2016_loader(wt, off):
    """Calculate the CFD score for a given sgRNA and off-target sequence.
    input:
      wt:  23 bp on-target nucleotide sequence. Contains 20 bp sgRNA + 3 bp PAM sequence.
      off: 23 bp off-target nucleotide sequence (20 bp + 3 bp PAM).
      Sequences must contain ONLY A, C, G, and T.
    output:
      CFD score
    """
    # Only calculate Doench score with cannonical nucleotides
    # Otherwise return a score of 0
    m_wt = regex.search('[^ATCG]', wt)
    m_off = regex.search('[^ATCG]', off)
    if (m_wt == None) and (m_off == None):
        pam = off[-2:] # Get the PAM site, excluding the 'N'
        sg = off[:-3]
        cfd_score = on_target_score_2016(wt, sg, pam)
        return cfd_score
    return 0.0

def on_target_score_2016(wt, sg, pam):
    """Also called CFD score.
    Cutting Frequency Determination (CFD) score calculator
    
    The on-target score is from Doench, Fusi et al. (2016)
    doi:10.1038/nbt.3437 and measures activity at the target location.
    Higher scores give higher confidence that the guide will be active.
    Scores range from 0-100, and should be used to rank guides relative to
    each other.
    
    The on-target score represents the cleavage efficiency of Cas96.
    You can think of the score as the probability a given gRNA will be in
    top 20% of cleavage activity. Note that the scoring system is not linear,
    and only 5% of gRNAs receive a score of 60 or higher.
    """
    
   # Calculates CFD score
    score = 1
    sg = sg.replace('T','U')
    wt = wt.replace('T','U')
    s_list = list(sg) # convert sg string to a list
    wt_list = list(wt)
    for i, sl in enumerate(s_list):
        #if wt_list[i] == sl:
        #    score*=1
        #else:
        if (wt_list[i] != sl):
            key = 'r'+wt_list[i]+':d'+rc(sl, kind="rna")+','+str(i+1)
            score*= MISMATCH_SCORES[key]
    score*=PAM_SCORES[pam]
    return (score*100)

def test():
    test_on_target_score_2014()
    
    a = 'CGATGGCTTGGATCGATTGAC'
    b = 'CGATGGCTTCGATCGATTGAC'
    c = 'CGATGGCTTCGAGCGATTGAC'
    print(on_target_score_2016_loader(a, a))
    print(on_target_score_2016_loader(a, b))
    print(on_target_score_2016_loader(a, c))
    

# Define module variables
# Load scores from the data files
try:
    MISMATCH_SCORES = load_mismatch_scores(os.path.join(os.path.dirname(__file__), 'doench_mismatch_scores.txt'))
    PAM_SCORES = load_pam_scores(os.path.join(os.path.dirname(__file__), 'doench_pam_scores.txt'))
except: 
    raise Exception("Could not find file with mismatch scores or PAM scores")

if (__name__ == "__main__"):
    test()