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
# Some fancy code for running this script as stand-alone
if (__name__ == "__main__"):
    import sys
    # Relative path for package to import
    PACKAGE_PARENT = '..'
    # Obtain path of currently-running file
    SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
    # Convert to absolute path, and add to the pythonpath
    sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))
    
    from utils import rc
else:
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

def load_old_scores(file_path, sep="\t"):
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

def on_target_score_2014(seq, pam, upstream='', downstream=''):
    """Function to calculate the sgRNA on-target efficacy score,
    as defined in Doench et al 2014
      seq = the gRNA sequence
      pam = the PAM sequence
      upstream = 0-10 bp immediately upstream of the gRNA
      downstream = 0-10 bp immediately downstream of the PAM
    Should return a score between 0.0 & 100.0, with higher numbers being
    better.
    """
    # Code retrieved from:
    #  https://github.com/maximilianh/crisporWebsite/doenchScore.py
    # Calculates the sgRNA on-target efficacy score from the article
    # "Rational design of highly active sgRNAs for CRISPR-Cas9-mediated gene inactivation"
    # by J Doench et al. 2014
    # The authors' web tool is available at http://www.broadinstitute.org/rnai/public/analysis-tools/sgrna-design
    # Thanks to Cameron Mac Pherson at Pasteur Paris for fixing my original version. Maximilian Haeussler 2014
    
    # We anchor the scoring algorithm at the PAM sequence
    # Inefficient code
    new = ['-'] * 30
    for i, nt in enumerate(pam):
        new[24+i] = nt

    for i, nt in enumerate(seq[::-1]):
        new[24-1-i] = nt

    for i, nt in enumerate(upstream[::-1]):
        ind = 24-len(seq)-1-i
        if (ind < 0):
            break
        new[ind] = nt

    for i, nt in enumerate(downstream):
        ind = 24+len(pam)+i
        if (ind >= len(new)):
            break
        new[ind] = nt
    new = ''.join(new)
    
    score = 0.59763615 # also called  'intercept'
    
    #guideSeq = seq[4:24]
    guideSeq = new[4:24]
    gcCount = guideSeq.count("G") + guideSeq.count("C")
    if gcCount <= 10:
        gcWeight = -0.2026259 # gcLow
    else:
        gcWeight = -0.1665878 # gcHigh
    score += abs(10 - gcCount) * gcWeight
    
    # ...shouldn't it be pos-1, so the scores are 0-indexed???
    for pos, modelSeq, weight in OLD_SCORES:
        #if (seq[pos:pos + len(modelSeq)] == modelSeq):
        if (new[pos:pos + len(modelSeq)] == modelSeq):
            score += weight
    
    return (1.0/(1.0+math.exp(-score))) * 100

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

def on_target_score_2016(seq1, seq2, pam):
    """Calculate the CFD score for a given gRNA and off-target sequence.
    
    Input oligonucleotides must be aligned with no gaps, and contain
    ONLY A, C, G, & T residues:
      seq1 = 20 bp genome sequence
      seq2 = 20 bp candidate gRNA sequence
       pam = 3 bp PAM sequence downstream of both seq1 & seq2
    
    Output:
      Cutting Frequency Determination (CFD) score
    
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
    # Return a score of 0 if sequences are of different length
    if (len(seq1) != len(seq2)):
        return 0.0
    
    # Only calculate Doench score with cannonical nucleotides
    m1 = regex.search('[^ATCG]', seq1)
    m2 = regex.search('[^ATCG]', seq2)
    mp = regex.search('[^ATCG]', pam[-2:])
    if ((m1 != None) or (m2 != None) or (mp != None)):
        # Otherwise return a score of 0
        return 0.0
    
    # Start with a perfect score
    score = 1
    
    # Ensure inputs are RNA sequences
    seq2 = seq2.replace('T','U')
    seq1 = seq1.replace('T','U')
    
    # Reduce the score multiplicatively if there is a mismatch
    for i in range(len(seq2)):
        if (seq1[i] != seq2[i]):
            key = 'r'+seq1[i]+':d'+rc(seq2[i], kind="rna")+','+str(i+1)
            score *= MISMATCH_SCORES[key]
    
    # Modulate the aggregate mismatch score by the PAM score
    # Exclude the most upstream residue in PAM motif (usually 'N'),
    # as it does not significantly affect score calculation
    score *= PAM_SCORES[pam[-2:]]
    
    return score * 100

def test():
    # Test Doench 2014 score:
    print("=== Doench 2014 ===")
    test_sequences = {
        ('TATA', 'GCTGCGATCTGAGGTAGGGA', 'GGG', 'ACC'): 71.3089368437,
        ('TATAGCT', 'GCGATCTGAGGTAGGGA', 'GGG', 'ACC'): 71.3089368437,
        ('',     'GCTGCGATCTGAGGTAGGGA', 'GGG', 'ACC'): None,
        ('TATA', 'GCTGCGATCTGAGGTAGGGA', 'GGG', ''   ): None,
        ('',     'GCTGCGATCTGAGGTAGGGA', 'GGG', ''   ): None,
        ('TCCG', 'CACCTGTCACGGTCGGGGCT', 'TGG', 'CGC'): 1.89838463593,
        ('--C-', 'GC-----TAGAAA-TCCG-A', 'G--', 'T-G'): None,
    }
    # Old test code
    #test_sequences = {
    #    "TATAGCTGCGATCTGAGGTAGGGAGGGACC": 0.713089368437,
    #    "TCCGCACCTGTCACGGTCGGGGCTTGGCGC": 0.0189838463593,
    #}
    for upstream, seq, pam, downstream in test_sequences:
        #print("Observed:", on_target_score_2014(seq))
        print((upstream, seq, pam, downstream))
        print("Observed:", on_target_score_2014(seq, pam, upstream, downstream))
        print("Expected:", test_sequences[(upstream, seq, pam, downstream)])
    
    # Test Doench 2016 score
    print("=== Doench 2016 ===")
    pam = 'AGG'
    gRNAa = 'CGATGGCTTGGATCGATTGA'
    gRNAb = 'CGTTGGCTTGGATCGATTGA'
    gRNAc = 'CGATGGCTTCGATCGATTGA'
    gRNAd = 'CGATGGCTTCGAGCGATTGA'
    print(on_target_score_2016(gRNAa, gRNAa, pam))
    print(on_target_score_2016(gRNAa, gRNAb, pam))
    print(on_target_score_2016(gRNAa, gRNAc, pam))
    print(on_target_score_2016(gRNAa, gRNAd, pam))
    

# Define module variables
# Load scores from the data files
try:
    MISMATCH_SCORES = load_mismatch_scores(os.path.join(os.path.dirname(__file__), 'doench_mismatch_scores.txt'))
    PAM_SCORES = load_pam_scores(os.path.join(os.path.dirname(__file__), 'doench_pam_scores.txt'))
    OLD_SCORES = load_old_scores(os.path.join(os.path.dirname(__file__), 'doench_params.txt'))
except: 
    raise Exception("Could not find file with mismatch scores or PAM scores")

if (__name__ == "__main__"):
    test()