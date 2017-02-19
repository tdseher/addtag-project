#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/hsu/__init__.py

# List imports
import os
import regex

# Define module functions

# Off-target score
#  Hsu specificity score
#   The higher the specificity score, the lower are off-target effects in the genome.
#   The specificity score ranges from 0-100 and measures the uniqueness
#   of a guide in the genome. See Hsu et al. Nat Biotech 2013.
#   (http://dx.doi.org/10.1038/nbt.2647) We recommend values >50, where possible.

def off_target_score(guide, target):
    """The off-target score is from Hsu et al. (2013) doi:10.1038/nbt.2647
    and measures how specific the guide is to the target location. Guides with
    scores over 50 are good enough to be candidates. Higher scores are better.
    Scores range from 0-100, and should be used to rank guides relative
    to each other.
    
    The off-target score tells you the inverse probability of Cas9 off-target
    binding. A higher score means the sequence has less chance to bind to
    sequences in the rest of the genome.
    
    Returns the off-target score and the list of potential off-target sequences
    and their genomic coordinates.
    """
    pass

# This code pulled from crispor.py
# Code is python implementation of http://crispr.mit.edu/about

def load_scores(file_path, sep="\t"):
    """Function to open and parse the tab-delimited 'scores' file for
    Hsu et al (2013), and return a list"""
    scores = []
    with open(file_path, 'r') as flo:
        for line in flo:
            line = line.rstrip()
            if (len(line) > 0):
                m = regex.match(r'^\s*#', line)
                if not m:
                    sline = line.split(sep)
                    scores.append(float(sline[1]))
    return scores

def hsu_score(seq1, seq2, iupac=False):
    '''Calculate Hsu score between two oligonucleotide sequences.
    Ideally, seq1 and seq2 are each 20 bp long.
    Will only score the final 20nt. No penalization for less than 20nt.
    See 'Scores of single hits' (http://crispr.mit.edu/about)
    for algorithm used to score single offtargets.
        
    Within the first term, 'e' runs over the mismatch positions between guide
    and offtarget, with 'M' representing the experimentally-determined effect
    of mismatch position on targeting. (Hsu et al, Nature Biotechnology 2013)
    (http://www.nature.com/nbt/journal/v31/n9/full/nbt.2647.html)
    And terms two and three factoring in the effect of mean pairwise distance
    between mismatches 'd' and a dampening penalty for highly mismatched
    targets.
    '''
    
    # Truncate any nt on the left side. Keep the right side the same:
    #  input:  'ACGCTTTAGCGCCAGACTCAGT'
    #  output:   'GCTTTAGCGCCAGACTCAGT'
    seq1 = seq1[-20:]
    seq2 = seq2[-20:]
    
    dists = [] # distances between mismatches, for part 2
    mmCount = 0 # number of mismatches, for part 3
    lastMmPos = None # position of last mismatch, used to calculate distance

    score1 = 1.0
    if iupac:
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
        for pos in range(-min(len(seq1), len(seq2)), 0):
            if (len([x for x in iupac[seq1[pos]] if x in iupac[seq2[pos]]]) == 0):
                mmCount+=1
                if (lastMmPos != None):
                    dists.append(pos-lastMmPos)
                score1 *= 1-SCORES[pos]
                lastMmPos = pos
    else:
        for pos in range(-min(len(seq1), len(seq2)), 0):
            if (seq1[pos] != seq2[pos]):
                mmCount+=1
                if (lastMmPos != None):
                    dists.append(pos-lastMmPos)
                score1 *= 1-SCORES[pos]
                lastMmPos = pos
    
    # 2nd part of the score
    if (mmCount < 2): # special case, not shown in the paper
        score2 = 1.0
    else:
        avgDist = sum(dists)/len(dists)
        score2 = 1.0 / (((19-avgDist)/19.0) * 4 + 1)
    
    # 3rd part of the score
    if (mmCount == 0): # special case, not shown in the paper
        score3 = 1.0
    else:
        score3 = 1.0 / (mmCount**2)

    score = score1 * score2 * score3 * 100
    return score

def calcHitScore(string1, string2):
    '''Scores of Single Hits
    The actual algorithm used to score single offtargets is:
        \prod_{e\in{\mathcal{M}}}\left(1- W[e]\right)\times\frac{1}{\left(\frac{(19 - \bar{d})}{19}\times4 + 1\right)}\times\frac{1}{n^2_{mm}}
    Within the first term, e runs over the mismatch positions between guide and offtarget, with M:
        M = [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,\\ 0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]
    representing the experimentally-determined effect of mismatch position on targeting. (Hsu et al, Nature Biotechnology 2013) (http://www.nature.com/nbt/journal/v31/n9/full/nbt.2647.html)
    And terms two and three factoring in the effect of mean pairwise distance between mismatches (d) and a dampening penalty for highly mismatched targets.
    '''
    
    # Load scores from the data file
    #try:
    #    hitScoreM = load_scores(os.path.join(os.path.dirname(__file__), 'hsu_scores.txt'))
    #except: 
    #    raise Exception("Could not find file with Hsu scores")
    hitScoreM = SCORES
    
    # see 'Scores of single hits' on http://crispr.mit.edu/about
    # The Patrick Hsu weighting scheme
    # S. aureus requires 21bp long guides. We fudge by using only last 20bp
    if (len(string1) == len(string2) == 21):
        string1 = string1[-20:]
        string2 = string2[-20:]
    
    # If string lengths do not match, then return 0
    if not (len(string1)==len(string2)==20):
        return 0.0
    
    # Throw AssertionError if the lengths are not equal
    assert(len(string1)==len(string2)==20)

    dists = [] # distances between mismatches, for part 2
    mmCount = 0 # number of mismatches, for part 3
    lastMmPos = None # position of last mismatch, used to calculate distance

    score1 = 1.0
    for pos in range(0, len(string1)):
        if string1[pos]!=string2[pos]:
            mmCount+=1
            if lastMmPos!=None:
                dists.append(pos-lastMmPos)
            score1 *= 1-hitScoreM[pos]
            lastMmPos = pos
    
    # 2nd part of the score
    if mmCount<2: # special case, not shown in the paper
        score2 = 1.0
    else:
        avgDist = sum(dists)/len(dists)
        score2 = 1.0 / (((19-avgDist)/19.0) * 4 + 1)
    
    # 3rd part of the score
    if mmCount==0: # special case, not shown in the paper
        score3 = 1.0
    else:
        score3 = 1.0 / (mmCount**2)

    score = score1 * score2 * score3 * 100
    return score

def calcMitGuideScore(hitSum):
    """Aggregate Scores by Guide
    Once individual hits have been scored, each guide is assigned a score:
        S_{guide} = \frac{100}{100 + \sum_{i=1}^{n_{mm}}S_{hit}(h_i)}
    and colored according to a broad categorization of guide quality which,
    taken into account with the presence or absence of marked genes in
    high-scoring offtargets indicate the relative (un)favorability of using a
    particular guide for specific targeting in the query region.
    
    This MitGuideScore is defined on http://crispr.mit.edu/about 
    Input is the sum of all off-target hit scores.
    Returns the specificity of the guide.
    """
    score = 100 / (100+hitSum)
    score = int(round(score*100))
    return score

def test():
    # Test Hsu score
    print("=== Hsu 2013 ===")
    a = 'CGATGGCTWGGATCGATTGAC'
    b = 'AAGTGCTCTTAAGAGAAATTC'
    c = 'ATGSCTCGGATCGATTGAC'
    print(calcHitScore(a, a), hsu_score(a, a), hsu_score(a, a, True))
    print(calcHitScore(a, b), hsu_score(a, b), hsu_score(a, b, True))
    print(calcHitScore(a, c), hsu_score(a, c), hsu_score(a, c, True))

# Define module variables
# Load scores from the data file
try:
    SCORES = load_scores(os.path.join(os.path.dirname(__file__), 'hsu_scores.txt'))
except: 
    raise Exception("Could not find file with Hsu scores")

if (__name__ == "__main__"):
    test()