#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/first.py

# Import standard packages
import os

# Import non-standard packages
import regex

if (__name__ == "__main__"):
    from algorithm import PairedSequenceAlgorithm
else:
    from .algorithm import PairedSequenceAlgorithm

class HsuZhang(PairedSequenceAlgorithm):
    def __init__(self):
        super().__init__(
            name="Hsu-Zhang",
            authors=['Hsu, Patrick D.', 'Scott, David A.', 'Weinstein, Joshua A.', 'Ran, F. Ann', 'Konermann, Silvana', 'Agarwala, Vineeta', 'Li, Yinqing', 'Fine, Eli J.', 'Wu, Xuebing', 'Shalem, Ophir', 'Cradick, Thomas J.', 'Marraffini, Luciano A.', 'Bao, Gang', 'Zhang, Feng'],
            title='DNA targeting specificity of RNA-guided Cas9 nucleases',
            journal='Nature Biotechnology',
            issuing='31(9):827-832',
            year=2013,
            doi='https://doi.org/10.1038/nbt.2647',
            #citation="Hsu, et al. DNA targeting specificity of RNA-guided Cas9 nucleases. Nature Biotechnology 31, 827–832 (2013)",
            off_target=True,
            prefilter=False,
            postfilter=False,
            minimum=0.1,
            maximum=100.0,
            default=100.0,
            weight_str='Hsu-Zhang:90+1.8' # Severly penalize any score less than 95
        )
    
    # def weight(self, x):
    #     """Severely penalize any score less than 95"""
    #     return 1.0/(1+1.8**(90-x))
    
    def calculate(self, intended, potential, *args, **kwargs):
        on_sequence, on_target, on_pam, on_upstream, on_downstream = intended
        off_sequence, off_target, off_pam, off_upstream, off_downstream = potential
        
        if 'iupac' not in kwargs:
            kwargs['iupac'] = False
        
        return self.score(on_target, off_target, iupac=kwargs['iupac'])
    
    def score(self, seq1, seq2, iupac=False):
        '''Calculate Hsu score between two oligonucleotide sequences.
        (Also called 'MIT score') Should not include PAM.
        
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
        
        shorter, longer = sorted([seq1, seq2], key=len)
        
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
            for pos in range(-len(shorter), 0):
                if (len([x for x in iupac[seq1[pos]] if x in iupac[seq2[pos]]]) == 0):
                    # Consider 1/2, 1/3, and 1/4 mismatches for ambiguities...
                    mmCount+=1
                    if (lastMmPos != None):
                        dists.append(pos-lastMmPos)
                    score1 *= 1-SCORES[pos]
                    lastMmPos = pos
        else:
            for pos in range(-len(shorter), 0):
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

# Load scores from the data files when this module is imported
try:
    SCORES = load_scores(os.path.join(os.path.dirname(__file__), 'hsuzhang_scores.txt'))
except FileNotFoundError:
    raise Exception("Could not find file with Hsu scores")

def test():
    """Code to test the classes and functions in 'source/hsuzhang/__init__.py'"""
    
    print("=== HsuZhang ===")
    C = HsuZhang()
    a = ('', 'CGATGGCTWGGATCGATTGAC', '', '', '')
    b = ('', 'AAGTGCTCTTAAGAGAAATTC', '', '', '')
    c = ('',   'ATGSCTCGGATCGATTGAC', '', '', '')
    d = ('', 'AAAATTAACTATAGGTAAAG', '', '', '')
    e = ('', 'AATATAAAAAATAGGTAAAG', '', '', '')
    f = ('', 'AAAAGTAGCCATAGGAAAAG', '', '', '')
    print(C.calculate(a, a, iupac=True)) # 100.0
    print(C.calculate(a, b, iupac=True)) # 
    print(C.calculate(a, c, iupac=True)) # 100.0
    print(C.calculate(d, e)) # 0.46532338607757784
    print(C.calculate(d, f)) # 0.2341671161825727

if (__name__ == "__main__"):
    test()
