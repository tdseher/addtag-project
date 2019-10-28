#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/labuhn.py

# Import standard packages
import sys
import re
import os

# Treat modules in PACKAGE_PARENT as in working directory
if (__name__ == "__main__"):
    # Relative path for package to import
    PACKAGE_PARENT = '..'
    # Obtain path of currently-running file
    SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
    # Convert to absolute path, and add to the PYTHONPATH
    sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

    from algorithm import SingleSequenceAlgorithm
    from nucleotides import get_nt_freq, random_sequence
else:
    from .algorithm import SingleSequenceAlgorithm
    from ..nucleotides import get_nt_freq, random_sequence

class CRISPRater(SingleSequenceAlgorithm):
    def __init__(self):
        super().__init__(
            name="CRISPRater",
            authors=['Labuhn, Maurice', 'Adams, Felix F', 'Ng, Michelle', 'Knoess, Sabine', 'Schambach, Axel', 'Charpentier, Emmanuelle M', 'Schwarzer, Adrian', 'Mateo, Juan L', 'Klusmann, Jan-Henning', 'Heckl, Dirk'],
            title='Refined sgRNA efficacy prediction improves large- and small-scale CRISPR–Cas9 applications',
            journal='Nucleic Acids Research',
            issuing='46(3):1375-1385',
            year=2018,
            doi='https://doi.org/10.1093/nar/gkx1268',
            #citation="Labuhn, et al. Refined sgRNA efficacy prediction improves large- and small-scale CRISPR-Cas9 applications. Nucleic Acids Research, 46(3), (2018).",
            off_target=True, # The CRISPRater paper asserts this can be used for off-target scores
            on_target=True,
            prefilter=False,
            postfilter=False,
            minimum=1.0,
            maximum=100.0,
            default=100.0
        )
#        self.patternCG = re.compile("[CG]")
        # nt position 1 is the 5' end of the spacer, and PAM-distal
        # nt position 20 is the 3' end of the spacer, and PAM-proximal
        #                   [GC:4–13,    G:20,       TA:3,       GA:12       G:6,        TA:4,        GA:18.       CA:5,        G:14,        A:15]
        self.model_weight = [0.14177385, 0.06966514, 0.04216254, 0.03303432, 0.02355430, -0.04746424, -0.04878001, -0.06981921, -0.07087756, -0.08160700]
        self.model_offset = 0.6505037

    def weight(self, x):
        """Severely penalize any score less than 50"""
        return 1.0/(1+1.22**(66-x))

    def calculate(self, intended, *args, **kwargs):
        '''
        Uses fixed set of arguments.
        Returns output of self.score()
        '''
        sequence, target, pam, upstream, downstream = intended
        
        return self.score(upstream, target)
    
    def score(self, us, seq):
        """
        Uses arbitrary arguments.
        Performs scoring calculation.
        """
        # Right-justify sequence, and pad left with 'N' characters
        full_seq = ('N'*20+us+seq)[-20:]

        # Parse sequence
        features = [0]*10
        features[0] = get_nt_freq(['C', 'G'], full_seq[3:13])
        features[1] = get_nt_freq(['G'], full_seq[19])
        features[2] = get_nt_freq(['A', 'T'], full_seq[2])
        features[3] = get_nt_freq(['A', 'G'], full_seq[11])
        features[4] = get_nt_freq(['G'], full_seq[5])
        features[5] = get_nt_freq(['A', 'T'], full_seq[3])
        features[6] = get_nt_freq(['A', 'G'], full_seq[17])
        features[7] = get_nt_freq(['A', 'C'], full_seq[4])
        features[8] = get_nt_freq(['G'], full_seq[13])
        features[9] = get_nt_freq(['A'], full_seq[14])

        # Calculate score
        score = 0
        for idx in range(0, len(features)):
            score = score + features[idx]*self.model_weight[idx]
        score = score + self.model_offset

        # Return result
        return score*100
        #return self.getScore(('N'*20+us+seq)[-20:])*100

    
    ########################################################################
    # CRISPRater score calculation
    ########################################################################
    # Retrieved from
    #   https://bitbucket.org/juanlmateo/cctop_standalone/src/95ea199ba2b65963adecdc1a2bf555c9171bb622/CCTop.py
    
#    def getGCFreq(self, seq):
#        cg = len(self.patternCG.findall(seq))
#        return float(cg)/len(seq)

#    def calcFeatures(self, seq):
#        feat = [0]*10
#        feat[0] = get_nt_freq(['C', 'G'], seq[3:13])
#        feat[1] = get_nt_freq(['G'], seq[19])
#        feat[2] = get_nt_freq(['A', 'T'], seq[2])
#        feat[3] = get_nt_freq(['A', 'G'], seq[11])
#        feat[4] = get_nt_freq(['G'], seq[5])
#        feat[5] = get_nt_freq(['A', 'T'], seq[3])
#        feat[6] = get_nt_freq(['A', 'G'], seq[17])
#        feat[7] = get_nt_freq(['A', 'C'], seq[4])
#        feat[8] = get_nt_freq(['G'], seq[13])
#        feat[9] = get_nt_freq(['A'], seq[14])
#
#        return feat

#    def calcFeatures_old(self, seq):
#        feat = [0]*10
#        feat[0] = self.getGCFreq(seq[3:13])
#
#        if (seq[19] == "G"):
#            feat[1] = 1
#        if ((seq[2] == "T") or (seq[2] == "A")):
#            feat[2] = 1
#        if ((seq[11] == "G") or (seq[11] == "A")):
#            feat[3] = 1
#        if (seq[5] == "G"):
#            feat[4] = 1
#        if ((seq[3] == "T") or (seq[3] == "A")):
#            feat[5] = 1
#        if ((seq[17] == "G") or (seq[17] == "A")):
#            feat[6] = 1
#        if ((seq[4] == "C") or (seq[4] == "A")):
#            feat[7] = 1
#        if (seq[13] == "G"):
#            feat[8] = 1
#        if (seq[14] == "A"):
#            feat[9] = 1
#
#        return feat
    
#    def getScore(self, seq):
#        features = self.calcFeatures(seq)
#
#        score = 0
#        for idx in range(0, len(features)):
#            score = score + features[idx]*self.model_weight[idx]
#        score = score + self.model_offset
#
#        return score
    
#    def getScoreText(self, score):
#        '''
#        Round score and bin as LOW, MEDIUM, or HIGH
#        '''
#        text = "{0:.2f}".format(score)
#        if score < 0.56:
#            return text + ' (LOW)' #low score -> red
#        elif score > 0.74:
#            return text + ' (HIGH)' #high score -> blue
#        else:
#            return text + ' (MEDIUM)' #medium score -> grey
    
    ########################################################################
    

def test():
    data = [
        # Sequence, Target, PAM, US, DS
        ('', 'AAAATTAACTATAGGTAAAG', 'TGG', '', ''),
        ('', 'AACATCAACTCTAGCTAACG', 'CGG', '', ''),
        ('', 'AACATCAACTNNNGCTAACG', 'CGG', '', ''),
        ('',  'ACATCAACTCTAGCTAACG', 'CGG', '', ''), # Should pad with Ns
        ('',   'CATCAACTCTAGCTAACG', 'CGG', '', ''), # Should pad with Ns
        ('',    'ATCAACTCTAGCTAACG', 'CGG', '', ''), # Should pad with Ns
        ('',     'TCAACTCTAGCTAACG', 'CGG', '', ''), # Should pad with Ns
    ]
    
    print("=== CRISPRater ===")
    C = CRISPRater()
    for d in data:
        print(C.calculate(d))




    import time

    start = time.time()
    seqs = []
    scores = []
    for I in range(100000):
        d = ('', random_sequence(20), random_sequence(3), random_sequence(4), random_sequence(3))
        seqs.append(d)
        scores.append(C.calculate(d, disambiguate=True))

    print('min={}, max={}, mean={}, time={}'.format(min(scores), max(scores), sum(scores)/len(scores), time.time()-start))

    def figure(scores, C, x_min=0, x_max=100, title=''):
        import matplotlib.pyplot as plt

        fig, ax1 = plt.subplots(figsize=(8, 4))
        num_bins = 50
        n, bins, patches = ax1.hist(scores, num_bins, color='orange', histtype='step', label='Histogram') # The histogram
        ax2 = ax1.twinx()
        n, bins, patches = ax2.hist(scores, num_bins, density=True, histtype='step', cumulative=True, label='Cumulative') # plot the cumulative histogram

        x = list(range(x_min, x_max+1))
        y = []
        for i in x:
            y.append(C.weight(i))
        ax2.plot(x, y, 'k--', linewidth=1.5, label='Weight')

        ax2.grid(True)

        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper left')

        ax2.set_title(title)
        ax1.set_xlabel('Score')
        ax1.set_ylabel('Count')
        ax2.set_ylabel('Frequency')
        fig.tight_layout()
        plt.show()

    figure(scores, C, title='Random CRISPRater scores')

if (__name__ == "__main__"):
    test()
