#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/labuhn.py

# Import standard packages
import sys
import re

if (__name__ == "__main__"):
    from algorithm import SingleSequenceAlgorithm
else:
    from .algorithm import SingleSequenceAlgorithm

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
        self.patternCG = re.compile("[CG]")
        self.model_weight = [0.14177385, 0.06966514, 0.04216254, 0.03303432, 0.02355430, -0.04746424, -0.04878001, -0.06981921, -0.07087756, -0.08160700]
        self.model_offset = 0.6505037
    
    def calculate(self, intended, *args, **kwargs):
        sequence, target, pam, upstream, downstream = intended
        
        return self.score(upstream, target)
    
    def score(self, us, seq):
        """
        """
        
        return self.getScore((us+seq)[-20:])*100
    
    ########################################################################
    # CRISPRater score calculation
    ########################################################################
    # Retrieved from
    #   https://bitbucket.org/juanlmateo/cctop_standalone/src/95ea199ba2b65963adecdc1a2bf555c9171bb622/CCTop.py
    
    def getGCFreq(self, seq):
        cg = len(self.patternCG.findall(seq))
        return float(cg)/len(seq)
    
    def calcFeatures(self, seq):
        # Need to add ambiguous nt processing.
        feat = [0]*10
        feat[0] = self.getGCFreq(seq[3:13])
        if(seq[19] == "G"):
            feat[1] = 1
        if((seq[2] == "T") or (seq[2] == "A")):
            feat[2] = 1
        if((seq[11] == "G") or (seq[11] == "A")):
            feat[3] = 1
        if(seq[5] == "G"):
            feat[4] = 1
        if((seq[3] == "T") or (seq[3] == "A")):
            feat[5] = 1
        if((seq[17] == "G") or (seq[17] == "A")):
            feat[6] = 1
        if((seq[4] == "C") or (seq[4] == "A")):
            feat[7] = 1
        if(seq[13] == "G"):
            feat[8] = 1
        if(seq[14] == "A"):
            feat[9] = 1
        
        return(feat)
    
    def getScore(self, seq):
        features = self.calcFeatures(seq)
        
        score = 0
        for idx in range(0, len(features)):
            score = score + features[idx]*self.model_weight[idx]
        score = score + self.model_offset
        
        return(score)
    
    def getScoreText(self, score):
        '''
        Round score and bin as LOW, MEDIUM, or HIGH
        '''
        text = "{0:.2f}".format(score)
        if score < 0.56:
            return text + ' (LOW)' #low score -> red
        elif score > 0.74:
            return text + ' (HIGH)' #high score -> blue
        else:
            return text + ' (MEDIUM)' #medium score -> grey
    
    ########################################################################
    

def test():
    data = [
        ('', 'AAAATTAACTATAGGTAAAG', 'TGG', '', ''),
        ('', 'AACATCAACTCTAGCTAACG', 'CGG', '', '')
    ]
    
    print("=== CRISPRater ===")
    C = CRISPRater()
    for d in data:
        print(C.calculate(d))

if (__name__ == "__main__"):
    test()
