#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/chari/__init__.py

# Metric for "efficiency"
#  Chari
#   Range: 0-100. Support Vector Machine, converted to rank-percent,
#   trained on data from 1235 guides targeting sequences that were also
#   transfected with a lentivirus into human 293T cells. See Chari et al.:
#   http://www.nature.com/nmeth/journal/v12/n9/abs/nmeth.3473.html
#   https://www.nature.com/articles/nmeth.3473
#   Chari R, Mali P, Moosburner M, et al. Unraveling CRISPR-Cas9 genome engineering parameters via a library-on-library approach. Nat Methods. 2015;12:823–826. DOI: 10.1038/nmeth.3473

# Requires support vector machine program 'svm_classify'
# see:
#   crisporWebsite-5dc1de3\crisporEffScores.py
#   crisporWebsite-5dc1de3\bin\src\svm_light
def calcChariScores(seqs):
    """ return dict with chari 2015 scores, returns two lists (rawScores, rankPercent)
    input seqs have lengths 21bp: 20 bp guide + 1bp first from PAM
    >>> calcChariScores(["CTTCTTCAAGGTAACTGCAGA", "CTTCTTCAAGGTAACTGGGGG"])
    ([0.54947621, 0.58604487], [80, 81])
    >>> calcChariScores(["CTTCTTCAAGGNAACTGCAGA"])
    ([0.9025848], [88])
    """
    pass