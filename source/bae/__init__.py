#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/bae/__init__.py

# Calculates the microhomology-associated 'out-of-frame' score
# that is positively correlated with the frequency of frame-shift mutations
# See Cas-OFFinder and Cas-Designer
# 
# Park, J, Bae, S and Kim, J-S (2015) Cas-Designer: a web-based tool for 
#  choice of CRISPR-Cas9 target sites. Bioinformatics, 31(24), 2015, 4014–4016.
#  doi: 10.1093/bioinformatics/btv537
# Bae,S. et al. (2014) Cas-OFFinder: A fast and versatile algorithm that
#  searches for potential off-target sites of Cas9 RNA-guided endonucleases.
#  Bioinformatics, 30, 1473–1475.
# Bae,S. et al. (2014) Microhomology-based choice of Cas9 nuclease target
#  sites. Nat. Methods, 11, 705–706.


def bae_score(seq):
    return 0