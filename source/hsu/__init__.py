#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/hsu/__init__.py

# Off-target score
#  Hsu specificity score
#   The higher the specificity score, the lower are off-target effects in the genome.
#   The specificity score ranges from 0-100 and measures the uniqueness
#   of a guide in the genome. See Hsu et al. Nat Biotech 2013.
#   (http://dx.doi.org/10.1038/nbt.2647) We recommend values >50, where possible.

def off_target_score():
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