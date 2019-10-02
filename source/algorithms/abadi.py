#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/abadi.py

# Import standard packages
import sys
import os
import logging

logger = logging.getLogger(__name__)

# Import non-standard packages
import regex

# Treat modules in PACKAGE_PARENT as in working directory
if (__name__ == "__main__"):
    # Relative path for package to import
    PACKAGE_PARENT = '..'
    # Obtain path of currently-running file
    SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
    # Convert to absolute path, and add to the PYTHONPATH
    sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

    from algorithm import SingleSequenceAlgorithm
    #from nucleotides import rc, disambiguate_iupac, random_sequence
    #from utils import which
else:
    from .algorithm import SingleSequenceAlgorithm
    #from ..nucleotides import rc, disambiguate_iupac
    #from ..utils import which

# Paper:
#  https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005807
#  A machine learning approach for predicting CRISPR-Cas9 cleavage efficiencies and patterns underlying its mechanism of action
#  Shiran Abadi, Winston X. Yan, David Amar, Itay Mayrose
#  PLOS Computational Biology. 2017; 13(10): e1005807.
#  October 16, 2017 (https://doi.org/10.1371/journal.pcbi.1005807)

# Source code to integrate:
#  https://www.tau.ac.il/~itaymay/software.html
#  http://crista.tau.ac.il/download.html
#  http://crista.tau.ac.il/CRISTA.zip

# A copy of the code is available here:
#  https://github.com/bm2-lab/CRISPR-off-target-data-imbalance
#  https://github.com/bm2-lab/CRISPR-off-target-data-imbalance/tree/master/scripts_for_improve_CRISTA/detailed
# Related repo:
#  https://github.com/MichaelChuai/CRISTA-comments

class Crista(SingleSequenceAlgorithm):
    def __init__(self):
        super().__init__("CRISTA", "Abadi, et al", 2017,
            citation=("Abadi, et al. A machine learning approach for predicting CRISPR-Cas9 cleavage efficiencies "
                      "and patterns underlying its mechanism of action. "
                      "PLOS Computational Biology 13(10): e1005807 (2017)."),
            off_target=False,
            on_target=True,
            prefilter=False,
            postfilter=False,
            minimum=1.0,
            maximum=100.0,
            default=None
        )

