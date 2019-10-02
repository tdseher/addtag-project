#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/wong.py

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
#  https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0784-0
#  WU-CRISPR: characteristics of functional guide RNAs for the CRISPR/Cas9 system
#  Nathan Wong, Weijun Liu, and Xiaowei Wang
#  Genome Biology volume 16, Article number: 218 (2015)
#  (https://doi.org/10.1186/s13059-015-0784-0)

# Source code to integrate:
#  https://github.com/wang-lab/sgDesigner   <-- Which of these is the algorithm?
#  https://github.com/wang-lab/WU-CRISPR    <--/

class SgDesigner(SingleSequenceAlgorithm):
    def __init__(self):
        super().__init__("SgDesigner", "Wong, et al", 2018,
            citation=("Wong, et al. WU-CRISPR: characteristics of functional guide RNAs for the CRISPR/Cas9 system. "
                      "Genome Biology 16, Article number 218 (2015)"),
            off_target=False,
            on_target=True,
            prefilter=False,
            postfilter=False,
            minimum=1.0,
            maximum=100.0,
            default=None
        )

