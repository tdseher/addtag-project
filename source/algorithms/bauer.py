#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/bauer.py

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


class Tuscan(SingleSequenceAlgorithm): # Need to pick an Algorithm subclass
    # Paper:
    #  https://www.liebertpub.com/doi/full/10.1089/crispr.2017.0021
    #  High Activity Target-Site Identification Using Phenotypic Independent CRISPR-Cas9 Core Functionality
    #  Laurence O.W. Wilson, Daniel Reti, Aidan R. O'Brien, Robert A. Dunne, and Denis C. Bauer
    #  The CRISPR Journal Vol. 1, No. 2
    #  1 Apr 2018 (https://doi.org/10.1089/crispr.2017.0021)

    # Source code to integrate:
    #  https://github.com/BauerLab/TUSCAN

    def __init__(self):
        super().__init__("TUSCAN", "Wilson, et al", 2018,
            citation="Wilson, et al. High Activity Target-Site Identification Using Phenotypic Independent CRISPR-Cas9 Core Functionality. The CRISPR Journal 1(2) (2018)",
            off_target=False,
            on_target=True,
            prefilter=False,
            postfilter=False,
            minimum=1.0,
            maximum=100.0,
            default=None
        )

class Varscot(SingleSequenceAlgorithm): # Need to pick an Algorithm subclass
    # Paper:
    #  https://bmcbiotechnol.biomedcentral.com/articles/10.1186/s12896-019-0535-5

    # Source code to integrate:
    #  https://github.com/BauerLab/VARSCOT
    
    def __init__(self):
        super().__init__("VARSCOT", "Wilson, et al", 2018,
            citation=("Wilson, et al. VARSCOT: variant-aware detection and scoring enables sensitive and "
                      "personalized off-target detection for CRISPR-Cas9. BMC Biotechnologyvolume 19, "
                      "Article number: 40 (2019)"),
            off_target=True,
            on_target=False,
            prefilter=False,
            postfilter=False,
            minimum=1.0,
            maximum=100.0,
            default=None
        )

