#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/chuai.py

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
#  https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1459-4
#  DeepCRISPR: optimized CRISPR guide RNA design by deep learning
#  Guohui Chuai, Hanhui Ma, Jifang Yan, Ming Chen, Nanfang Hong, Dongyu Xue, Chi Zhou, Chenyu Zhu, Ke Chen, Bin Duan, Feng Gu, Sheng Qu, Deshuang Huang, Jia Wei, and Qi Liu
#  Genome Biology. 2018;19(1):80. doi: 10.1186/s13059-018-1459-4.

# Source code to integrate:
#  http://www.deepcrispr.net/
#  https://github.com/bm2-lab/DeepCRISPR


class DeepCrispr(SingleSequenceAlgorithm):
    def __init__(self):
        super().__init__(
            name="DeepCRISPR",
            authors=['Chuai, Guohui', 'Ma, Hanhui', 'Yan, Jifang', 'Chen, Ming', 'Hong, Nanfang', 'Xue, Dongyu', 'Zhou, Chi', 'Zhu, Chenyu', 'Chen, Ke', 'Duan, Bin', 'Gu, Feng', 'Qu, Sheng', 'Huang, Deshuang', 'Wei, Jia', 'Liu, Qi'],
            title='DeepCRISPR: optimized CRISPR guide RNA design by deep learning',
            journal='Genome Biology',
            issuing='19(1):80',
            year=2018,
            doi='https://doi.org/10.1186/s13059-018-1459-4',
            #citation=("Chuai, et al. DeepCRISPR: optimized CRISPR guide RNA design by deep learning. "
            #          "Genome Biology. 19(1):80 (2018)."),
            off_target=False,
            on_target=True,
            prefilter=False,
            postfilter=False,
            minimum=1.0,
            maximum=100.0,
            default=None
        )

