#!/usr/bin/env python

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/azimuth/azimuth.py

import sys

import azimuth.model_comparison
import numpy as np

# Raise a useful error if utilizing an incompatible version of Python
if ((sys.version_info < (2, 7, 12)) or (sys.version_info >= (3, 0, 0))):
    raise Exception("AddTag Azimuth wrapper requires 3.0.0 > Python >= 2.7.12. You are running with " + ".".join(map(str, sys.version_info[:3])) + ".")

def main():
    
    sequences = np.array(['ACAGCTGATCTCCAGATATGACCATGGGTT', 'CAGCTGATCTCCAGATATGACCATGGGTTT', 'CCAGAAGTTTGAGCCACAAACCCATGGTCA'])
    amino_acid_cut_positions = np.array([2, 2, 4])
    percent_peptides = np.array([0.18, 0.18, 0.35])
    predictions = azimuth.model_comparison.predict(sequences, amino_acid_cut_positions, percent_peptides)
    
    for i, prediction in enumerate(predictions):
        print sequences[i], prediction
    
    
    azimuth.model_comparison.predict(np.array(['GGGAGGCTGCTTTACCCGCTGTGGGGGCGC']), np.array([-1]), np.array([-1]))
    # No model file specified, using V3_model_nopos
    # array([ 0.58477438])

if (__name__ == '__main__'):
    main()
