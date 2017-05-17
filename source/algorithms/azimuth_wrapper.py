#!/usr/bin/env python

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/azimuth_wrapper.py

# Import standard packages
import sys

# Import non-standard packages
import azimuth.model_comparison
import numpy as np

# Raise a useful error if utilizing an incompatible version of Python
if ((sys.version_info < (2, 7, 10)) or (sys.version_info >= (3, 0, 0))):
    raise Exception("AddTag Azimuth wrapper requires 3.0.0 > Python >= 2.7.10. You are running with " + ".".join(map(str, sys.version_info[:3])) + ".")

def test():
    sequences = np.array(['ACAGCTGATCTCCAGATATGACCATGGGTT', 'CAGCTGATCTCCAGATATGACCATGGGTTT', 'CCAGAAGTTTGAGCCACAAACCCATGGTCA'])
    amino_acid_cut_positions = np.array([2, 2, 4])
    percent_peptides = np.array([0.18, 0.18, 0.35])
    predictions = azimuth.model_comparison.predict(sequences, amino_acid_cut_positions, percent_peptides)
    
    for i, prediction in enumerate(predictions):
        print sequences[i], prediction
    
    
    azimuth.model_comparison.predict(np.array(['GGGAGGCTGCTTTACCCGCTGTGGGGGCGC']), np.array([-1]), np.array([-1]))
    # No model file specified, using V3_model_nopos
    # array([ 0.58477438])

def main():
    # All nucleotides must have length of 30
    for arg in sys.argv[1:]:
        if (len(arg) != 30):
            exit(1)
    
    sequences = np.array(sys.argv[1:])
    #amino_acid_cut_positions = np.array([-1]*len(sys.argv[1:]))
    #percent_peptides = np.array([-1]*len(sys.argv[1:]))
    # Help on function predict in module azimuth.model_comparison:
    # predict(seq, aa_cut=None, percent_peptide=None, model=None, model_file=None, pam_audit=True, length_audit=False, learn_options_override=None)
    # Args:
    #     seq: numpy array of 30 nt sequences.
    #     aa_cut: numpy array of amino acid cut positions (optional).
    #     percent_peptide: numpy array of percent peptide (optional).
    #     model: model instance to use for prediction (optional).
    #     model_file: file name of pickled model to use for prediction (optional).
    #     pam_audit: check PAM of each sequence.
    #     length_audit: check length of each sequence.
    #     learn_options_override: a dictionary indicating which learn_options to override (optional).
    # Returns: a numpy array of predictions.
    predictions = azimuth.model_comparison.predict(sequences, pam_audit=False)
    #predictions = azimuth.model_comparison.predict(sequences, amino_acid_cut_positions, percent_peptides)
    
    for i, prediction in enumerate(predictions):
        print sequences[i], prediction

if (__name__ == '__main__'):
    main()
