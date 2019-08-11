#!/usr/bin/env python

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/azimuth_wrapper.py

# Import standard packages
import sys
import platform

# Import non-standard packages
if (platform.system() == 'Darwin'): # 'Darwin', 'Linux', 'Windows', 'SunOs'
    # Code to fix Mac OS X errors/warnings, especially when freezing/packing
    
    # Explicitly set Matplotlib rendering engine to 'Agg' because Darwin
    # platforms have trouble selecting an appropriate one automatically
    import matplotlib
    matplotlib.use('Agg')
    
    # In Darwin platforms, Matplotlib has a negligible font error/warning
    # that should be ignored. Importing 'azimuth.model_comparison' will
    # trigger it, so it is wrapped in a 'catch_warnings()' context manager
    # to suppress the warnings.
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import azimuth.model_comparison
else:
    import azimuth.model_comparison
from numpy import array

# Raise a useful error if utilizing an incompatible version of Python
#if ((sys.version_info < (2, 7, 10)) or (sys.version_info >= (3, 0, 0))):
#    raise Exception("AddTag Azimuth wrapper requires 3.0.0 > Python >= 2.7.10. You are running with " + ".".join(map(str, sys.version_info[:3])) + ".")

def test():
    sequences = array(['ACAGCTGATCTCCAGATATGACCATGGGTT', 'CAGCTGATCTCCAGATATGACCATGGGTTT', 'CCAGAAGTTTGAGCCACAAACCCATGGTCA'])
    amino_acid_cut_positions = array([2, 2, 4])
    percent_peptides = array([0.18, 0.18, 0.35])
    predictions = azimuth.model_comparison.predict(sequences, amino_acid_cut_positions, percent_peptides)
    
    for i, prediction in enumerate(predictions):
        print(sequences[i] + ' ' + str(prediction))
    
    
    azimuth.model_comparison.predict(array(['GGGAGGCTGCTTTACCCGCTGTGGGGGCGC']), array([-1]), array([-1]))
    # No model file specified, using V3_model_nopos
    # array([ 0.58477438])

def main():
    # All nucleotides must have length of 30
    # Otherwise, exit the script with an error code
    for arg in sys.argv[1:]:
        if (len(arg) != 30):
            sys.exit(1)
    
    sequences = array(sys.argv[1:])
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
    predictions = azimuth.model_comparison.predict(sequences, aa_cut=None, percent_peptide=None, model=None, pam_audit=False)
    #predictions = azimuth.model_comparison.predict(sequences, amino_acid_cut_positions, percent_peptides)
    
    for i, prediction in enumerate(predictions):
        print(sequences[i] + ' ' + str(prediction))

if (__name__ == '__main__'):
    main()
