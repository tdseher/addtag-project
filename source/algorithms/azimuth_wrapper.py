#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/azimuth_wrapper.py

# Import standard packages
from __future__ import print_function
import sys
import os
import platform
import re

# Import non-standard packages
#import regex

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

PROGRAM = os.path.basename(sys.argv[0])
AUTHOR = 'Thaddeus D. Seher (@tdseher)'
DATE = '2017'
COPYRIGHT = ('copyright:' '\n'
             '  {PROGRAM} Copyright (c) {DATE} {AUTHOR}. All rights reserved.' '\n'
            ).format(**globals())
DESCRIPTION = ('description:' '\n'
               '  This is the description' '\n'
              )
USAGE = ('usage: {PROGRAM} [-h] [-i *.seqs] [-o *.scores] [SEQ [SEQ ...]]' '\n'
        ).format(**globals())

POSITIONAL_ARGS = ('STDIN/positional arguments:' '\n'
                   '  SEQ                             Sequence to analyze' '\n'
                   )
OPTIONAL_ARGS = ('optional arguments:' '\n'
                 '  -h, --help                      Show this help message and exit' '\n'
                 '  -i *.seqs, --input *.seqs       Path to file containing sequences to analyze' '\n'
                 '  -o *.scores, --output *.scores  Send output analysis to this file instead of STDOUT' '\n'
                 )
EPILOG = ('examples:' '\n'
          '  echo ACAGCTGATCTCCAGATATGACCATGGGTT CAGCTGATCTCCAGATATGACCATGGGTTT | {PROGRAM} -o output.txt' '\n'
          '  {PROGRAM} CCAGAAGTTTGAGCCACAAACCCATGGTCA ATTAACCCCCGATAAGGCTGAATCGATAGG > output.txt' '\n'
          '  {PROGRAM} GGCTTAGAGCTAGCTTTAGCTCTAGAGAAA ACCCCGTAAGGGCTCGCCTAGGATCCAAAA -i input1.txt -i input2.txt -o output.txt' '\n'
          ).format(**globals())


def parse_arguments():
    # Check arguments for help flag
    for pp in sys.argv[1:]:
        if pp in ['-h', '--help']:
            print(USAGE)
            print(DESCRIPTION)
            print(COPYRIGHT)
            print(POSITIONAL_ARGS)
            print(OPTIONAL_ARGS)
            print(EPILOG)
            sys.exit(1)
    
    # Make empty list to store sequences
    seqs = []
    
    # Add sequences from STDIN
    #seqs += re.split(r'\s+', sys.stdin.read().strip())
    
    if not sys.stdin.isatty():
    #if not os.isatty(sys.stdin.fileno()):
        seqs += re.split(r'\s+', sys.stdin.read().strip())
    elif (len(sys.argv[1:]) == 0):
        print(USAGE)
        sys.exit(1)
    
    # If a positional parameter has '-i' or '--input', then treat it as an input file
    # If a positional parameter has '-o' or '--output', then treat it as the output file
    # Otherwise, treat the positional parameter as a sequence
    infiles = []
    outfile = sys.stdout
    ipp = False
    opp = False
    for pp in sys.argv[1:]:
        if ipp:
            infiles.append(pp)
            ipp = False
        elif opp:
            outfile = open(pp, 'w')
            opp = False
        else:
            if pp in ['-i', '--input']:
                ipp = True
            elif pp in ['-o', '--output']:
                opp = True
            else:
                seqs.append(pp)
    
    # Add sequences from the files
    for f in infiles:
        with open(f, 'r') as flo:
            for line in flo:
                seqs += re.split(r'\s+', line.strip())
    
    #print('OUTPUT', file=outfile)
    #print(seqs, file=outfile)
    
    return seqs, outfile

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

def chunk(l, n):
    """
    Generator to split list 'l' into chunks of at most size 'n'
    """
    for i in range(0, len(l), n):  
        yield l[i:i + n]

def main():
    seqs, outfile = parse_arguments()
    
    # All nucleotide sequences must have length of 30
    # Otherwise, exit the script with an error code
    for seq in seqs:
        if (len(seq) != 30):
            sys.exit(1)
    
    # Benchmarks for 20,000 random sequences
    # ver  chunks  chunksize  outlines  time(s)
    # 3.6   1      20000      20001     129.6287431716919
    # 3.6   2      10000      20002     105.1144630908966
    # 3.6   3       6667      20003      90.43231105804443
    # 3.6   4       5000      20004      94.43034148216248
    # 3.6   5       4000      20005      97.77741646766663
    # 3.6   6       3334      20006      97.67567348480225
    # 3.6   7       2858      20007      99.47095489501953
    # 3.6   8       2500      20008     108.14812755584717
    # 3.6   9       2223      20009     119.17225742340088
    # 3.6  10       2000      20010     121.30107951164246
    #
    # 2.7   1      20000      20001     144.61829352378845
    # 2.7   2      10000      20002     108.74389004707336
    # 2.7   3       6667      20003      93.69686198234558
    # 2.7   4       5000      20004      98.56164050102234
    # 2.7   5       4000      20005      98.97256851196289
    # 2.7   6       3334      20006      92.09488368034363
    # 2.7   7       2858      20007      92.62291359901428
    # 2.7   8       2500      20008      94.20428347587585
    # 2.7   9       2223      20009      94.46121072769165
    # 2.7  10       2000      20010      97.0748233795166    
    
    chunk_size = 6000 # Maximum number of sequences to put into Azimuth at a time
    
    for s in chunk(seqs, chunk_size):
        # Convert to numpy array
        sequences = array(s)
        
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
            print(sequences[i] + ' ' + str(prediction), file=outfile)

if (__name__ == '__main__'):
    main()
