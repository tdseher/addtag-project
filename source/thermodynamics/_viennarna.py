#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/thermodynamics/_viennarna.py

# Import standard packages
import sys
import os
import logging
import math

# Import non-standard packages
import regex

logger = logging.getLogger(__name__)

# import included AddTag-specific modules

# Treat modules in PACKAGE_PARENT as in working directory
if (__name__ == "__main__"):
    # Relative path for package to import
    PACKAGE_PARENT = '..'
    # Obtain path of currently-running file
    SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
    # Convert to absolute path, and add to the PYTHONPATH
    sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))
    
    from nucleotides import rc
    from oligo import Oligo, Structure
    from utils import which
else:
    from ..nucleotides import rc
    from ..utils import which
    from .oligo import Oligo, Structure

class ViennaRNA(Oligo):
    logger = logger.getChild(__qualname__)
    
    loaded = False
    def __init__(self):
        super().__init__("ViennaRNA", "", 0,
            citation=""
        )
    
    @classmethod
    def load(cls):
        global RNA
        
        import RNA
        # If ViennaRNA is not installed correctly, then this error will be raised:
        #   ModuleNotFoundError: No module named 'RNA'
        
        try:
            RNA.params_load_DNA_Mathews2004()
            # If the ViennaRNA version does not support DNA parameters, then this error will be raised:
            #   AttributeError: module 'RNA' has no attribute 'params_load_DNA_Mathews2004'
        except AttributeError:
            # We try to load the DNA parameters manually
            
            # A list to hold the most-likely paths of the DNA parameters file
            dna_paths = []
            
            dna_paths.append(
                os.path.abspath(os.path.join(os.path.dirname(RNA.__file__), '..', '..', '..', '..', 'share', 'ViennaRNA', 'dna_mathews2004.par'))
            )
            
            prog_path = which('RNAfold')
            if prog_path:
                dna_paths.append(
                    os.path.abspath(os.path.join(os.path.dirname(prog_path), 'Misc', 'dna_mathews2004.par'))
                )
            
            for dp in dna_paths:
                if os.path.exists(dp):
                    RNA.params_load(dp)
                    break
            else:
                raise FileNotFoundError('Cannot find file encoding DNA thermodynamics parameters required by ViennaRNA: {}'.format(repr(dna_paths)))
        cls.loaded = True
    
    @classmethod
    def find_structures(cls, folder, seq1, seq2=None, sodium=0.05, magnesium=0.0, temperature=25, concentration=0.00000025, **kwargs):
        """
        Should return the list of 'Structure' objects with delta-G, deltaH, deltaS, and Tm values.
        
        Accepts 1 or 2 input sequences. Automatically runs either:
         * Hairpin     (1 input sequence: A=seq1, UNAFold run on A)
         * Homodimer   (2 identical input sequences: A=seq1=seq2, UNAFold run on A & A)
         * Heterodimer (2 input sequences: A=seq1 B=seq2, UNAFold run on A & B)
        """
        if not cls.loaded:
            cls.load()
        
        mfe = None
        if (seq1 == seq2): # Homodimer calculation
            duplex = RNA.duplexfold(seq1, seq1)
            mfe = duplex.energy
        elif (seq2 == None): # Hairpin calculation
            ss, mfe = RNA.fold(seq1)
        else: # Heterodimer calculation, Tm calculation [seq1, rc(seq1)]
            duplex = RNA.duplexfold(seq1, seq2)
            mfe = duplex.energy
        
        if mfe != None:
            s = Structure(seq1, seq2, mfe, None, None, None, sodium, magnesium, temperature, concentration)
        else:
            s = Structure(seq1, seq2, math.inf, math.inf, math.inf, math.inf, sodium, magnesium, temperature, concentration)
        
        return [s]

def test():
    """Code to test the classes and functions in 'source/oligos/_unafold.py'"""
    
    C = ViennaRNA()
    print("===", C.name, "===")
    
    a = 'GAAATCGCTTAGCGCGAACTCAGACCAT'
    b = 'CCTAGCTATTTAATAAATC'
    c = 'TTCTCCACTTCCATCACCGT'
    
    tempdir = '.'
    
    print('Hairpin: {}'.format(repr(a)))
    for s in C.find_structures(tempdir, a, None):
        print('', s)
    print('Homodimer: {} {}'.format(repr(a), repr(a)))
    for s in C.find_structures(tempdir, a, a):
        print('', s)
    print('Heterodimer: {} {}'.format(repr(a), repr(b)))
    for s in C.find_structures(tempdir, a, b):
        print('', s)
    print('Reverse-complements: {} {}'.format(repr(a), repr(rc(a))))
    for s in C.find_structures(tempdir, a, rc(a)):
        print('', s)

if (__name__ == "__main__"):
    test()
