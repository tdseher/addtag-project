#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/thermodynamics/_viennarna.py

# Import standard packages
import sys
import os
import logging
import math
import subprocess
from collections import OrderedDict
import importlib.util

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
    parameters = None
    
    def __init__(self):
        super().__init__(
            name='ViennaRNA',
            authors=['Lorenz, Ronny', 'Bernhart, Stephan H.', 'Höner zu Siederdissen, Christian', 'Tafer, Hakim', 'Flamm, Christoph', 'Stadler, Peter F.', 'Hofacker, Ivo L.'],
            title='ViennaRNA Package 2.0',
            journal='Algorithms for Molecular Biology',
            issuing='6(1):26',
            year=2011,
            doi='https://doi.org/10.1186/1748-7188-6-26'
            #citation="Lorenz, et al. ViennaRNA Package 2.0. Algorithms for Molecular Biology 6:1, 26 (2011)."
        )
    
    def is_available(self):
        """
        Determines if the prerequisites for ViennaRNA have been met.
        :return: True or False
        """
        spec = importlib.util.find_spec('RNA')
        if spec:
            #self.available = True
            return True
        else:
            return False
    
    @classmethod
    def parameter_path(cls):
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
                return dp
        
        raise FileNotFoundError('Cannot find file encoding DNA thermodynamics parameters required by ViennaRNA: {!r}'.format(dna_paths))
    
    @classmethod
    def load(cls):
        if not cls.loaded:
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
                
                cls.parameters = cls.parameter_path()
                RNA.params_load(cls.parameters)
                
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
    
    @classmethod
    def find_tms(cls, sequences, sodium=0.05, magnesium=0.0, temperature=25, concentration=0.00000025, **kwargs):
        cls.load()
        
        def flatten(iterable, remove_none=False, add_equal=False):
            """Make a flat list out of a list of lists"""
            # Will remove None
            if remove_none:
                if add_equal:
                    return ['{}={}'.format(par, val) if (val != None) else par for (par, val) in iterable]
                else:
                    return [item for sublist in iterable for item in sublist if item != None]
            else:
                if add_equal:
                    return ['{}={}'.format(par, val) for (par, val) in iterable]
                else:
                    return [item for sublist in iterable for item in sublist]

        # RNAplex in 'probe mode' only calculates the reverse-complement Tm, and cannot
        # calculate the hairpin, homodimer, or heterodimer Tms.
        options = OrderedDict([
            ('--paramFile', cls.parameters),
            ('--probe-mode', None),
            ('--probe-concentration', Oligo.float_to_str(concentration)),
            ('--na-concentration', Oligo.float_to_str(sodium)),
            ('--mg-concentration', Oligo.float_to_str(magnesium)),
            ('--tris-concentration', Oligo.float_to_str(0.0)),
            ('--k-concentration', Oligo.float_to_str(0.0)),
            ('--temp', temperature),
        ])
        flat_options = flatten(options.items(), remove_none=True, add_equal=True)
        command_list = ['RNAplex'] + list(map(str, flat_options))
        command_str = ' '.join(command_list)
        
        cls.logger.info('command: {!r}'.format(command_str))
        
        cp = subprocess.run(command_list, input=bytes('\n'.join(sequences), 'utf-8'), shell=False, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        # The output typically looks like this:
        #   Probe mode
        #   Concentration K:0.000 TNP:0.000 Mg:0.000 Na:0.050 probe:0.000
        #   
        #                                                                                               sequence  DDSL98  DDSL04  DRSU95  RRXI98 CURRENT
        #                                                                                   AGGCTTTAGGGCTATAGGAA   51.78  50.76   50.74   62.06   52.13
        #                                                                                    CGAATTTAGAGCCTATAAT   43.60  42.90   39.52   48.07   43.46
        #                                                                                      GGCTATGAGATAGCTAA   43.14  42.68   41.24   52.81   43.44
        
        out_lines = cp.stdout.decode().splitlines()
        data_found = False
        
        tm_list = []
        for line in out_lines:
            line = line.rstrip()
            if not data_found:
                m = regex.search(r'^\s+sequence\s+DDSL98\s+DDSL04\s+DRSU95\s+RRXI98\s+CURRENT', line)
                if m:
                    data_found = True
            else:
                m = regex.match(r'^\s*(\S+)(?:\s+(\S+)){5}$', line)
                if m:
                    seq = m.captures(1)[0]
                    tm = float(m.captures(2)[-1]) # Use the 'CURRENT' column
                    tm_list.append(tm)
        
        return tm_list

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

    print('Reverse-complement Tms:')
    seqs = [a, b, c]
    tm_list = C.find_tms(seqs)
    for s, tm in zip(seqs, tm_list):
        print('', s, tm)
    
    # Expected output:
    #   === ViennaRNA ===
    #   Hairpin: 'GAAATCGCTTAGCGCGAACTCAGACCAT'
    #    Structure(dG=-3.5, dH=None, dS=None, Tm=None)
    #   Homodimer: 'GAAATCGCTTAGCGCGAACTCAGACCAT' 'GAAATCGCTTAGCGCGAACTCAGACCAT'
    #    Structure(dG=-10.399999618530273, dH=None, dS=None, Tm=None)
    #   Heterodimer: 'GAAATCGCTTAGCGCGAACTCAGACCAT' 'CCTAGCTATTTAATAAATC'
    #    Structure(dG=-4.800000190734863, dH=None, dS=None, Tm=None)
    #   Reverse-complements: 'GAAATCGCTTAGCGCGAACTCAGACCAT' 'ATGGTCTGAGTTCGCGCTAAGCGATTTC'
    #    Structure(dG=-37.900001525878906, dH=None, dS=None, Tm=None)
    #   Tms:
    #    GAAATCGCTTAGCGCGAACTCAGACCAT 63.51
    #    CCTAGCTATTTAATAAATC 39.73
    #    TTCTCCACTTCCATCACCGT 54.09

if (__name__ == "__main__"):
    test()
