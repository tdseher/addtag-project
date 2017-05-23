#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/distance.py

# Import standard packages
import sys
import os

# Import non-standard packages

# Treat modules in PACKAGE_PARENT as in working directory
if (__name__ == "__main__"):
    # Relative path for package to import
    PACKAGE_PARENT = '..'
    # Obtain path of currently-running file
    SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
    # Convert to absolute path, and add to the PYTHONPATH
    sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))
    
    from algorithm import PairedSequenceAlgorithm
    from nucleotides import rc
else:
    from .algorithm import PairedSequenceAlgorithm
    from ..nucleotides import rc

class Substitutions(PairedSequenceAlgorithm):
    def __init__(self):
        super().__init__("Substitutions", "Seher", 2017,
            citation="AddTag",
            off_target=False,
            prefilter=False,
            postfilter=False,
            minimum=0.0,
            maximum=5.0,
            default=0.0
        )
    
    def calculate(self, intended, potential, *args, **kwargs):
        on_sequence, on_target, on_pam, on_upstream, on_downstream = intended
        off_sequence, off_target, off_pam, off_upstream, off_downstream = potential
        
        return self.score(on_target, off_target)
    
    def score(self, seq1, seq2):
        """"""
        return 0.0

class Insertions(PairedSequenceAlgorithm):
    def __init__(self):
        super().__init__("Insertions", "Seher", 2017,
            citation="AddTag",
            off_target=False,
            prefilter=False,
            postfilter=False,
            minimum=0.0,
            maximum=2.0,
            default=0.0
        )
    
    def calculate(self, intended, potential, *args, **kwargs):
        on_sequence, on_target, on_pam, on_upstream, on_downstream = intended
        off_sequence, off_target, off_pam, off_upstream, off_downstream = potential
        
        return self.score(on_target, off_target)
    
    def score(self, seq1, seq2):
        """"""
        return 0.0

class Deletions(PairedSequenceAlgorithm):
    def __init__(self):
        super().__init__("Deletions", "Seher", 2017,
            citation="AddTag",
            off_target=False,
            prefilter=False,
            postfilter=False,
            minimum=0.0,
            maximum=2.0,
            default=0.0
        )
    
    def calculate(self, intended, potential, *args, **kwargs):
        on_sequence, on_target, on_pam, on_upstream, on_downstream = intended
        off_sequence, off_target, off_pam, off_upstream, off_downstream = potential
        
        return self.score(on_target, off_target)
    
    def score(self, seq1, seq2):
        """"""
        return 0.0

class Errors(PairedSequenceAlgorithm):
    def __init__(self):
        super().__init__("Errors", "Seher", 2017,
            citation="AddTag",
            off_target=False,
            prefilter=False,
            postfilter=False,
            minimum=0.0,
            maximum=5.0,
            default=0.0
        )
    
    def calculate(self, intended, potential, *args, **kwargs):
        on_sequence, on_target, on_pam, on_upstream, on_downstream = intended
        off_sequence, off_target, off_pam, off_upstream, off_downstream = potential
        
        return self.score(on_target, off_target)
    
    def score(self, seq1, seq2):
        """"""
        return 0.0

def test():
    """Code to test the functions and classes"""
    
    a = ('', 'CGATGGCTAGGATCGATTGA', 'TGG', '', '')
    b = ('', 'RYMKWSACGTbDHVNACGTA', 'TGG', '', '')
    c = ('',   'ATGSCTCGGATCGATTGA', 'AGG', '', '')
    
    print("=== Substitutions ===")
    C = Substitutions()
    print(C.calculate(a, b))
    print(C.calculate(a, c))
    
    print("=== Insertions ===")
    C = Insertions()
    print(C.calculate(a, b))
    print(C.calculate(a, c))
    
    print("=== Deletions ===")
    C = Deletions()
    print(C.calculate(a, b))
    print(C.calculate(a, c))
    
    print("=== Errors ===")
    C = Errors()
    print(C.calculate(a, b))
    print(C.calculate(a, c))

if (__name__ == "__main__"):
    test()