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
    pass

class Insertions(PairedSequenceAlgorithm):
    pass

class Deletions(PairedSequenceAlgorithm):
    pass

class Errors(PairedSequenceAlgorithm):
    pass