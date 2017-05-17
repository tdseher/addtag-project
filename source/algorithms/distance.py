#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/distance.py

# Import non-standard packages
if (__name__ == "__main__"):
    from algorithm import PairedSequenceAlgorithm
    from nucleotides import rc # to change
else:
    from .algorithm import PairedSequenceAlgorithm
    from .nucleotides import rc # to change

class Substitutions(PairedSequenceAlgorithm):
    pass

class Insertions(PairedSequenceAlgorithm):
    pass

class Deletions(PairedSequenceAlgorithm):
    pass

class Errors(PairedSequenceAlgorithm):
    pass