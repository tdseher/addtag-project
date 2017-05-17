#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/aligners/__init__.py

# Import standard packages
import os
from importlib import import_module

# Import included AddTag-specific modules
from .aligner import Aligner

# Import all Algorithm subclasses defined in python files within this same folder
exclusions = [
    '__init__.py',
    'blastplus.py',
    'blat.py',
    'bowtie.py',
    'bwa.py',
    'cctop.py',
]

path = os.path.dirname(os.path.abspath(__file__))
files = [f.rstrip(".py") for f in os.listdir(path) if (f.endswith('.py') and (f not in exclusions))]
for f in files:
    module = import_module('.'.join([__name__, f]))

# Create an instance of each SingleSequenceAlgorithm subclass, and add to this list
aligners = []
for C in Aligner.__subclasses__():
    aligners.append(C())