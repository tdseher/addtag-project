#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/oligos/__init__.py

# Import standard packages
import os
from importlib import import_module

# Import included AddTag-specific modules
from .oligo import Oligo

# Import all Algorithm subclasses defined in python files within this same folder
exclusions = [
    '__init__.py',
    #'_primer3.py',
    '_primer3_draft1',
    #'_unafold.py',
    'addtagprimer.py',
    'bioprimer3.py',
    'oligos_draft1.py',
    'primer.py',
    'unafold_draft1.py',
]

path = os.path.dirname(os.path.abspath(__file__))
files = [f.rstrip(".py") for f in os.listdir(path) if (f.endswith('.py') and (f not in exclusions))]

for f in files:
    module = import_module('.'.join([__name__, f]))

# Create an instance of each Aligner subclass, and add to this list
oligos = []
for C in Oligo.__subclasses__():
    oligos.append(C())
