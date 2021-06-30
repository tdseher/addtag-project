#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/aligners/__init__.py

# Import standard packages
import os
from importlib import import_module

# Import included AddTag-specific modules
from .aligner import PairwiseAligner, MultipleSequenceAligner

def suffix_strip(text, suffix):
    return text[:-len(suffix)] if text.endswith(suffix) else text

# Import all Algorithm subclasses defined in python files within this same folder
exclusions = [
    '__init__.py',
    'aligner.py'
#    'blastplus.py',
    'blat.py',
    'bowtie.py',
#    'bowtie2.py',
#    'bwa.py',
#    'casoffinder.py',
    'casot.py',
    'dsnickfury.py',
    'hyperscan.py',
    'inbuilt.py',
    'kalign3.py',
#    'mafft.py',
    'python-bwa.py',
    'python-nmslib.py',
    'seqmap.py',
    'star.py',
    'usearch.py',
    'varscot.py',
]

path = os.path.dirname(os.path.abspath(__file__))
#files = [f.rstrip(".py") for f in os.listdir(path) if (f.endswith('.py') and (f not in exclusions))]
files = [suffix_strip(f, '.py') for f in os.listdir(path) if (f.endswith('.py') and (f not in exclusions))]

for f in files:
    module = import_module('.'.join([__name__, f]))

# Create an instance of each Aligner subclass, and add to this list
pw_aligners = []
for C in PairwiseAligner.__subclasses__():
    pw_aligners.append(C())

ms_aligners = []
for C in MultipleSequenceAligner.__subclasses__():
    ms_aligners.append(C())
