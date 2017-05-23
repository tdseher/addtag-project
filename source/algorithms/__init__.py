#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/__init__.py

# Import standard packages
import os
from importlib import import_module

# Import included AddTag-specific modules
from .algorithm import Algorithm, SingleSequenceAlgorithm, PairedSequenceAlgorithm, BatchedSingleSequenceAlgorithm

# Import all Algorithm subclasses defined in python files within this same folder
exclusions = [
    '__init__.py',
    'nucleotides.py',
    'azimuth_wrapper.py',
    'bae.py',
    'chari.py',
    'oof.py',
    'proxgc.py',
    'wang.py',
    'xu.py',
]
path = os.path.dirname(os.path.abspath(__file__))
files = [f.rstrip(".py") for f in os.listdir(path) if (f.endswith('.py') and (f not in exclusions))]
for f in files:
    module = import_module('.'.join([__name__, f]))

# Create an instance of each SingleSequenceAlgorithm subclass, and add to this list
single_algorithms = []
for C in SingleSequenceAlgorithm.__subclasses__():
    single_algorithms.append(C())

# Create an instance of each PairedSequenceAlgorithm subclass, and add to this list
paired_algorithms = []
#paired_algorithms_dict = {}
for C in PairedSequenceAlgorithm.__subclasses__():
    paired_algorithms.append(C())
    #obj = C()
    #paired_algorithms.append(obj)
    #paired_algorithms_dict[obj.name] = obj

# Create an instance of each SingleSequenceAlgorithm subclass, and add to this list
batched_single_algorithms = []
for C in BatchedSingleSequenceAlgorithm.__subclasses__():
    batched_single_algorithms.append(C())
