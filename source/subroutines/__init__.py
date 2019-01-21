#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/subroutines/__init__.py

# Import standard packages
import os
from importlib import import_module

# Import included AddTag-specific modules
from .subroutine import Subroutine, CustomHelpFormatter, ValidateFlanktags, ValidateKodDNA

def prefix_strip(text, prefix):
    return text[len(prefix):] if text.startswith(prefix) else text

def suffix_strip(text, suffix):
    return text[:-len(suffix)] if text.endswith(suffix) else text

# Import all Algorithm subclasses defined in python files within this same folder
exclusions = [
    '__init__.py',
    'subroutine.py'
]
path = os.path.dirname(os.path.abspath(__file__))
#files = [f.rstrip(".py") for f in os.listdir(path) if (f.endswith('.py') and (f not in exclusions))]
files = [suffix_strip(f, '.py') for f in os.listdir(path) if (f.endswith('.py') and (f not in exclusions))]

for f in files:
    module = import_module('.'.join([__name__, f]))

# Create an instance of each Subroutine subclass, and add to this list
subroutines = []

def make_subroutines(subparsers):
    # Clear contents of 'subroutines' list, but keep the memory address of the object the same
    subroutines.clear() # Same as 'del subroutines[:]'
    
    # Re-populate the 'subroutines' list with all 'Subroutine' objects
    for S in Subroutine.__subclasses__():
        subroutines.append(S(subparsers))
