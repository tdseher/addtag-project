#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/aligners/blat.py

# List general Python imports
import sys
import os
import subprocess

# import non-standard package
import regex

# import AddTag-specific packages
from .. import utils

def align(outfile, queryfile, subjectfile):
    """Aligns sequences using BLAT"""
    # blat subject.fasta query.fasta -t=dna -q=dna -tileSize=8 -stepSize=4 -oneOff=1 -minMatch=1 -minIdentity=90 -out=pslx alignment.psl
    return None

def test():
    """Code to test the classes and functions in 'source/blat.py'"""
    
    print("=== align ===")
    queryfile = 'alignment.fasta'
    subjectfile = 'subject.fasta'
    outfile = 'alignment.pslx'
    pslx = align(outfile, queryfile, subjectfile)

if (__name__ == '__main__'):
    test()