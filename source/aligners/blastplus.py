#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/aligners/blastplus.py

# List general Python imports
import sys
import os
import subprocess

# import non-standard package
import regex

# import AddTag-specific packages
from .. import utils

# Requires BLAST+ >= 2.6.0
#  because it has SAM format output
#  $ blastn -query query.fasta -subject subject.fasta -outfmt 17 > output.sam

def align(outfile, queryfile, subjectfile):
    """Aligns sequences using BLAST+"""
    return None

def test():
    """Code to test the classes and functions in 'source/blat.py'"""
    
    print("=== align ===")
    queryfile = 'alignment.fasta'
    subjectfile = 'subject.fasta'
    outfile = 'alignment.blastn'
    pslx = align(outfile, queryfile, subjectfile)

if (__name__ == '__main__'):
    test()