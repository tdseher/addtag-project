#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/azimuth/__init__.py

# List general Python imports
import sys
import os
import subprocess

# Azimuth requires Python 2.7
# as well as specific versions of pandas, scipy, scikit-learn, etc
# So rather than call it natively, we will open it using subprocess

def azimuth_score(seq, pam, upstream='', downstream=''):
    #error_file = os.path.join(folder, name + '.err')
    error_file = os.path.join('azimuth.err')
    
    command_list = ['python', 'azimuth.py', seq]
    
    with open(error_file, 'w+') as flo:
        cp = subprocess.run(command_list, shell=False, check=True, stdout=flo, stderr=subprocess.STDOUT)
    


def test():
    """Code to test the classes and functions in 'source/azimuth/__init__.py'"""
    
    pass

if (__name__ == "__main__"):
    test()
