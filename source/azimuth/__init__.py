#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/azimuth/__init__.py

# List general Python imports
import sys
import os
import subprocess

import regex

# Azimuth requires Python 2.7
# as well as specific versions of pandas, scipy, scikit-learn, etc
# So rather than call it natively, we will open it using subprocess

def batch_azimuth_score(sequence_list):
    queries = []
    
    for query in sequence_list:
        seq, pam, upstream, downstream = query
        queries.append(build_query(seq, pam, upstream, downstream))
    
    # Obtain path of currently-running file
    SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
    WRAPPER_PATH = os.path.join(SCRIPT_DIR, "wrapper.py")
    # Should use args.python2_path
    command_list = ['python', WRAPPER_PATH] + queries
    
    #with open(error_file, 'w+') as flo:
    #    cp = subprocess.run(command_list, shell=False, check=True, stdout=flo, stderr=subprocess.STDOUT)
    cp = subprocess.run(command_list, shell=False, check=True, stdout=subprocess.PIPE)
    output = cp.stdout.decode()
    
    o_scores = []
    for line in output.splitlines():
        if not line.startswith('No model file specified'):
            i_seq, i_score = line.split(" ")
            o_scores.append(100*float(i_score))
    
    return o_scores
    
    
    
def build_query(seq, pam, upstream='', downstream=''):
    # us     seq                pam  ds
    # ACAG CTGATCTCCAGATATGACCA|TGG GTT
    # CAGC TGATCTCCAGATATGACCAT|GGG TTT
    # CCAG AAGTTTGAGCCACAAACCCA|TGG TCA
    
    new = ['-'] * 30
    for i, nt in enumerate(pam):
        new[24+i] = nt

    for i, nt in enumerate(seq[::-1]):
        new[24-1-i] = nt

    for i, nt in enumerate(upstream[::-1]):
        ind = 24-len(seq)-1-i
        if (ind < 0):
            break
        new[ind] = nt

    for i, nt in enumerate(downstream):
        ind = 24+len(pam)+i
        if (ind >= len(new)):
            break
        new[ind] = nt
    query = ''.join(new)
    return query

def azimuth_score(seq, pam, upstream='', downstream=''):
    #error_file = os.path.join(folder, name + '.err')
    #error_file = os.path.join('azimuth.err')
    
    query = build_query(seq, pam, upstream, downstream)
    
    if (len(query) != 30):
        return 0.0
    
    m = regex.search('[^ATCGatcg]', query)
    if m:
        return 0.0
    
    # Obtain path of currently-running file
    SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
    WRAPPER_PATH = os.path.join(SCRIPT_DIR, "wrapper.py")
    # Should use args.python2_path
    command_list = ['python', WRAPPER_PATH, query]
    
    #with open(error_file, 'w+') as flo:
    #    cp = subprocess.run(command_list, shell=False, check=True, stdout=flo, stderr=subprocess.STDOUT)
    cp = subprocess.run(command_list, shell=False, check=True, stdout=subprocess.PIPE)
    output = cp.stdout.decode()
    
    o_seq = None
    o_score = 0.0
    for line in output.splitlines():
        if not line.startswith('No model file specified'):
            o_seq, o_score = line.split(" ")
    
    return 100*float(o_score)

def test():
    """Code to test the classes and functions in 'source/azimuth/__init__.py'"""
    
    us, seq, pam, ds = 'ACAG', 'CTGATCTCCAGATATGACCA', 'TGG', 'GTT'
    print(us+' '+seq+'|'+pam+' '+ds, azimuth_score(seq, pam, upstream=us, downstream=ds))

if (__name__ == "__main__"):
    test()
