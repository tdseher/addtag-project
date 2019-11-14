#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# tests/test_cigarstrings.py

# List general Python imports
import os
import math

# Import AddTag-specific packages
from source.cigarstrings import *
from source import evalues

def test_specificate_cigar():
    q = 'AAAAAATCCC'
    s = 'AGACGAAACCC'
    c = '3M2D3MI3M'
    assert specificate_cigar(c, q, s) == '=X=DD===I==='

def test_collapse_cigar():
    c = '=X=DD===I==='
    assert collapse_cigar(c, abbreviate=False) == '1=1X1=2D3=1I3='
    assert collapse_cigar(c, abbreviate=True) == '=X=2D3=I3='

def test_expand_cigar():
    assert expand_cigar('1=1X1=2D3=1I3=') == '=X=DD===I==='
    assert expand_cigar('=X=2D3=I3=') == '=X=DD===I==='

def test_cigar_errors():
    assert cigar_errors('=X=2D3=I3=') == (8, 1, 1, 2, 4)

def test_matrix_averages():
    data_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
    score_path = os.path.join(data_folder, 'source', 'algorithms', 'distance_scores.txt')
    score_matrix = evalues.load_scores(score_path)
    a = matrix_averages(score_matrix)
    
    assert a[0] == 1.0
    assert math.isclose(a[1], -1.133333333333333)
    assert a[2] == -3.0
    assert a[3] == -3.0
    assert a[4] == -1.8

def test_cigar2score():
    data_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
    
    score_path = os.path.join(data_folder, 'source', 'algorithms', 'distance_scores.txt')
    score_matrix = evalues.load_scores(score_path)
    print(cigar2score('1=2X2=1I1=1X10=1X4=', score_matrix))
    print(cigar2score('23=', score_matrix))
    
    score_path = os.path.join(data_folder, 'source', 'aligners', 'bowtie2_scores.txt')
    score_matrix = evalues.load_scores(score_path)
    print(cigar2score('1=2X2=1I1=1X10=1X4=', score_matrix))
    print(cigar2score('23=', score_matrix))

def test_cigar2query_position():
    pass

def test_cigar2query_aligned_length():
    pass

def test_cigar2subject_aligned_length():
    pass

def test_reverse_cigar():
    pass

def test_match2cigar():
    pass

def test_alignment2cigar():
    q = "AAA--AAATCCC"
    s = "AGACGAAA-CCC"
    assert alignment2cigar(q, s) == '3M2D3M1I3M'
    assert alignment2cigar(q, s, abbreviated=True) == '3M2D3MI3M'
    assert alignment2cigar(q, s, specific=True) == '1=1X1=2D3=1I3='
    assert alignment2cigar(q, s, specific=True, abbreviated=True) == '=X=2D3=I3='

def test_sam_orientation():
    pass

def test_decode_sam_flags():
    pass
