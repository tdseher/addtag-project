#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# tests/test_nucleotides.py

from source import nucleotides

def test_rc():
    seq = 'AaCcGgTtRrYyMmKkWwSsBbDdHhVvNn'
    rcseq = 'nNbBdDhHvVsSwWmMkKrRyYaAcCgGtT'
    assert nucleotides.rc(seq) == rcseq

