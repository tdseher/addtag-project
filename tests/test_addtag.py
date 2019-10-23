#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# tests/test_addtag.py

import sys
import subprocess

import pytest

def test_no_args():
    # Linux: exit status = 1
    # Windows: exit status = 0
    cp = subprocess.run([sys.executable, 'addtag'], shell=False, check=False, stdout=subprocess.PIPE)
    output = cp.stdout.decode().rstrip()
    assert output == 'usage: addtag [-h] [-v] action ...'

def test_helps():
    from source import subroutines
    for sub in [x.name for x in subroutines.subroutines]:
        cp = subprocess.run([sys.executable, 'addtag', sub, '-h'], shell=False, check=True, stdout=subprocess.PIPE)
        #output = cp.stdout.decode().rstrip()
        assert cp.returncode == 0

