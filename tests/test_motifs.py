#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# tests/test_motifs.py

from source.motifs import *

def test_motif_class():
    m = Motif(r'N{17,20}|N{3,4}|AAA>NR{0,1}GG')
    print('motif_string = {}'.format(m.motif_string))
    print('parsed_list = {}'.format(m.parsed_list))
    print('spacer_sense_cuts = {}'.format(m.spacer_sense_cuts))
    print('spacer_antisense_cuts = {}'.format(m.spacer_antisense_cuts))
    print('pam_sense_cuts = {}'.format(m.pam_sense_cuts))
    print('pam_antisense_cuts = {}'.format(m.pam_antisense_cuts))
    print('regex_string = {}'.format(m.regex_string))

