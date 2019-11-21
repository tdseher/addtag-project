#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# tests/test_thermodynamics.py

from source import thermodynamics
from source import nucleotides

def around(x, center, spread):
    return center-spread <= x <= center+spread

def test_oligo():
    for oligo in thermodynamics.oligos:
        if oligo.available:
            print("=== ", oligo.name, " ===")

            a = 'GAAATCGCTTAGCGCGAACTCAGACCAT'
            b = 'CCTAGCTATTTAATAAATC'
            c = 'TTCTCCACTTCCATCACCGT'

            tempdir = '.'

            print('Hairpin: {}'.format(repr(a)))
            hairpins = oligo.find_structures(tempdir, a, None)
            for s in hairpins:
                print('', s)

            print('Homodimer: {} {}'.format(repr(a), repr(a)))
            homodimers = oligo.find_structures(tempdir, a, a)
            for s in homodimers:
                print('', s)

            print('Heterodimer: {} {}'.format(repr(a), repr(b)))
            heterodimers = oligo.find_structures(tempdir, a, b)
            for s in heterodimers:
                print('', s)

            print('Reverse-complements: {} {}'.format(repr(a), repr(nucleotides.rc(a))))
            revcomps = oligo.find_structures(tempdir, a, nucleotides.rc(a))
            for s in revcomps:
                if (s.melting_temperature == None):
                    tm_list = oligo.find_tms([s.seq1])
                    s.melting_temperature = tm_list[0]
                print('', s)

            if (oligo.name == 'Primer3'):
                assert around(hairpins[0].delta_G, -3.8094209545944215, 0.1)
                assert around(homodimers[0].delta_G, -9.363839560726309, 0.1)
                assert around(heterodimers[0].delta_G, -4.013193909188849, 0.1)
                assert around(revcomps[0].delta_G, -35.37314759135284, 0.1)
                assert around(revcomps[0].melting_temperature, 62.994719476938315, 1)
            elif (oligo.name == 'UNAFold'):
                assert around(hairpins[0].delta_G, -3.43067, 0.1)
                assert around(homodimers[0].delta_G, -8.09575, 0.1)
                assert around(heterodimers[0].delta_G, -3.54803, 0.1)
                assert around(revcomps[0].delta_G, -35.2562, 0.1)
                assert around(revcomps[0].melting_temperature, 62.798723427739105, 1)
            elif (oligo.name == 'ViennaRNA'):
                assert around(hairpins[0].delta_G, -3.5, 0.1)
                assert around(homodimers[0].delta_G, -10.399999618530273, 0.1)
                assert around(heterodimers[0].delta_G, -4.800000190734863, 0.1)
                assert around(revcomps[0].delta_G, -37.900001525878906, 0.1)
                assert around(revcomps[0].melting_temperature, 63.51, 1)
            else:
                assert False, 'No thermodyanmics calculator detected'
        else:
            print("=== ", oligo.name, " ===")
            print("    UNAVAILABLE")
