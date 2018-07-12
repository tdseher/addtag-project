#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/oligos/primer.py

# Import standard packages
#import sys
#import os
#import random

# import non-standard package
#import regex

class Primer(object):
    def __init__(self, sequence, position, strand, o_hairpin, o_self_dimer, o_reverse_complement, gc, checks=None):
        self.sequence = sequence
        self.position = position
        self.strand = strand
        self.o_hairpin = o_hairpin
        self.o_self_dimer = o_self_dimer
        self.o_reverse_complement = o_reverse_complement
        self.gc = gc
        if (checks == None):
            self.checks = []
        else:
            self.checks = checks
    
    def get_tm(self):
        if (self.o_reverse_complement.__class__.__name__ == 'ThermoResult'):
            return self.o_reverse_complement.tm
        else:
            return min(self.o_reverse_complement).melting_temperature
    
    def get_min_delta_G(self):
        if (self.o_hairpin.__class__.__name__ == 'ThermoResult'):
            return min(self.o_hairpin.dg, self.o_self_dimer.dg)/1000
        else:
            return min(self.o_hairpin + self.o_self_dimer).delta_G
    
    def __repr__(self):
        labs = ['seq', 'pos', 'strand', 'Tm', 'GC', 'min(dG)']
        vals = [
            self.sequence,
            self.position,
            self.strand,
            round(self.get_tm(), 2),
            round(self.gc, 2),
            round(self.get_min_delta_G(), 2)
        ]
        return self.__class__.__name__ + '(' + ', '.join('='.join(map(str, x)) for x in zip(labs, vals)) + ')'

class PrimerPair(object):
    def __init__(self, forward_primer_object, reverse_primer_object, o_heterodimer=None, checks=None):
        self.forward_primer = forward_primer_object
        self.reverse_primer = reverse_primer_object
        self.o_heterodimer = o_heterodimer
        if (checks == None):
            self.checks = []
        else:
            self.checks = checks
    
    def get_min_delta_G(self):
        if self.o_heterodimer:
            if (self.o_heterodimer.__class__.__name__ == 'ThermoResult'):
                return min(self.forward_primer.get_min_delta_G(), self.reverse_primer.get_min_delta_G(), self.o_heterodimer.dg/1000)
            else:
                return min(
                    self.forward_primer.o_hairpin +
                    self.forward_primer.o_self_dimer +
                    self.reverse_primer.o_hairpin +
                    self.reverse_primer.o_self_dimer +
                    self.o_heterodimer).delta_G
        else:
            return min(self.forward_primer.get_min_delta_G(), self.reverse_primer.get_min_delta_G())
    
    def get_amplicon_size(self):
        start = self.forward_primer.position
        end = self.reverse_primer.position+len(self.reverse_primer.sequence)
        return end-start
    
    def get_tms(self):
        return (self.forward_primer.get_tm(), self.reverse_primer.get_tm())
    
    def get_gcs(self):
        return (round(self.forward_primer.gc, 2), round(self.reverse_primer.gc, 2))
    
    def __repr__(self):
        labs = ['seq', 'amplicon_size', 'Tm', 'GC', 'min(dG)']
        vals = [
            (self.forward_primer.sequence, self.reverse_primer.sequence),
            self.get_amplicon_size(),
            tuple(round(x, 2) for x in self.get_tms()),
            self.get_gcs(),
            round(self.get_min_delta_G(), 2)
        ]
            
        return self.__class__.__name__ + '(' + ', '.join('='.join(map(str, x)) for x in zip(labs, vals)) + ')'