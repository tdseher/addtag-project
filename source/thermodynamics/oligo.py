#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/thermodynamics/oligo.py

# Import standard packages
import sys
import os
import math
import inspect

# Import non-standard packages
import regex

# Treat modules in PACKAGE_PARENT as in working directory
if ((__name__ == "__main__") or (os.path.basename(inspect.stack()[-1][1]) in ['_primer3.py', '_unafold.py'])): # or __name__ == 'unafold'):
    # Relative path for package to import
    PACKAGE_PARENT = '..'
    # Obtain path of currently-running file
    SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
    # Convert to absolute path, and add to the PYTHONPATH
    sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))
    
    from nucleotides import rc
else:
    from ..nucleotides import rc

def logistic_up(x, upslope=8, up=0, height=1.0):
    return height/(1+upslope**(-x+up))

def logistic_down(x, downslope=8, down=0, height=1.0):
    return height/(1+downslope**(x-down))

def logistic_updown(x, upslope, up, downslope, down, height=1.0):
    return height * logistic_up(x, upslope, up, 1.0) * logistic_down(x, downslope, down, 1.0)

def normal_pdf(x, mean=0, std=1):
    return math.exp((x-mean)**2/(-2*std**2))/(2*math.pi*std**2)**0.5

def gamma(z):
    return math.factorial(z-1)

def gamma_pdf(x, shape, scale=1):
    return (x**(shape-1) * math.exp(-x/scale))/(scale**shape * gamma(shape))

class Oligo(object): # Name of the subclass
    """
    General class that should be subclassed when adding a new thermodynamics
    calculation program.
    """
    def __init__(self, 
        name,
        author,
        year,
        citation=None,
    ):
        """
        Specify general information regarding this new instance of the 
        Oligo class.
        """
        self.name = name             # Unique name for the Oligo subprocess (str). No other Oligo objects should have this name.
        self.author = author         # Author of the algorithm (str)
        self.year = year             # Year algorithm published (int)
        self.citation = citation     # Citation (None/str)
    
    def scan_sequence(self, *args, **kwargs):
        """
        Defines interface for calculation.
        
        Overload this method.
        """
        return None
    
    def scan(self, seq, side, *args, **kwargs):
        """
        Pass in a sequence, then either design left/forward primers
        or design right/reverse primers.
        Returns list of decent primers.
        
        Overload this method.
        """
        if (side in ['left', 'forward']):
            pass
        elif (side in ['right', 'reverse']):
            pass
        
        return None
    
    def pair(self, forward_list, reverse_list, *args, **kwargs):
        """
        Builds valid sets of left/forward and right/reverse primers.
        Typically, it checks that the Tm of each pair are within a few degrees,
        the amplicon length is within desired parameters, and there are no
        low energy secondary structures.
        
        Overload this method.
        """
        return None
    
    def group_weight(self, primers, *args, **kwargs):
        w = 1.0
        
        # Weigh by Tm
        temps = [p.get_tm() for p in primers if p] # will not include any p == None
        mean_temp = sum(temps)/len(temps) # mean Tm
        #print('temps =', temps, file=sys.stderr)
        for i, tm in enumerate(temps):
            w *= logistic_updown(tm-mean_temp, 5, -2.5, 5, 2.5)
        
        return w
    
    def get_3prime_complementation_length(self, seq1, seq2, max_3prime_search_length=5):
        """
        Calculates the length of 3' complementation for both sequences up to
        max_3prime_length.
        
        Counts length of complementation where the 3' end of one sequence could
        bind to another sequence.
        
        Finds sequences like this, where the terminal 3' end of both sequences
        are complements of each other:
          seq1 5'-ACAATACGAC-3'
                        ||||
          seq2       3'-GCTGTTAAG-5'
        
        And also sequences like this, where the terminal 3' end of one sequence
        is complemented on the inside of the other sequence:
          seq1 5'-GCTCTAAGATCACA-3'
                            ||||
          seq2      3'-CCTGGGTGTGAACT-5'
        """
        # Reverse--NOT reverse complement
        #rev2 = seq2[::-1]
        rev2 = rc(seq2)
        match_length = 0
        
        for i in range(max_3prime_search_length, 0, -1):
            if regex.search(seq1[-i:], rev2):
                match_length = i
                break
            if regex.search(rev2[:i], seq1):
                match_length = i
                break
        
        return match_length
    
    def simulate_amplification(self, primer_pair, contigs, amplicon_range=(10,5000)):
        """
        Simulate a PCR experiment, and return all the predicted amplicons.
        Does not predict chimeric (template-jumping) amplification.
        """
        regex_1F = regex.compile('('+primer_pair.forward_primer.sequence+'){s<=1}', flags=regex.IGNORECASE)
        regex_1R = regex.compile('('+rc(primer_pair.reverse_primer.sequence)+'){s<=1}', flags=regex.IGNORECASE)
        regex_2F = regex.compile('('+primer_pair.reverse_primer.sequence+'){s<=1}', flags=regex.IGNORECASE)
        regex_2R = regex.compile('('+rc(primer_pair.forward_primer.sequence)+'){s<=1}', flags=regex.IGNORECASE)
        
        amplicons = []
        for name, sequence in contigs.items():
            starts, stops = [], []
            for m in regex_1F.finditer(sequence):
                starts.append(m)
            for m in regex_1R.finditer(sequence):
                stops.append(m)
            for m in regex_2F.finditer(sequence):
                starts.append(m)
            for m in regex_2R.finditer(sequence):
                stops.append(m)
            
            for start_m in starts:
                for stop_m in stops:
                    start = start_m.start()
                    stop = stop_m.end()
                    size = stop - start
                    if (amplicon_range[0] <= size <= amplicon_range[1]):
                        amplicons.append((name, start, stop, size)) # 0-indexed
        return amplicons
    
    def __repr__(self):
        """
        Return the string representation of the Oligo
        """
        return self.__class__.__name__ + '(' + ', '.join(['name='+repr(self.name), 'author='+repr(self.author), 'year='+repr(self.year)]) + ')'


class Primer(object):
    """
    Class to store general information on an oligonucleotide sequence.
    Specifically, it stores the program-specific conformations and deltaGs
    as attributes.
    """
    def __init__(self, sequence, position, template_length, strand, o_hairpin, o_self_dimer, o_reverse_complement, gc=None, checks=None, name=None):
        self.sequence = sequence
        self.position = position
        self.template_length = template_length
        self.strand = strand
        self.o_hairpin = o_hairpin
        self.o_self_dimer = o_self_dimer
        self.o_reverse_complement = o_reverse_complement
        if (gc == None):
            self.gc = self.get_gc()
        else:
            self.gc = gc
        if (checks == None):
            self.checks = []
        else:
            self.checks = checks
        self.weight = self.get_weight()
        self.name = name
    
    def get_gc(self):
        if (self.gc != None):
            return self.gc
        else:
            C_count = seq.count('C')
            G_count = seq.count('G')
            return (C_count+G_count)/len(seq)
    
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
    
    def get_weight(self, *args, **kwargs):
        """
        Use logistic multiplier to weigh each Primer component
        """
        # Begin weight calculation with maximum score
        w = 1.0
        
        # Weigh by sequence length (logarithmic, with optimal at 25)
        #w *= gamma_pdf(len(self.sequence), shape=25, scale=1)
        w *= logistic_updown(len(self.sequence), upslope=3, up=17, downslope=1.7, down=30)
        
        # weigh by %GC
        #w *= normal_pdf(self.gc, mean=0.5, std=0.07)
        w *= logistic_updown(self.gc*100, upslope=1.7, up=40, downslope=1.7, down=60)
        
        # Weigh by minimum delta-G
        w *= logistic_up(self.get_min_delta_G(), upslope=3, up=-4)
        
        # Weigh by Tm
        # skip
        
        return w
    
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
    """
    Class that stores 2 Primer objects and their heterodimer thermodynamic calculations.
    """
    def __init__(self, forward_primer_object, reverse_primer_object, o_heterodimer=None, checks=None, *args, intervening=0, **kwargs):
        self.forward_primer = forward_primer_object
        self.reverse_primer = reverse_primer_object
        self.o_heterodimer = o_heterodimer
        self.intervening = intervening
        if (checks == None):
            self.checks = []
        else:
            self.checks = checks
        self.weight = self.get_weight(**kwargs)
    
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
        #                      missing vvvv
        #  forward     ............===>....
        #  intervening                     .........
        #  reverse                                  ........<===...
        
        # Assuming forward_primer and reverse_primer are separated by self.intervening (i.e. have different templates)
        # If forward_primer and reverse_primer are on the same template, then this works if template_length=0
        return (self.forward_primer.template_length - self.forward_primer.position) + self.intervening + (self.reverse_primer.position + len(self.reverse_primer.sequence))
    
    def get_tms(self):
        return (self.forward_primer.get_tm(), self.reverse_primer.get_tm())
    
    def get_gcs(self):
        return (round(self.forward_primer.gc, 2), round(self.reverse_primer.gc, 2))
    
    def get_formatted(self):
        ####### This needs to be re-done (it isn't correct) ####### <---- I made changes, so it needs to be re-checked
        #interspace = self.reverse_primer.position - (self.intervening + self.forward_primer.position + len(self.forward_primer.sequence))
        interspace = self.forward_primer.template_length - self.forward_primer.position - len(self.forward_primer.sequence) + self.intervening + self.reverse_primer.position
        return self.forward_primer.sequence + '·'*interspace + rc(self.reverse_primer.sequence)
        #print(' '*self.forward_primer.position + fp.sequence)
        #print(' '*self.reverse_primer.position + rc(rp.sequence))
    
    def get_weight(self, *args, **kwargs):
        """
        Use logistic multiplier to weigh the PrimerPair
        """
        w = 1.0
        
        # Weigh by amplicon size
        # If this is supposed to be a 'sF' 'sR' pair, then we want to minimize the amplicon size
        # Otherwise, we want the amplicon size to be around 500 nt
        if (("minimize" in kwargs) and kwargs["minimize"]):
            w *= logistic_down(self.get_amplicon_size(), downslope=1.03, down=150)
        else:
            w *= logistic_updown(self.get_amplicon_size(), 1.03, 300, 1.03, 700)
        
        # Weigh by minimum delta-G
        w *= logistic_up(self.get_min_delta_G(), upslope=3, up=-4)
        
        # Weigh by Tm difference
        def diff(x,y):
            return x-y
        w *= logistic_updown(diff(*self.get_tms()), 5, -2.5, 5, 2.5)
        
        return w
    
    def get_joint_weight(self):
        # If the sequences are identical, then use only one weight
        if (self.forward_primer.sequence == self.reverse_primer.sequence):
            return self.weight * ((self.forward_primer.weight + self.reverse_primer.weight)/2)
        else:
            return self.weight * self.forward_primer.weight * self.reverse_primer.weight
    
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

