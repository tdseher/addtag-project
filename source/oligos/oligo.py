#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/oligos/oligo.py

# Import standard packages
import sys
import os
import inspect

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
    """
    Class that stores 2 Primer objects and their heterodimer thermodynamic calculations.
    """
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
    
    def get_formatted(self):
        interspace = self.reverse_primer.position - self.forward_primer.position - len(self.forward_primer.sequence)
        return self.forward_primer.sequence + '-'*interspace + rc(self.reverse_primer.sequence)
        #print(' '*self.forward_primer.position + fp.sequence)
        #print(' '*self.reverse_primer.position + rc(rp.sequence))
    
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

def lr_justify(left_text, right_text, width=140):
    if (len(left_text)+len(right_text) < width):
        text = ('{:>'+str(width)+'}').format(right_text)
        return left_text + text[len(left_text):]
    else:
        half = width//2 - 1
        if (len(right_text) > half):
            return left_text[:half] + '.'*(width-half-half) + right_text[-half:]
        else:
            return left_text[:width-len(right_text)-3] + '.'*3 + right_text
