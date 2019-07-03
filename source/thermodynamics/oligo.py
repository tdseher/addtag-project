#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/thermodynamics/oligo.py

# Import standard packages
import sys
import os
import math
import inspect
import time
import logging
import copy
from collections import defaultdict
import random

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
    
    from nucleotides import rc, get_gc_freq
    import utils
else:
    from ..nucleotides import rc, get_gc_freq
    from .. import utils

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
    
    # def scan_sequence(self, *args, **kwargs):
    #     """
    #     Defines interface for calculation.
    #     
    #     Overload this method.
    #     """
    #     return None
    
    # def scan(self, seq, side, *args, **kwargs):
    #     """
    #     Pass in a sequence, then either design left/forward primers
    #     or design right/reverse primers.
    #     Returns list of decent primers.
    #     
    #     Overload this method.
    #     """
    #     if (side in ['left', 'forward', '+']):
    #         pass
    #     elif (side in ['right', 'reverse', '-']):
    #         pass
    #     
    #     return None
    
    def scan(self, sequence, locations=None):
        '''
        Use this function to scan substrings that have (shift) coordinates
          scan(sequence, contig, start, end, strand='+', external=True)
            # scan = sequence
            # primer_pos = (s, e)
            # location = (contig, start+s, start+e, strand)
          
        OR
        
        Use this function to scan a substring BY its coordinates
          scan(sequence, contig, start, end, strand='+', internal=True)
            # scan = sequence[start:end]
            # primer_pos = (s, e)
            # location = (contig, start+s, start+e, strand)
        
        
        '''
        pass
    
    def scan(self, seq, side, *args, primer_size=(18,26), us_seq='', ds_seq='', min_junction_overlap=(4,8), time_limit=60, **kwargs):
        """
        Pass in a sequence, then either design left/forward/+ primers
        or design right/reverse/- primers.
        Returns list of decent primers.
        """
        # 'min_junction_overlap' allows for primers to span into the us/ds sequences
        #           uuuuuIIIIddddd
        #              └primer┘
        # 
        # min_junction_overlap = (
        #    minimum nt of primer that CANNOT be on the 5' flanking sequence (us_seq for forward primer),
        #    minimum nt of primer that CANNOT be on the 3' flanking sequence (ds_seq for forward primer)
        # )
        # 
        # Example of min_junction_overlap=(1,3) (for forward primer)
        #     us_seq    GTAAAGAAGG             10 nt
        #     seq                 CC            2 nt
        #     ds_seq                CCAAATTA    8 nt
        #                         ┌1 nt overlap
        #     primers         FFFFF 
        #                      FFFFF
        #                       FFFFF
        #     4 primers          FFFFF
        #                        └─┴3 nt overlap
        # Example:
        #       xxxxxxxxxx len(seq) = 10
        #   FFFFF FFFFF    len(primer) = 5
        #    FFFFF FFFFF   min_junction_overlap = (1,3)
        #     FFFFF FFFFF  12 primers
        #      FFFFF FFFFF
        #       FFFFF FFFFF
        #        FFFFF FFFFF
        #    ffff ffff ffff len(primer) = 4
        #     ffff ffff     11 primers
        #      ffff ffff
        #       ffff ffff
        #        ffff ffff
        # Code to temporarily limit number of primers computed set to None to disable
        #subset_size = 5000 # 1000 # 9000
        
        time_expired = False
        start_time = time.time() # Time in seconds
        
        # Total number of primers that will be scanned
        est_a = lambda x: max(0, min(x-min_junction_overlap[0], len(us_seq))) # number primers overlapping with us_seq
        est_b = lambda x: max(0, min(x-min_junction_overlap[1], len(ds_seq))) # number primers overlapping with ds_seq
        est_c = lambda x: len(seq)-x+1 # number primers contained within seq
        
        n_est = sum(max(0, est_a(x)+est_b(x)+est_c(x)) for x in range(primer_size[0], primer_size[1]+1))
        n_max = 0
        n_count = 0
        
        logging.info("  Primer scan: Scanning '{}'".format(side))
        
        # First, create the primers
        # Second, filter the Primer objects
        # Third, return only the good Primer objects
        good_primers = []
        if (side in ['left', 'L', 'forward', 'F', '+']):
            for p_len in range(primer_size[0], primer_size[1]+1):
                # Include the flanking regions
                l_seq = us_seq[len(us_seq) - (p_len-min_junction_overlap[0]):]
                r_seq = ds_seq[:max(0, p_len-min_junction_overlap[1])]
                
                j_seq = l_seq + seq + r_seq
                
                if (len(j_seq) != len(seq)):
                    logging.info('l_seq = ' + l_seq)
                    logging.info('  seq = ' + ' '*len(l_seq) + seq)
                    logging.info('r_seq = ' + ' '*len(l_seq+seq) + r_seq)
                
                for i in range(len(j_seq)-p_len+1):
                    if not time_expired:
                        pos = i-len(l_seq) # Position on sequence (can take negative number if hanging off the edge)
                        pf = j_seq[i:i+p_len]
                        
                        primer_object = Primer.create(sequence=pf, position=pos, template_length=len(seq), strand='+', o_oligo=self)
                        if primer_object.summarize(primer_object.checks):
                            good_primers.append(primer_object)
                        
                        n_count += 1
                    n_max += 1
                
                logging.info('  Primer scan: {}/~{}, time_expired: {}'.format(n_count, n_est, time_expired))
                #if ((subset_size != None) and (n_count >= subset_size)): ##### Temporary early break for testing purposes ##### Should remove eventually
                if ((time.time()-start_time) > time_limit):
                    time_expired = True
        
        elif (side in ['right', 'R', 'reverse', '-']):
            for p_len in range(primer_size[0], primer_size[1]+1):
                # Include the flanking regions
                l_seq = us_seq[len(us_seq) - (p_len-min_junction_overlap[1]):]
                r_seq = ds_seq[:max(0, p_len-min_junction_overlap[0])]
                
                j_seq = l_seq + seq + r_seq
                
                if (len(j_seq) != len(seq)):
                    logging.info('l_seq = ' + l_seq)
                    logging.info('  seq = ' + ' '*len(l_seq) + seq)
                    logging.info('r_seq = ' + ' '*len(l_seq+seq) + r_seq)
                
                for i in range(len(j_seq)-p_len+1):
                    if not time_expired:
                        pos = i-len(l_seq) # Position on sequence (can take negative number if hanging off the edge)
                        pr = rc(j_seq[i:i+p_len])
                    
                        primer_object = Primer.create(sequence=pr, position=pos, template_length=len(seq), strand='-', o_oligo=self)
                        if primer_object.summarize(primer_object.checks):
                            good_primers.append(primer_object)
                        
                        n_count += 1
                    n_max += 1
                
                logging.info('  Primer scan: {}/~{}, time_expired: {}'.format(n_max, n_est, time_expired))
                #if ((subset_size != None) and (n_count >= subset_size)): ##### Temporary early break for testing purposes ##### Should remove eventually
                if ((time.time()-start_time) > time_limit):
                    time_expired = True
        
        logging.info('  Primer scan: analyzed {}/{} potential primers'.format(n_count, n_max))
        logging.info('  Primer scan: found {} good primers'.format(len(good_primers)))
        logging.info('  Primer scan: time elapsed: {}s'.format(time.time()-start_time))
        
        return good_primers
    
    # def pair(self, forward_list, reverse_list, *args, **kwargs):
    #     """
    #     Builds valid sets of left/forward/+ and right/reverse/- primers.
    #     Typically, it checks that the Tm of each pair are within a few degrees,
    #     the amplicon length is within desired parameters, and there are no
    #     low energy secondary structures.
    #     
    #     Overload this method.
    #     """
    #     return None
    
    def pair(self, good_forwards, good_reverses, *args, intervening=0, same_template=False, time_limit=60, **kwargs):
        """
        Builds valid sets of left/forward/+ and right/reverse/- primers.
        Typically, it checks that the Tm of each pair are within a few degrees,
        the amplicon length is within desired parameters, and there are no
        low energy secondary structures.
        """
        # Time how long the process takes
        start_time = time.time()
        time_expired = False
        
        # Filter pairs by amplicon size
        n_count = 0
        n_max = 0
        n_est = len(good_forwards) * len(good_reverses)
        
        good_pairs = []
        queue = []
        for gf in good_forwards:
            for gr in good_reverses:
                queue.append((gf.get_weight()*gr.get_weight(), gf, gr))
        
        for weight, gf, gr in sorted(queue, reverse=True):
            if not time_expired:
                cgf, cgr = copy.deepcopy(gf), copy.deepcopy(gr)
                # Set their template lengths to 0 because they come from the same template
                if same_template:
                    cgf.template_length, cgr.template_length = 0, 0
                
                pp_object = PrimerPair.create(cgf, cgr, *args, intervening=intervening, o_oligo=self, **kwargs) # 'pp_object' may not have desired 'intervening' value
                if pp_object.forward_primer.summarize(pp_object.checks):
                    good_pairs.append(pp_object)
            
                n_count += 1
            n_max += 1
            if (n_max % 1000 == 0):
                logging.info('  Primer pair: {}/~{}, time_expired: {}'.format(n_max, n_est, time_expired))
                if ((time.time()-start_time) > time_limit):
                    time_expired = True
        
        logging.info('  Primer pair: analyzed {}/{} potential pairs'.format(n_count, n_max))
        logging.info('  Primer pair: found {} good primer pairs'.format(len(good_pairs)))
        logging.info('  Primer pair: time elapsed: {}s'.format(time.time()-start_time))
        
        return good_pairs # unsorted
    
    def find_structures(self, *args, **kwargs):
        """
        Should return the list of structures with delta-G values.
        
        Overload this method.
        """
        logging.info("THIS CODE SHOULD NOT BE EXECUTED")
        return []
    
    # Only in this class because of the 'logistic_updown()' function
    def p_group_weight(self, primers, merge_duplicates=True, *args, **kwargs):
        """
        Weigh by difference of each primer's Tm and the mean Tm
        """
        w = 1.0
        
        # First, we weed out any None, and take average Tm if a Primer is not unique
        if merge_duplicates:
            pdict = {}
            for p in primers:
                if p:
                    pdict.setdefault(p.sequence, []).append(p.get_tm())
            temps = []
            for pseq, tms in pdict.items():
                temps.append(sum(tms)/len(tms))
        else:
            temps = [p.get_tm() for p in primers if p] # will not include any p == None
        
        if (len(temps) > 0):
            mean_temp = sum(temps)/len(temps) # mean Tm
            #print('temps =', temps, file=sys.stderr)
            for i, tm in enumerate(temps):
                w *= logistic_updown(tm-mean_temp, 5, -2.5, 5, 2.5)
            
            return w
        else:
            return 0.0
    
    def pp_group_weight(self, primer_pairs, merge_duplicates=True, *args, **kwargs):
        """
        Weigh a list of PrimerPair objects.
        """
        weight = 1.0
        
        if merge_duplicates:
            ppdict = {}
            for pp in primer_pairs:
                if pp:
                    k = (pp.forward_primer.sequence, pp.reverse_primer.sequence)
                    ppdict.setdefault(k, []).append(pp.get_joint_weight())
            w_list = []
            for k, pp_w in ppdict.items():
                w_list.append(sum(pp_w)/len(pp_w))
            
        else:
            w_list = [x.get_joint_weight() for x in primer_pairs if x]
        
        for i, w in enumerate(w_list):
            weight *= w
        
        return weight
    
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
    
    sequences = {} # key: nucleotide sequence, value: 'Primer' object
    
    def __init__(self, sequence, gene, locus, genome, region, contig, strand, start, end, name=None):
        self.sequence = sequence
        
        # Build the location tuple
        location = (gene, locus, genome, region, contig, strand, start, end)
        
        # Make empty set to store locations where this primer originates from (NOT where it is aligned to)
        self.locations = set()
        
        # Keep track of names for this 'Primer'
        self.names = set()
        
        # Add this to the non-redundant dict of 'Primer' objects if it doesn't already exist.
        # Otherwise, add this locus to the pre-existing 'Primer' object
        the_primer = self.sequences.setdefault(self.sequence, self)
        the_primer.locations.add(location)
        
        if name:
            the_primer.names.add(name)
        
        #self.position = position
        #self.template_length = template_length
        #self.strand = strand
        
        # Set default values
        self.gc = math.nan # use NaN as undefined
        self.checks = [None]*9
        self.o_hairpin = None
        self.o_self_dimer = None
        self.o_reverse_complement = None
        self.weight = 0.0
        
        # Only do expensive calculations if this is a sequence that hasn't
        # been observed yet
        if (self == the_primer):
            self.gc = get_gc_freq(self.sequence)
            #self.check()
    
    def __lt__(self, other):
        '''
        Sort 'Primer' objects by their weights
        '''
        #if isinstance(other, Primer):
        return self.weight < other.weight
        #else:
        #    return False
        # __gt__(), __le__(), __ne__(), __ge__(), __eq__()
    
    def set_name(self, name):
        self.names.add(name)
    
    def get_name(self):
        return ','.join(sorted(self.names))
    
    def get_specificity(self):
        # location = (gene, locus, genome, region, contig, strand, start, end)
        return set(x[1] for x in self.locations)
    
    # def is_masked(self, args): # Code incomplete
    #     """
    #     Returns 'True' if the sequence should be masked
    #     (e.g. contains lower-case or ambiguous characters).
    #     """
    #     return False
    
    @classmethod
    def check_case(cls, case, seq):
        '''
        Returns 'True' if the sequence passes the case requirements,
        Otherwise, returns 'False'
        '''
        
        # Check the case of the potential gRNA sequence
        if (case == "upper-only"):
            if regex.search('[a-z]', seq):
                return False # Reject this sequence because it has lower-case characters
        elif (case == "lower-only"):
            if regex.search('[A-Z]', seq):
                return False # Reject this sequence because it has upper-case characters
        elif (case == "mixed-lower"):
            if not regex.search('[a-z]', seq):
                return False
        elif (case == "mixed-upper"):
            if not regex.search('[A-Z]', seq):
                return False
        elif (case == "mixed-only"):
            if not (regex.search('[a-z]', seq) and regex.search('[A-Z]', seq)):
                return False # Reject this sequence because it does not have both lower-case and upper-case characters
        #elif (args.case == "ignore") # then do nothing
        #    pass
        return True
    
    @classmethod
    def scan(cls, sequence, gene, locus, genome, region, contig, orientation, start, end, name, primer_size=(18,26), min_junction_overlap=(4,8), case='ignore', **kwargs):
        """
        Pass in a sequence, then either design left/forward/+ primers
        or design right/reverse/- primers.
        Returns list of decent primers.
        """
        # 'min_junction_overlap' allows for primers to span into the us/ds sequences
        #           uuuuuIIIIddddd
        #              └primer┘
        # 
        # min_junction_overlap = (
        #    minimum nt of primer that CANNOT be on the 5' flanking sequence (us_seq for forward primer),
        #    minimum nt of primer that CANNOT be on the 3' flanking sequence (ds_seq for forward primer)
        # )
        # min_junction_overlap = (
        #    number nt on 3' end of F primer that MUST be on seq (not us_seq),
        #    number nt on 5' end of F primer that MUST be on seq (not ds_seq)
        # )
        # 
        # Example of min_junction_overlap=(1,3) (for forward primer)
        #     us_seq    GTAAAGAAGG             10 nt
        #     seq                 CC            2 nt
        #     ds_seq                CCAAATTA    8 nt
        #                         ┌1 nt overlap
        #     primers         FFFFF 
        #                      FFFFF
        #                       FFFFF
        #     4 primers          FFFFF
        #                        └─┴3 nt overlap
        # Example:
        #       xxxxxxxxxx len(seq) = 10
        #   FFFFF FFFFF    len(primer) = 5
        #    FFFFF FFFFF   min_junction_overlap = (1,3)
        #     FFFFF FFFFF  12 primers
        #      FFFFF FFFFF
        #       FFFFF FFFFF
        #        FFFFF FFFFF
        #    ffff ffff ffff len(primer) = 4
        #     ffff ffff     11 primers
        #      ffff ffff
        #       ffff ffff
        #        ffff ffff
        
        
        seq = sequence[start:end]
        us_extend_start = max(0, start-(max(primer_size)-1))
        ds_extend_end = min(len(sequence), end+(max(primer_size)-1))
        us_seq = sequence[us_extend_start:start]
        ds_seq = sequence[end:ds_extend_end]
        
        # Total number of primers that will be scanned
        if ((min_junction_overlap == None) or (min_junction_overlap == (None, None))):
            est_a = lambda x: 0 # number primers overlapping with us_seq
            est_b = lambda x: 0 # number primers overlapping with ds_seq
            est_c = lambda x: len(seq)-x+1 # number primers contained within seq
        else:
            est_a = lambda x: max(0, min(x-min_junction_overlap[0], len(us_seq))) # number primers overlapping with us_seq
            est_b = lambda x: max(0, min(x-min_junction_overlap[1], len(ds_seq))) # number primers overlapping with ds_seq
            est_c = lambda x: len(seq)-x+1 # number primers contained within seq
        
        n_est = sum(max(0, est_a(x)+est_b(x)+est_c(x)) for x in range(primer_size[0], primer_size[1]+1))
        n_max = 0
        n_count = 0
        
        logging.info("  Primer scan: Scanning '{}'".format(orientation))
        
        # First, create the primers
        # Second, filter the Primer objects
        # Third, return only the good Primer objects
        
        mjo = min_junction_overlap
        
        if (orientation in ['left', 'L', 'l', 'forward', 'F', 'f', '+']):
            strand = '+'
            for p_len in range(primer_size[0], primer_size[1]+1):
                # Include the flanking regions
                
                if ((min_junction_overlap == None) or (min_junction_overlap == (None, None))):
                    mjo = (p_len, p_len)
                
                #l_seq = us_seq[len(us_seq) - (p_len-min_junction_overlap[0]):]
                if (p_len > mjo[0]):
                    l_start = -(p_len-mjo[0]) # Should be a negative number most of the time
                else:
                    l_start = len(us_seq)
                l_seq = us_seq[l_start:]
                
                r_end = max(0, p_len-mjo[1])
                r_seq = ds_seq[:r_end]
                
                j_seq = l_seq + seq + r_seq
                
                if (len(j_seq) != len(seq)):
                    logging.info('    l_seq = ' + l_seq)
                    logging.info('      seq = ' + ' '*len(l_seq) + seq)
                    logging.info('    r_seq = ' + ' '*len(l_seq+seq) + r_seq)
                
                for i in range(len(j_seq)-p_len+1):
                    #pos = i-len(l_seq) # Position on sequence (can take negative number if hanging off the edge)
                    pos = start-len(l_seq)+i
                    pf = j_seq[i:i+p_len]
                    
                    if cls.check_case(case, pf):
                        #primer_object = Primer.create(sequence=pf, position=pos, template_length=len(seq), strand='+', o_oligo=self)
                        #p = Primer(pf, gene, locus, genome, region, contig, strand, start+i, start+i+p_len, name=name)
                        p = Primer(pf, gene, locus, genome, region, contig, strand, pos, pos+p_len, name=name)
                    
                    n_count += 1
                    n_max += 1
                
                logging.info('    Primer scan: {}/~{}'.format(n_count, n_est))
        
        elif (orientation in ['right', 'R', 'r', 'reverse', '-']):
            strand = '-'
            for p_len in range(primer_size[0], primer_size[1]+1):
                # Include the flanking regions
                
                if ((min_junction_overlap == None) or (min_junction_overlap == (None, None))):
                    mjo = (p_len, p_len)
                
                #l_seq = us_seq[len(us_seq) - (p_len-min_junction_overlap[1]):]
                if (p_len > mjo[1]):
                    l_start = -(p_len-mjo[1]) # Should be a negative number most of the time
                else:
                    l_start = len(us_seq)
                l_seq = us_seq[l_start:]
                
                r_end = max(0, p_len-mjo[0])
                r_seq = ds_seq[:r_end]
                
                j_seq = l_seq + seq + r_seq
                
                if (len(j_seq) != len(seq)):
                    logging.info('    l_seq = ' + l_seq)
                    logging.info('      seq = ' + ' '*len(l_seq) + seq)
                    logging.info('    r_seq = ' + ' '*len(l_seq+seq) + r_seq)
                
                for i in range(len(j_seq)-p_len+1):
                    #pos = i-len(l_seq) # Position on sequence (can take negative number if hanging off the edge)
                    pos = start-len(l_seq)+i
                    pr = rc(j_seq[i:i+p_len])
                    
                    if cls.check_case(case, pr):
                        #primer_object = Primer.create(sequence=pr, position=pos, template_length=len(seq), strand='-', o_oligo=self)
                        #p = Primer(pr, gene, locus, genome, region, contig, strand, start+i, start+i+p_len, name=name)
                        p = Primer(pr, gene, locus, genome, region, contig, strand, pos, pos+p_len, name=name)
                    
                    n_count += 1
                    n_max += 1
                
                logging.info('    Primer scan: {}/~{}'.format(n_count, n_est))
        
        logging.info('    Primer scan: added {}/{} potential primers'.format(n_count, n_max))
    
    @classmethod
    def scan_seqs_only(cls, seq, side, primer_size=(18,26), us_seq='', ds_seq='', min_junction_overlap=(4,8)):
        """
        Pass in a sequence, then either design left/forward/+ primers
        or design right/reverse/- primers.
        Returns full list of primer sequences (not checked for quality).
        """
        # 'min_junction_overlap' allows for primers to span into the us/ds sequences
        #           uuuuuIIIIddddd
        #              └primer┘
        # 
        # min_junction_overlap = (
        #    minimum nt of primer that CANNOT be on the 5' flanking sequence (us_seq for forward primer),
        #    minimum nt of primer that CANNOT be on the 3' flanking sequence (ds_seq for forward primer)
        # )
        # 
        # Example of min_junction_overlap=(1,3) (for forward primer)
        #     us_seq    GTAAAGAAGG             10 nt
        #     seq                 CC            2 nt
        #     ds_seq                CCAAATTA    8 nt
        #                         ┌1 nt overlap
        #     primers         FFFFF 
        #                      FFFFF
        #                       FFFFF
        #     4 primers          FFFFF
        #                        └─┴3 nt overlap
        # Example:
        #       xxxxxxxxxx len(seq) = 10
        #   FFFFF FFFFF    len(primer) = 5
        #    FFFFF FFFFF   min_junction_overlap = (1,3)
        #     FFFFF FFFFF  12 primers
        #      FFFFF FFFFF
        #       FFFFF FFFFF
        #        FFFFF FFFFF
        #    ffff ffff ffff len(primer) = 4
        #     ffff ffff     11 primers
        #      ffff ffff
        #       ffff ffff
        #        ffff ffff
        # Code to temporarily limit number of primers computed set to None to disable
        #subset_size = 5000 # 1000 # 9000
        
        logging.info("  Preliminary primer sequence scan: Scanning '{}'".format(side))
        
        primer_list = []
        if (side in ['left', 'L', 'forward', 'F', '+']):
            for p_len in range(primer_size[0], primer_size[1]+1):
                # Include the flanking regions
                l_seq = us_seq[len(us_seq) - (p_len-min_junction_overlap[0]):]
                r_seq = ds_seq[:max(0, p_len-min_junction_overlap[1])]
                
                j_seq = l_seq + seq + r_seq
                
                for i in range(len(j_seq)-p_len+1):
                    pos = i-len(l_seq) # Position on sequence (can take negative number if hanging off the edge)
                    pf = j_seq[i:i+p_len]
                    primer_list.append(pf)
        
        elif (side in ['right', 'R', 'reverse', '-']):
            for p_len in range(primer_size[0], primer_size[1]+1):
                # Include the flanking regions
                l_seq = us_seq[len(us_seq) - (p_len-min_junction_overlap[1]):]
                r_seq = ds_seq[:max(0, p_len-min_junction_overlap[0])]
                
                j_seq = l_seq + seq + r_seq
                
                for i in range(len(j_seq)-p_len+1):
                    pos = i-len(l_seq) # Position on sequence (can take negative number if hanging off the edge)
                    pr = rc(j_seq[i:i+p_len])
                    primer_list.append(pr)
        
        logging.info('  Preliminary primer sequence scan: found {} primers'.format(len(primer_list)))
        
        return primer_list
    
    @classmethod
    def create(cls, sequence, position, template_length, strand, checks=True, name=None, o_oligo=None):
        """
        Add this to the non-redundant list of 'Primer.sequences' objects if it
        doesn't already exist.
        
        Returns the object
        """
        # Check to see if a 'Primer' object exists with this specific sequence
        p = cls.sequences.get(sequence, None)
        
        # If no 'Primer' object exists with this sequence, then create a new one, and add it
        if (p == None):
            p = cls(sequence, position, template_length, strand, checks=checks, name=name, o_oligo=o_oligo)
            cls.sequences[sequence] = p
        # Otherwise, since a 'Primer' object already exists with this specific sequence, add the new location to it
        else:
            #p.locations.add((gene, locus, genome, region, contig, strand, start, end))
            pass
        return p
    
    def get_tm(self):
        if (self.o_reverse_complement.__class__.__name__ == 'ThermoResult'):
            return self.o_reverse_complement.tm
        else:
            return min(self.o_reverse_complement).melting_temperature
    
    def get_min_delta_G(self):
        if ((self.o_hairpin != None) and (self.o_self_dimer != None)):
            if (self.o_hairpin.__class__.__name__ == 'ThermoResult'):
                return min(self.o_hairpin.dg, self.o_self_dimer.dg)/1000
            else:
                return min(self.o_hairpin + self.o_self_dimer).delta_G
        else:
            return -math.inf
    
    @classmethod
    def get_3prime_complementation_length(cls, seq1, seq2, max_3prime_search_length=5):
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
    
    def check(self, o_oligo):
        # Checks all attributes from scratch
        self.checks, self.o_hairpin, self.o_self_dimer, self.o_reverse_complement = self.check_primer_attributes(o_oligo=o_oligo)
        
        # Re-calculate weight
        self.weight = self.get_weight()
    
    def progressive_check(self, cutoffs):
        '''
        Checks current primer against the current cutoff parameters
        can resume from where it left off if the cutoffs were changed
        '''
        calc_weights = False
        for i, c in enumerate(self.checks):
            # If the check passed, then go to the next check
            if c:
                continue
            # If the current check was one that failed (False)
            # or hasn't been attempted (None), then do it
            # (because cutoff parameters may have changed)
            else: # ((c == False) or (c == None)):
                cc = self.check_index(i, cutoffs)
                self.checks[i] = cc
                
                # If the check succeeds, then mark weights to be re-calculated
                if cc:
                    calc_weights = True
                # If the check fails, then stop checking
                else:
                    break
        
        # Re-calculate weights if necessary
        if calc_weights:
            self.weight = self.get_weight()
    
    def check_index(self, index, args):
        '''
        args must be a dict containing the following:
            args = {
                mask                               'none',
                length:                            (17, 28),
                last5gc_count:                     (1, 3),
                gc_clamp_length:                   (0, 2),
                gc:                                (0.25, 0.75),
                max_run_length:                    4,
                max_3prime_complementation_length: 3,
                min_delta_g:                       -5.0, 
                folder:                            '/dev/shm/addtag',
                tm:                                (52,65),
                o_oligo:                           None,
            }
        '''
        functions = [
            self.check_mask,
            self.check_length,
            self.check_last5gc_count,
            self.check_gc_clamp_length,
            self.check_gc_freq,
            self.check_max_run_length,
            self.check_3prime_complementation_length,
            self.check_min_delta_g,
            self.check_tm,
        ]
        return functions[index](**args)
    
    def check_primer_attributes(self, thorough=False, o_oligo=None):
        o1, o2, o3 = None, None, None
        if thorough:
            checks = [
                self.check_length(),
                self.check_last5gc_count(),
                self.check_gc_clamp_length(),
                self.check_gc_freq(),
                self.check_max_run_length(),
                self.check_3prime_complementation_length(),
                None,
                None
            ]
            checks[6], o1, o2 = self.check_min_delta_g(o_oligo=o_oligo)
            checks[7], o3 = self.check_tm(o_oligo=o_oligo)
        else:
            checks = [None]*8
        
            checks[0] = self.check_length()
            if checks[0]:
                checks[1] = self.check_last5gc_count()
                if checks[1]:
                    checks[2] = self.check_gc_clamp_length()
                    if checks[2]:
                        checks[3] = self.check_gc_freq()
                        if checks[3]:
                            checks[4] = self.check_max_run_length()
                            if checks[4]:
                                checks[5] = self.check_3prime_complementation_length()
                                if checks[5]:
                                    checks[6], o1, o2 = self.check_min_delta_g(o_oligo=o_oligo)
                                    if checks[6]:
                                        checks[7], o3 = self.check_tm(o_oligo=o_oligo)
        
        return checks, o1, o2, o3
    
    def check_primer_attributes_old(self, thorough=False):
        length_passed = self.check_length()
        last5gc_count_passed = self.check_last5gc_count()
        gc_clamp_length_passed = self.check_gc_clamp_length()
        gc_freq_passed = self.check_gc_freq()
        max_run_length_passed = self.check_max_run_length()
        max_3prime_complementation_length_passed = self.check_3prime_complementation_length()
        
        results = [
            length_passed,
            last5gc_count_passed,
            gc_clamp_length_passed,
            gc_freq_passed,
            max_run_length_passed,
            max_3prime_complementation_length_passed,
            None, # min_delta_g_passed,
            None # tm_passed
        ]
        
        if (thorough or summarize(results)):
            # Only perform these two if the above 5 passed
            min_delta_g_passed, o1, o2 = self.check_min_delta_g()
            tm_passed, o3 = self.check_tm()
            
            results = [
                length_passed,
                last5gc_count_passed,
                gc_clamp_length_passed,
                gc_freq_passed,
                max_run_length_passed,
                max_3prime_complementation_length_passed,
                min_delta_g_passed,
                tm_passed
            ]
        return results, o1, o2, o3
    
    def check_mask(self, mask='none', **kwargs): # Code incomplete (does not deal with ambiguous characters)
        # Remove primers that have "masked" characters (depending on upper/lower-cases)
        if (mask != None):
            if (mask == 'none'):
                mask_passed = True
            elif (mask == 'lower'):
                mask_passed = self.sequence.isupper()
            elif (mask == 'upper'):
                mask_passed = self.sequence.islower()
        else:
            mask_passed = None
        return mask_passed
    
    def check_length(self, length=(17, 28), **kwargs):
        # Check length of primer
        # primer length should be 17-28 nt long
        if (length != None): # 'None' means skip this calculation
            length_passed = length[0] <= len(self.sequence) <= length[1]
        else:
            length_passed = None
        return length_passed
    
    def check_last5gc_count(self, last5gc_count=(1, 3), **kwargs):
        # Does the last 5 nt of the sequence have 1-3 C/G bases?
        # should avoid runs of 3-or-more Cs and Gs in 3' end
        if (last5gc_count != None): # 'None' means skip this calculation
            C_count = self.sequence[-5:].count('C')
            G_count = self.sequence[-5:].count('G')
            last5gc_count_passed = last5gc_count[0] <= C_count+G_count <= last5gc_count[1]
        else:
            last5gc_count_passed = None
        
        return last5gc_count_passed
    
    def check_gc_clamp_length(self, gc_clamp_length=(0, 2), **kwargs): # (1,2)
        # NOT iupac-aware
        # Does the 3' end of the sequence have the sequence 
        # 3' GC clamp
        if (gc_clamp_length != None): # 'None' means skip this calculation
            i = -1
            gcl = 0
            while (self.sequence[i] in ['G', 'C']):
                i -= 1
                gcl += 1
            
            gc_clamp_length_passed = gc_clamp_length[0] <= gcl <= gc_clamp_length[1]
        else:
            gc_clamp_length_passed = None
        
        return gc_clamp_length_passed
    
    def check_gc_freq(self, gc=(0.25, 0.75), **kwargs):
        # %GC should be between 40-60
        if (gc != None): # 'None' means skip this calculation
            gc_passed = gc[0] <= self.gc <= gc[1]
        else:
            gc_passed = None
            
        return gc_passed
    
    def check_max_run_length(self, max_run_length=4, **kwargs):
        # NOT iupac-aware
        
        # Primers with long runs of a single base should generally be avoided as they can misprime
        if (max_run_length != None): # 'None' means skip this calculation
            matches = list(regex.finditer(r'(.)\1*', self.sequence))
            max_run = 0
            for m in matches:
                max_run = max(max_run, len(m.group()))
            max_run_length_passed = max_run <= max_run_length
        else:
            max_run_length_passed = None
        
        return max_run_length_passed
    
    def check_3prime_complementation_length(self, max_3prime_complementation_length=3, **kwargs):
        # 3' ends of primers should NOT be complementary
        # As even a little dimerization will inhibit target annealing
        if (max_3prime_complementation_length != None):
            max_3prime_complementation_length_passed = self.get_3prime_complementation_length(self.sequence, self.sequence, max_3prime_search_length=max_3prime_complementation_length+1) <= max_3prime_complementation_length
        else:
            max_3prime_complementation_length_passed = None
            
        return max_3prime_complementation_length_passed
    
    def check_min_delta_g(self, min_delta_g=-5.0, folder='/dev/shm/addtag', o_oligo=None, **kwargs):
        #o1 = None
        #o2 = None
        
        # min delta-G should be -3.0
        if ((min_delta_g != None) and (o_oligo != None)): # 'None' means skip this calculation
            # Calculate hairpin delta-G
            if (self.o_hairpin == None):
                #o1 = Structure.new_calculate_simple(folder, self.sequence)
                self.o_hairpin = o_oligo.find_structures(folder, self.sequence)
            
            # Calculate homodimer delta-G
            if (self.o_self_dimer == None):
                #o2 = Structure.new_calculate_simple(folder, self.sequence, self.sequence)
                self.o_self_dimer = o_oligo.find_structures(folder, self.sequence, self.sequence)
            
            min_delta_g_passed = min_delta_g <= min(min(self.o_hairpin).delta_G, min(self.o_self_dimer).delta_G)
        else:
            min_delta_g_passed = None
        
        #self.o_hairpin, self.o_self_dimer = o1, o2
        #return min_delta_g_passed, o1, o2
        return min_delta_g_passed
    
    def check_tm(self, tm=(52,65), folder='/dev/shm/addtag', o_oligo=None, **kwargs):
        #o3 = None
        #logging.info('      check_tm()')
        # Tm (against reverse_complement) should be between 50-70 or 55-65
        if ((tm != None) and (o_oligo != None)): # 'None' means skip this calculation
            #logging.info('        if(True)')
            # Calculate reverse-complement delta-G and Tm
            if (self.o_reverse_complement == None):
                #logging.info('          calculating o_reverse_complement')
                #o3 = Structure.new_calculate_simple(folder, self.sequence, rc(self.sequence))
                self.o_reverse_complement = o_oligo.find_structures(folder, self.sequence, rc(self.sequence))
        
            tm_passed = tm[0] <= min(self.o_reverse_complement).melting_temperature <= tm[1]
        else:
            tm_passed = None
        
        #self.o_reverse_complement = o3
        #return tm_passed, o3
        return tm_passed
    
    @classmethod
    def summarize(cls, results):
        return all(x for x in results if (x != None))
    
    def align(self, sequence):
        """
        Align the primer input 'sequence' on both strands.
        Returns list of all alignments:
          [(start, end strand), ...]
        """
        # Search for primer alignments to 'sequence' on both '+' and '-' strands
        locations = []
        
        # Search for '+' primer
        for m in regex.finditer(self.sequence, sequence, flags=regex.IGNORECASE, overlapped=True):
            locations.append((m.start(), m.end(), '+'))
        
        # Search for '-' primer
        for m in regex.finditer(rc(self.sequence), sequence, flags=regex.IGNORECASE, overlapped=True):
            locations.append((m.start(), m.end(), '-'))
        
        return locations
    
    def __repr__(self):
        labs = ['name', 'seq', 'pos', 'strand', 'Tm', 'GC', 'min(dG)']
        vals = [
            self.get_name(),
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
    pairs = {}
    
    def __init__(self, forward_primer, reverse_primer):
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        
        self.weight = 0.0
        self.checks = [None]*4
        self.o_heterodimer = None
        
        # Only add to dict if this PrimerPair has correct orientation/positioning
        if (len(self.get_amplicon_sizes()) > 0):
            # Add this to the non-redundant 'PrimerPair.pairs' dict if it isn't already present.
            the_pair = self.pairs.setdefault((self.forward_primer.sequence, self.reverse_primer.sequence), self)
        
        # Only do expensive calculations if this is a pair that hasn't been added yet
        #if (self == the_pair):
        #    self.weight = self.get_weight(*args, **kwargs)
        #    # self.check()
    
#    def get_locations(self):
#        self.forward_primer.locations
#        self.reverse_primer.locations
#        if any((locus, f_reg, f_ori) == loc[0:2]+loc[3:4] for loc in p.locations):
    
    def __lt__(self, other):
        '''
        Sort 'PrimerPair' objects by their joint weights
        '''
        return self.get_joint_weight() < other.get_joint_weight()
    
    @classmethod
    def create(cls, forward_primer_object, reverse_primer_object, *args, checks=True, intervening=0, o_oligo=None, **kwargs):
        """
        Add this to the non-redundant list of 'PrimerPair.pairs' objects if it
        doesn't already exist.
        
        Returns the object
        """
        ii = (forward_primer_object.sequence, reverse_primer_object.sequence)
        pp = cls.pairs.get(ii, None)
        if (pp == None):
            pp = cls(forward_primer_object, reverse_primer_object, *args, checks=checks, intervening=intervening, o_oligo=o_oligo, **kwargs)
            cls.pairs[ii] = pp
        
        return pp
    
    def check(self, o_oligo, locus, f_region, r_region, contig, minimize):
        self.checks, self.o_heterodimer = self.check_pair_attributes(o_oligo=o_oligo)
        self.weight = self.get_weight(locus, f_region, r_region, contig, minimize=minimize)
    
    def progressive_check(self, cutoffs):
        '''
        Checks PrimerPair against the current cutoff parameters
        can resume from where it left off if the cutoffs were changed.
        
        Assumes that 'self.forward_primer' and 'self.reverse_primer' have completed their checks.
        '''
        calc_weights = False
        for i, c in enumerate(self.checks):
            # If the check passed, then go to the next check
            if c:
                continue
            # If the current check was one that failed (False)
            # or hasn't been attempted (None), then do it
            # (because cutoff parameters may have changed)
            else: # ((c == False) or (c == None)):
                cc = self.check_index(i, cutoffs)
                self.checks[i] = cc
                
                # If the check succeeds, then mark weights to be re-calculated
                if cc:
                    calc_weights = True
                # If the check fails, then stop checking
                else:
                    break
        
        # Re-calculate weights if necessary
        if calc_weights:
            self.weight = self.get_weight()
    
    def calc_weight(self, locus, f_region, r_region, contig, kwargs):
        self.weight = self.get_weight(locus, f_region, r_region, contig, **kwargs) # pass in the 'minimize=True', 'locus' and 'contig' arguments
    
    def check_index(self, index, args):
        '''
        args must be a dict containing the following:
            args = {
                amplicon_size_range:               (300,700),
                max_tm_difference:                 3.0,
                max_3prime_complementation_length: 3,
                min_delta_g:                       -5.0, 
                folder:                            '/dev/shm/addtag',
                o_oligo:                           None,
            }
        '''
        functions = [
            self.check_amplicon_size,
            self.check_max_tm_difference,
            self.check_3prime_complementation_length,
            self.check_min_delta_g,
        ]
        return functions[index](**args)
    
    def check_pair_attributes(self, thorough=False, o_oligo=None, **kwargs):
        o4 = None
        if thorough:
            checks = [
                self.check_amplicon_size(),
                self.check_max_tm_difference(),
                self.check_3prime_complementation_length(),
                None
            ]
            checks[3], o4 = self.check_min_delta_g(o_oligo=o_oligo)
        else:
            checks = [None]*4
        
            checks[0] = self.check_amplicon_size(**kwargs)
            if checks[0]:
                checks[1] = self.check_max_tm_difference()
                if checks[1]:
                    checks[2] = self.check_3prime_complementation_length()
                    if checks[2]:
                        checks[3], o4 = self.check_min_delta_g(o_oligo=o_oligo)
        
        return checks, o4
    
    def check_amplicon_size(self, amplicon_size_range=(300,700), **kwargs): # locus='', contig=''
        if (amplicon_size_range != None):
            amplicon_size_passed = any(amplicon_size_range[0] <= x[4] <= amplicon_size_range[1] for x in self.get_amplicon_sizes()) # (locus, contig, size)
        else:
            amplicon_size_passed = None
        return amplicon_size_passed
    
    def check_max_tm_difference(self, max_tm_difference=3.0, **kwargs):
        # The difference in Tms should be as small as possible
        if (max_tm_difference != None):
            tm1 = min(self.forward_primer.o_reverse_complement).melting_temperature
            tm2 = min(self.reverse_primer.o_reverse_complement).melting_temperature
            
            max_tm_difference_passed = abs(tm1-tm2) <= max_tm_difference
        else:
            max_tm_difference_passed = None
        return max_tm_difference_passed
    
    def check_3prime_complementation_length(self, max_3prime_complementation_length=3, **kwargs):
        # 3' ends of primers should NOT be complementary
        # As even a little dimerization will inhibit target annealing
        if (max_3prime_complementation_length != None):
            max_3prime_complementation_length_passed = self.forward_primer.get_3prime_complementation_length(self.forward_primer.sequence, self.reverse_primer.sequence, max_3prime_search_length=max_3prime_complementation_length+1) <= max_3prime_complementation_length
        else:
            max_3prime_complementation_length_passed = None
            
        return max_3prime_complementation_length_passed
    
    def check_min_delta_g(self, min_delta_g=-5.0, folder='/dev/shm/addtag', o_oligo=None, **kwargs):
        #o4 = None
        
        # min delta-G should be -3.0
        if ((min_delta_g != None) and (o_oligo != None)): # 'None' means skip this calculation
            # Calculate heterodimer delta-G
            if (self.o_heterodimer == None):
                #o4 = Structure.new_calculate_simple(folder, self.forward_primer.sequence, self.reverse_primer.sequence)
                self.o_heterodimer = o_oligo.find_structures(folder, self.forward_primer.sequence, self.reverse_primer.sequence)
            
            min_delta_g_passed = min_delta_g <= min(self.o_heterodimer).delta_G
        else:
            min_delta_g_passed = None
        
        return min_delta_g_passed #, o4
    
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
    
    def get_real_amplicon_sizes(self, template):
        """
        Search template for amplification of forward and reverse primers.
        """
        return [len(x) for x in self.amplify(template)]
    
    def get_amplicon_size(self, gene, locus, genome, region_f, region_r, contig):
        # Go through the locations within each 'Primer' object
        # Return the desired amplicon size
        #location = (gene, locus, genome, region, contig, strand, start, end)
        for f_gene, f_locus, f_genome, f_region, f_contig, f_strand, f_start, f_end in self.forward_primer.locations:
            for r_gene, r_locus, r_genome, r_region, r_contig, r_strand, r_start, r_end in self.reverse_primer.locations:
                if (gene == f_gene == r_gene):
                    if (locus == f_locus == r_locus):
                        if (genome == f_genome == r_genome):
                            if (contig == f_contig == r_contig):
                                if ((region_f == f_region) and (region_r == r_region)):
                                    if ((f_strand == '+') and (r_strand == '-') and (f_start < r_end)):
                                        return r_end - f_start
                                    elif ((f_strand == '-') and (r_strand == '+') and (r_start < f_end)):
                                        return f_end - r_start
        return 0
    
    def get_amplicon_sizes(self):
        # Go through the locations within each 'Primer' object
        # Return the desired amplicon size
        sizes = []
        for f_gene, f_locus, f_genome, f_region, f_contig, f_strand, f_start, f_end in self.forward_primer.locations:
            for r_gene, r_locus, r_genome, r_region, r_contig, r_strand, r_start, r_end in self.reverse_primer.locations:
                if ((f_gene == r_gene) and (f_locus == r_locus) and (f_genome == r_genome) and (f_contig == r_contig)): # and (f_region == r_region) # Region doesn't matter here
                    if ((f_strand == '+') and (r_strand == '-') and (f_start < r_end)):
                        sizes.append((f_gene, f_locus, f_genome, f_contig, r_end - f_start))
                    elif ((f_strand == '-') and (r_strand == '+') and (r_start < f_end)):
                        sizes.append((f_gene, f_locus, f_genome, f_contig, f_end - r_start))
        return sizes
    
    def get_amplicon_size_old(self):
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
        seqs = []
        for f_gene, f_locus, f_genome, f_region, f_contig, f_strand, f_start, f_end in self.forward_primer.locations:
            for r_gene, r_locus, r_genome, r_region, r_contig, r_strand, r_start, r_end in self.reverse_primer.locations:
                if ((f_gene == r_gene) and (f_locus == r_locus) and (f_genome == r_genome) and (f_contig == r_contig)): # and (f_region == r_region) # Region doesn't matter here
                    if ((f_strand == '+') and (r_strand == '-') and (f_start < r_end)):
                        interspace = r_start-f_end
                        seqs.append((f_locus, f_contig, self.forward_primer.sequence + '·'*interspace + rc(self.reverse_primer.sequence)))
                    elif ((f_strand == '-') and (r_strand == '+') and (r_start < f_end)):
                        interspace = f_start-r_end
                        seqs.append((f_locus, f_contig, self.reverse_primer.sequence + '·'*interspace + rc(self.forward_primer.sequence)))
        return seqs
        
    #    ####### This needs to be re-done (it isn't correct) ####### <---- I made changes, so it needs to be re-checked
    #    #interspace = self.reverse_primer.position - (self.intervening + self.forward_primer.position + len(self.forward_primer.sequence))
    #    interspace = self.forward_primer.template_length - self.forward_primer.position - len(self.forward_primer.sequence) + self.intervening + self.reverse_primer.position
    #    return self.forward_primer.sequence + '·'*interspace + rc(self.reverse_primer.sequence)
    #    #print(' '*self.forward_primer.position + fp.sequence)
    #    #print(' '*self.reverse_primer.position + rc(rp.sequence))
    
    def get_weight(self, gene=None, locus=None, genome=None, f_region=None, r_region=None, contig=None, minimize=False, **kwargs):
        """
        Use logistic multiplier to weigh the PrimerPair
        """
        w = 1.0
        
        # Weigh by amplicon size
        # If this is supposed to be a 'sF' 'sR' pair, then we want to minimize the amplicon size
        # Otherwise, we want the amplicon size to be around 500 nt
        #if (("minimize" in kwargs) and kwargs["minimize"]):
        
        if ((gene != None) and (locus != None) and (genome != None) and (f_region != None) and (r_region != None) and (contig != None)):
            if minimize:
                w *= logistic_down(self.get_amplicon_size(gene, locus, genome, f_region, r_region, contig), downslope=1.03, down=150)
            else:
                w *= logistic_updown(self.get_amplicon_size(gene, locus, genome, f_region, r_region, contig), 1.03, 300, 1.03, 700)
        else:
            amps = self.get_amplicon_sizes()
            amp_sizes = [x[4] for x in amps]
            if minimize:
                f_sizes = []
                for a in amps:
                    m_list = minimize.get(a[:4], None)
                    if m_list:
                        for m in m_list:
                            f_sizes.append(a[4]-m)
                # No guard against division by zero
                w *= logistic_down(sum(f_sizes)/len(f_sizes), downslope=1.03, down=150)
            else:
                w *= logistic_updown(sum(amp_sizes)/len(amp_sizes), 1.03, 300, 1.03, 700)
            
            # Previous code
            #    w *= logistic_down(sum(sizes)/len(sizes), downslope=1.03, down=150)
            #else:
            #    w *= logistic_updown(sum(sizes)/len(sizes), 1.03, 300, 1.03, 700)
        
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
    
    def amplify(self, sequence):
        """
        Search input sequence for locations where 'forward' and 'reverse'
        primers can bind, and return all putative amplification sequences.
        
        Primer matches must be exact.
        
        Does not work with circular sequences
        """
        # Search 'forward' primer on both '+' and '-' strands
        locations = self.forward_primer.align(sequence) + self.reverse_primer.align(sequence)
        
        plus_list = []
        minus_list = []
        for loc in locations:
            if (loc[2] == '+'):
                plus_list.append(loc)
            else:
                minus_list.append(loc)
        
        results = []
        for pl in plus_list:
            for ml in minus_list:
                if (pl[0] < ml[1]):
                    results.append(sequence[pl[0]:ml[1]])
        
        return results
        
    # INCOMPLETE
    # def in_silico_pcr_from_file(self, list_of_fasta_templates=(), basename=False):
    #     results = []
    #     for filename in list_of_fasta_templates:
    #         contigs = utils.old_load_fasta_file(filename)
    #         
    #         fn = filename
    #         if basename:
    #             fn = os.path.basename(filename)
    #         
    #         for cname, cseq in contigs.items():
    #             amp_list = self.amplify(cseq)
    #             results.extend(amp_list)
    #     
    #     return [len(r) for r in results]
    
    def in_silico_pcr(self, contigs):
        """
        Returns list of expected amplicon sizes if the PrimerPair
        was used to amplify 'genome' as a template
        """
        
        results = []
        
        # contigs should be a dict with key=contig name, value=sequence
        if contigs:
            for cname, cseq in contigs.items():
                amp_list = self.amplify(cseq)
                results.extend([len(a) for a in amp_list])
        
        return results
    
    def __repr__(self):
        labs = ['names', 'seqs', 'amplicon_sizes', 'Tms', 'GCs', 'min(dG)']
        vals = [
            (self.forward_primer.get_name(), self.reverse_primer.get_name()),
            (self.forward_primer.sequence, self.reverse_primer.sequence),
            self.get_amplicon_sizes(),
            tuple(round(x, 2) for x in self.get_tms()),
            self.get_gcs(),
            round(self.get_min_delta_G(), 2)
        ]
            
        return self.__class__.__name__ + '(' + ', '.join('='.join(map(str, x)) for x in zip(labs, vals)) + ')'

class PrimerSet(object):
    """
    Holds a list of primers and their weight
    """
    def __init__(self, primer_pair_list=()):
        self.pp_list = [pp for pp in primer_pair_list]
        self.weight = self.calculate_weight()
    
    def calculate_weight(self):
        '''
        Calculate joint weight of this specific pp_list
        '''
        p_group_weight = self.p_group_weight(self.flatten_primer_pair_list(self.pp_list))
        pp_group_weight = self.pp_group_weight(self.pp_list)
        joint_weight = p_group_weight * pp_group_weight
        return joint_weight
    
    def p_group_weight(self, primers, merge_duplicates=True):
        """
        Weigh by difference of each primer's Tm and the mean Tm
        """
        w = 1.0
        
        # First, we weed out any None, and take average Tm if a Primer is not unique
        if merge_duplicates:
            pdict = {}
            for p in primers:
                if p:
                    pdict.setdefault(p.sequence, []).append(p.get_tm())
            temps = []
            for pseq, tms in pdict.items():
                temps.append(sum(tms)/len(tms))
        else:
            temps = [p.get_tm() for p in primers if p] # will not include any p == None
        
        if (len(temps) > 0):
            mean_temp = sum(temps)/len(temps) # mean Tm
            #print('temps =', temps, file=sys.stderr)
            for i, tm in enumerate(temps):
                w *= logistic_updown(tm-mean_temp, 5, -2.5, 5, 2.5)
            
            return w
        else:
            return 0.0
    
    def pp_group_weight(self, primer_pairs, merge_duplicates=True):
        """
        Weigh a list of PrimerPair objects.
        """
        weight = 1.0
        
        if merge_duplicates:
            ppdict = {}
            for pp in primer_pairs:
                if pp:
                    k = (pp.forward_primer.sequence, pp.reverse_primer.sequence)
                    ppdict.setdefault(k, []).append(pp.get_joint_weight())
            w_list = []
            for k, pp_w in ppdict.items():
                w_list.append(sum(pp_w)/len(pp_w))
            
        else:
            w_list = [x.get_joint_weight() for x in primer_pairs if x]
        
        for i, w in enumerate(w_list):
            weight *= w
        
        return weight
    
    def flatten_primer_pair_list(self, pair_list):
        p_set = []
        for pp in pair_list:
            if pp:
                p_set.append(pp.forward_primer)
                p_set.append(pp.reverse_primer)
            else:
                p_set.append(None)
                p_set.append(None)
        return p_set
    
    def __repr__(self):
        return self.__class__.__name__ + '(weight=' + str(self.weight) + ', ' + str(self.pp_list) + ')'

class PrimerDesign(object):
    """
    The 'PrimerDesign' class is to make the simulated annealing (Markov chain Monte-Carlo)
    more simple to execute.
    
    The argument to '__init__()' is the 2D list of 'PrimerPair' objects.
    Then the 'optimize()' command runs the simulated annealing until either
    convergence or a certain number of iterations have elapsed.
    
    The final design is a single list of 'PrimerPair' objects that has the
    greatest 'weight' of all simulated annealings.
    
    Example usage:
      pp_2d_list = [
        [PrimerPair(), PrimerPair(), PrimerPair()],
        [PrimerPair(), PrimerPair(), PrimerPair()]
      ]
      design = PrimerDesign(pp_2d_list)
      design.optimize()
      optimal = design.optimal
    """
    def __init__(self, primer_pair_list):
        self.primer_pair_list = primer_pair_list
        self.optimal = None # Store 'PrimerSet' object
        self.done = False
    
    def random_choices(self, population, weights, k=1):
        """
        Return a k sized list of population elements chosen with replacement.
        """
        weight_sum = sum(weights)
        choices = zip(population, weights)
        values = []
        for i in range(k):
            r = random.uniform(0, weight_sum)
            upto = 0
            for c, w in choices:
                if upto + w >= r:
                    values.append(c)
                    break
                upto += w
            else:
                values.append(random.choice(population))
        return values
    
    def random_primer_by_weight(self, pairs):
        if (len(pairs) == 0):
            return None
        elif ('choices' in random.__all__):
            return random.choices(pairs, [x.get_joint_weight() for x in pairs])[0]
        else:
            return self.random_choices(pairs, [x.get_joint_weight() for x in pairs])[0]
    
    def rank_order(self, x, reverse=False, shift=0):
        return [y+shift for y in sorted(range(len(x)), key=x.__getitem__, reverse=reverse)]
    
    def anneal(self, current=None, mode='ranked'):
        """
        Returns 'PrimerSet' object which contains the weight and pp_list.
        
        This will just be a simple, random sample of the potential PrimerPair objects.
        """
        weight = 1
        
        if (mode == 'uniform'):
            # Get random set based on uniform distribution
            pp_set = [random.choice(pp_list) if (len(pp_list) > 0) else None for pp_list in self.primer_pair_list]
        elif (mode == 'ranked'):
            # Get ranked-random set
            pp_set = [self.random_primer_by_weight(pp_list) for pp_list in self.primer_pair_list]
        elif (mode == 'best'):
            # Get top-ranked set, assuming they are sorted
            pp_set = [pp_list[0] if (len(pp_list) > 0) else None for pp_list in self.primer_pair_list]
        elif (mode == 'direct'):
            # Replace the worst-weighted 'PrimerPair'
            # If no replacement would improve, then do nothing (return the current)
            
            # If no swap would improve the weight, then choose a random one
            pp_set = self.anneal(mode='ranked').pp_list
            
            if (current != None):
                if not self.done:
                    # Get list of indices from smallest weight to largest weight (excluding 'sF' and 'sR')
                    # If a primer is 'None', then it is right-most in the list.
                    # Left-most element is the smallest weight. Weights increase as the index of the list increases.
                    # If primer is 'None', then its weight will be 'math.inf', and the index
                    # corresponding to it will be the last element in wi_order:
                    #   rank_order([1, 2, 2.5, 0.001, math.inf]) # [3, 0, 1, 2, 4]
                    wi_order = self.rank_order([pp.get_joint_weight() if pp else math.inf for pp in current.pp_list])
                    
                    # Start at the worst index (wi), and go to the next-worst, then next, then next
                    for wi in wi_order:
                        # Shallow copy the list of primer pairs
                        pp_set_d = PrimerSet(current.pp_list[:])
                        
                        # Swap the worst-performing primer with an alternative.
                        # If the alternative gives a better weight, then keep it
                        # and break out of the loop
                        for source_pp in self.primer_pair_list[wi]:
                            pp_set_d.pp_list[wi] = source_pp
                            pp_set_d.weight = pp_set_d.calculate_weight()
                            
                            weight_delta = pp_set_d.weight - current.weight
                            
                            if (weight_delta > 0):
                                pp_set = pp_set_d.pp_list
                                break
                        else:
                            # If no pp swap is higher-weighted, then return 0 rather than the last-calculated 'weight_delta'
                            weight_delta = 0
                        
                        # If the pp swap already produces a higher-weighted primer set, then skip the rest of the 'wi'
                        if (weight_delta > 0):
                            break
                    
                    else:
                        # If no swap would improve the weight, then return the
                        # current, and ignore further calls to this function
                        self.done = True
                        return current
                else:
                    return current
        
        design = PrimerSet(pp_set)
        
        return design
    
    def optimize(self, mode='ranked', iterations=1000):
        # Define some benchmarking variables for parameter tuning
        accepts = 0
        rejects = 0
        improves = 0
        retrogresses = 0
        optimizes = 0
        ignores = 0
        
        # Copy the current optimal
        optimal = self.optimal
        
        self.done = False
        
        # Start out at a random place
        current = self.anneal(mode=mode)
        
        # Set minimum and maximum weight values for use with exponential cooling/heating
        weight_max = 1.0 # Starting value
        weight_min = 0.0001 # Ending value
        
        weight_factor = -math.log(weight_max/weight_min) # (log base 'e' by default)
        
        #logging.info('\t'.join(['i', 'threshold', 'dW', 'weight_threshold']))
        
        for i in range(iterations):
            # As 'i' approaches 'iterations', W decreases
            weight_threshold = weight_max * math.exp(weight_factor * i / iterations)
            
            # Perform the anneal
            a = self.anneal(current=current, mode=mode)
            
            dW = a.weight - current.weight
            
            # If the new combination gives a lower weight,
            # and the majority of the time (random)
            # Reject it
        #    try:
        #        threshold = math.exp(-dW/weight_threshold)
        #    except OverflowError:
        #        # print('   weight_min = {}'.format(weight_min), file=sys.stderr)
        #        # print('   weight_max = {}'.format(weight_max), file=sys.stderr)
        #        # print('weight_factor = {}'.format(weight_factor), file=sys.stderr)
        #        # print('           cW = {}'.format(cW), file=sys.stderr)
        #        # print('           dW = {}'.format(dW), file=sys.stderr)
        #        # sys.exit('ERROR')
        #        threshold = math.inf
            
            # Update performance variables
            if (dW < 0.0):
                retrogresses += 1
            elif (dW > 0.0):
                improves += 1
            elif (dW == 0.0):
                ignores += 1
            
            if ((dW < 0.0) and (random.random() > weight_threshold)):
                # Then reject this new state, and keep the current state unchanged
                rejects += 1
            else: # A small percent of the retrogressions will be accepted
                # Accept this new state
                accepts += 1
                
                # Replace current with new annealing
                current = a
                
                # Compare to the best state
                # If current is optimal, then replace the optimal
                if optimal:
                    if (current.weight > optimal.weight):
                        optimal = current
                        optimizes += 1
                else:
                    optimal = current
                    optimizes += 1
            
            #logging.info('\t'.join(map(str, [i, threshold, dW, weight_threshold])))
        
        logging.info('     accepts = {}'.format(accepts))
        logging.info('     rejects = {}'.format(rejects))
        logging.info('     ignores = {}'.format(ignores))
        logging.info('    improves = {}'.format(improves))
        logging.info('retrogresses = {}'.format(retrogresses))
        logging.info('   optimizes = {}'.format(optimizes))
        logging.info('  iterations = {}'.format(iterations))
        logging.info('     optimal = {}'.format(optimal))
        
        self.optimal = optimal
    
    def log_pp_union(self, pp_2d_list):
        """
        Write to log file which PP are shared
        
              0     1     2     3    F
          0  ---   0/5   0/2   3/10
          1  0/5   ---   0/4   0/6
          2  0/2   0/4   ---   0/3
          3  4/9   0/6   0/3   ---
          R
          
        The upper-right half correspond to the forward/+ Primer within each pair
        The lower-left half are reverse/- Primers within each pair
        Displayed values are the len(intersect)/len(union)
        """
        outputs = {}
        
        for i1, pp1 in enumerate(pp_2d_list):
            f1_set = set([x.forward_primer.sequence for x in pp1 if x])
            r1_set = set([x.reverse_primer.sequence for x in pp1 if x])
            
            for i2, pp2 in enumerate(pp_2d_list):
                f2_set = set([x.forward_primer.sequence for x in pp2 if x])
                r2_set = set([x.reverse_primer.sequence for x in pp2 if x])
                
                if (i1 < i2): # Analyze 'F' sets
                    top = len(f1_set.intersection(f2_set))
                    bot = len(f1_set.union(f2_set))
                    outputs[(i1, i2)] = '{:>3}/{:<3}'.format(top, bot)
                elif (i1 > i2): # Analyze 'R' sets
                    top = len(r1_set.intersection(r2_set))
                    bot = len(r1_set.union(r2_set))
                    outputs[(i1, i2)] = '{:>3}/{:<3}'.format(top, bot)
                else: # Diagonal
                    outputs[(i1, i2)] = '{:>3}/{:<3}'.format('-', '-')
        
        lines = []
        lines.append('{:>3}'.format(''))
        for i in range(len(pp_2d_list)):
            lines[-1] += ' {:^7}'.format(i)
        lines[-1] += ' F'
        
        for i1 in range(len(pp_2d_list)):
            lines.append('{:>3}'.format(i1))
            for i2 in range(len(pp_2d_list)):
                lines[-1] += ' ' + outputs[(i1, i2)]
        lines.append('{:>3}'.format('R'))
        
        logging.info('\n'.join(lines))
    
            # # Report the number of PrimerPairs that are shared among the inserts
            # t_num = len(pp_sources)//2
            # for t1 in range(t_num):
            #     t1L_set = set([(x.forward_primer.sequence, x.reverse_primer.sequence) for x in pp_sources[2*t1] if x])
            #     t1R_set = set([(x.forward_primer.sequence, x.reverse_primer.sequence) for x in pp_sources[(2*t1)+1] if x])
            #     for t2 in range(t_num):
            #         if (t1 < t2):
            #             t2L_set = set([(x.forward_primer.sequence, x.reverse_primer.sequence) for x in pp_sources[2*t2] if x])
            #             t2R_set = set([(x.forward_primer.sequence, x.reverse_primer.sequence) for x in pp_sources[(2*t2)+1] if x])
            #             tL_num = len(t1L_set.intersection(t2L_set))
            #             tL_den = len(t1L_set.union(t2L_set))
            #             tR_num = len(t1R_set.intersection(t2R_set))
            #             tR_den = len(t1R_set.union(t2R_set))
            #             logging.info('number shared insert PrimerPairs ({} vs {}): L={}/{}, R={}/{}'.format(t1, t2, tL_num, tL_den, tR_num, tR_den))
        
        
    
    def log_pp_identical(self, pp_list):
        """
        Does all pairwise comparisons of PrimerPairs within list 'pp_set'
        to identify which pairs have identical sequences.
        
              0 1 2 3 F
              P n n P
          0 P - F F T  
          1 n F - N F
          2 n F N - F
          3 P T F F -
          R
        
        T/True means the primers are identical
        F/False means the primers are diffetent
        N/None means there are no primers to compare
        The upper-right half correspond to the forward/+ Primer within each pair
        The lower-left half are reverse/- Primers within each pair
        """
        
        outputs = {}
        
        for i1, pp1 in enumerate(pp_list):
            for i2, pp2 in enumerate(pp_list):
                if (i1 < i2):
                    if pp1 and pp2:
                        if (pp1.forward_primer.sequence == pp2.forward_primer.sequence):
                            outputs[(i1, i2)] = 'T'
                        else:
                            outputs[(i1, i2)] = 'F'
                    elif pp1 or pp2:
                        outputs[(i1, i2)] = 'f'
                    else:
                        outputs[(i1, i2)] = 'N'
                elif (i1 > i2):
                    if pp1 and pp2:
                        if (pp1.reverse_primer.sequence == pp2.reverse_primer.sequence):
                            outputs[(i1, i2)] = 'T'
                        else:
                            outputs[(i1, i2)] = 'F'
                    elif pp1 or pp2:
                        outputs[(i1, i2)] = 'f'
                    else:
                        outputs[(i1, i2)] = 'N'
                else:
                    outputs[(i1, i2)] = '-'
        
        lines = []
        
        # First row (pp indices)
        lines.append('{:>3} {}'.format('', ' '))
        for i in range(len(pp_list)):
            lines[-1] += ' {:>2}'.format(i)
        lines[-1] += ' F'
        
        # Second row ('P' or 'n')
        lines.append('{:>3} {}'.formaT('', ' '))
        for i, pp in enumerate(pp_list):
            if pp:
                lines[-1] += ' P'
            else:
                lines[-1] += ' n'
        
        for i1, pp1 in enumerate(pp_list):
            pp1s = 'n'
            if pp1:
                pp1s = 'P'
            lines.append('{:>3} {}'.format(i1, pp1s))
            for i2 in range(len(pp_list)):
                lines[-1] += ' ' + outputs[(i1, i2)]
        lines.append('{:>3}'.format('R'))
        
        logging.info('\n'.join(lines))
    
    @classmethod
    def make_non_redundant(cls, set_list, sort=True):
        # Make Amp-D pp list non-redundant
        ppdict = {}
        for pp_set in set_list:
            w, pp_list = pp_set
            k = []
            for pp in pp_list:
                if pp:
                    k.append((pp.forward_primer.sequence, pp.reverse_primer.sequence))
                else:
                    k.append((None, None))
            k = tuple(k)
            
            d_pp_set = ppdict.setdefault(k, pp_set)
            if (d_pp_set[0] < w):
                ppdict[k] = pp_set
        
        if sort:
            return sorted(ppdict.values())
        else:
            return list(ppdict.values())
