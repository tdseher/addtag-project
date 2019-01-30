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
        
        # Code to temporarily limit number of primers computed set to None to disable
        #subset_size = 5000 # 1000 # 9000
        
        time_expired = False
        start_time = time.time() # Time in seconds
        
        # Total number of primers that will be scanned
        n_est = max(0, sum(len(seq)-x+1 for x in range(primer_size[0], primer_size[1]+1)))
        n_max = 0
        n_count = 0
        
        logging.info("  Primer scan: Scanning '{}'".format(side))
        
        # First, create the primers
        # Second, filter the Primer objects
        # Third, return only the good Primer objects
        good_primers = []
        if (side in ['left', 'forward', '+']):
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
        
        elif (side in ['right', 'reverse', '-']):
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
                
                pp_object = PrimerPair.create(cgf, cgr, intervening=intervening, o_oligo=self)
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
        
        
        mean_temp = sum(temps)/len(temps) # mean Tm
        #print('temps =', temps, file=sys.stderr)
        for i, tm in enumerate(temps):
            w *= logistic_updown(tm-mean_temp, 5, -2.5, 5, 2.5)
        
        return w
    
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
    
    def __init__(self, sequence, position, template_length, strand, checks=True, name=None, o_oligo=None):
        self.sequence = sequence
        
        
        # Code to store where this Primer is on the different genomes.
    #    location = (contig, orientation, start, end)
    #    self.locations = set()
        # Add this to the non-redundant list of 'Primer' objects if it doesn't already exist
        # otherwise, add this location to the pre-existing one
        #the_primer = self.sequences.setdefault(self.sequence, self)
    #    the_primer.locations.add(location)
        
        # Only do expensive calculations if this is a sequence that hasn't
        # been observed yet
        #if (self == the_primer):
        self.position = position
        self.template_length = template_length
        self.strand = strand
        
        self.gc = get_gc_freq(self.sequence)
        
        if checks:
            self.checks, self.o_hairpin, self.o_self_dimer, self.o_reverse_complement = self.check_primer_attributes(o_oligo=o_oligo)
        else:
            self.checks = [None]*8
            self.o_hairpin = None
            self.o_self_dimer = None
            self.o_reverse_complement = None
        
        self.weight = self.get_weight()
        self.name = name
    
    def __lt__(self, other):
        #if isinstance(other, Primer):
        return self.weight < other.weight
        #else:
        #    return False
        # __gt__(), __le__(), __ne__(), __ge__(), __eq__()
    
    # def get_gc(self):
    #     if (self.gc != None):
    #         return self.gc
    #     else:
    #         C_count = seq.count('C')
    #         G_count = seq.count('G')
    #         return (C_count+G_count)/len(seq)
    
    @classmethod
    def create(cls, sequence, position, template_length, strand, checks=True, name=None, o_oligo=None):
        """
        Add this to the non-redundant list of 'Primer.sequences' objects if it
        doesn't already exist.
        
        Returns the object
        """
        p = cls.sequences.get(sequence, None)
        if (p == None):
            p = cls(sequence, position, template_length, strand, checks=checks, name=name, o_oligo=o_oligo)
            cls.sequences[sequence] = p
        
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
        
    def check_length(self, length=(17,28)):
        # Check length of primer
        # primer length should be 17-28 nt long
        if (length != None): # 'None' means skip this calculation
            length_passed = length[0] <= len(self.sequence) <= length[1]
        else:
            length_passed = None
        return length_passed
    
    def check_last5gc_count(self, last5gc_count=(1,3)):
        # Does the last 5 nt of the sequence have 1-3 C/G bases?
        # should avoid runs of 3-or-more Cs and Gs in 3' end
        if (last5gc_count != None): # 'None' means skip this calculation
            C_count = self.sequence[-5:].count('C')
            G_count = self.sequence[-5:].count('G')
            last5gc_count_passed = last5gc_count[0] <= C_count+G_count <= last5gc_count[1]
        else:
            last5gc_count_passed = None
        
        return last5gc_count_passed
    
    def check_gc_clamp_length(self, gc_clamp_length=(1,2)):
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
    
    def check_gc_freq(self, gc=(0.25,0.75)):
        # %GC should be between 40-60
        if (gc != None): # 'None' means skip this calculation
            gc_passed = gc[0] <= self.gc <= gc[1]
        else:
            gc_passed = None
            
        return gc_passed
    
    def check_max_run_length(self, max_run_length=4):
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
    
    def check_3prime_complementation_length(self, max_3prime_complementation_length=3):
        # 3' ends of primers should NOT be complementary
        # As even a little dimerization will inhibit target annealing
        if (max_3prime_complementation_length != None):
            max_3prime_complementation_length_passed = self.get_3prime_complementation_length(self.sequence, self.sequence, max_3prime_search_length=max_3prime_complementation_length+1) <= max_3prime_complementation_length
        else:
            max_3prime_complementation_length_passed = None
            
        return max_3prime_complementation_length_passed
    
    def check_min_delta_g(self, min_delta_g=-5.0, folder='/dev/shm/addtag', o_oligo=None):
        o1 = None
        o2 = None
        
        # min delta-G should be -3.0
        if ((min_delta_g != None) and (o_oligo != None)): # 'None' means skip this calculation
            # Calculate hairpin delta-G
            #o1 = Structure.new_calculate_simple(folder, self.sequence)
            o1 = o_oligo.find_structures(folder, self.sequence)
            
            # Calculate homodimer delta-G
            #o2 = Structure.new_calculate_simple(folder, self.sequence, self.sequence)
            o2 = o_oligo.find_structures(folder, self.sequence, self.sequence)
            
            min_delta_g_passed = min_delta_g <= min(min(o1).delta_G, min(o2).delta_G)
        else:
            min_delta_g_passed = None
        
        return min_delta_g_passed, o1, o2
    
    def check_tm(self, tm=(52,65), folder='/dev/shm/addtag', o_oligo=None):
        o3 = None
        
        # Tm (against reverse_complement) should be between 50-70 or 55-65
        if ((tm != None) and (o_oligo != None)): # 'None' means skip this calculation
            # Calculate reverse-complement delta-G and Tm
            #o3 = Structure.new_calculate_simple(folder, self.sequence, rc(self.sequence))
            o3 = o_oligo.find_structures(folder, self.sequence, rc(self.sequence))
        
            tm_passed = tm[0] <= min(o3).melting_temperature <= tm[1]
        else:
            tm_passed = None
        
        return tm_passed, o3
    
    def summarize(self, results):
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
    pairs = {}
    
    def __init__(self, forward_primer_object, reverse_primer_object, *args, checks=True, intervening=0, o_oligo=None, **kwargs):
        self.forward_primer = forward_primer_object
        self.reverse_primer = reverse_primer_object
        
        #the_pair1 = self.pairs.setdefault((self.forward_primer.sequence, self.reverse_primer.sequence), self)
        # Only do expensive calculations if this is a pair that hasn't been added yet
        #if (self == the_pair):
        self.intervening = intervening
        
        if checks:
            self.checks, self.o_heterodimer = self.check_pair_attributes(o_oligo=o_oligo)
        else:
            self.checks = [None]*4
            self.o_heterodimer = None
        
        self.weight = self.get_weight(**kwargs)
    
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
    
    def check_pair_attributes(self, thorough=False, o_oligo=None):
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
        
            checks[0] = self.check_amplicon_size()
            if checks[0]:
                checks[1] = self.check_max_tm_difference()
                if checks[1]:
                    checks[2] = self.check_3prime_complementation_length()
                    if checks[2]:
                        checks[3], o4 = self.check_min_delta_g(o_oligo=o_oligo)
        
        return checks, o4
    
    def check_amplicon_size(self, amplicon_size_range=(300,700)):
        if (amplicon_size_range != None):
            amplicon_size_passed = amplicon_size_range[0] <= self.get_amplicon_size() <= amplicon_size_range[1]
        else:
            amplicon_size_passed = None
        return amplicon_size_passed
    
    def check_max_tm_difference(self, max_tm_difference=3.0):
        # The difference in Tms should be as small as possible
        if (max_tm_difference != None):
            tm1 = min(self.forward_primer.o_reverse_complement).melting_temperature
            tm2 = min(self.reverse_primer.o_reverse_complement).melting_temperature
            
            max_tm_difference_passed = abs(tm1-tm2) <= max_tm_difference
        else:
            max_tm_difference_passed = None
        return max_tm_difference_passed
    
    def check_3prime_complementation_length(self, max_3prime_complementation_length=3):
        # 3' ends of primers should NOT be complementary
        # As even a little dimerization will inhibit target annealing
        if (max_3prime_complementation_length != None):
            max_3prime_complementation_length_passed = self.forward_primer.get_3prime_complementation_length(self.forward_primer.sequence, self.reverse_primer.sequence, max_3prime_search_length=max_3prime_complementation_length+1) <= max_3prime_complementation_length
        else:
            max_3prime_complementation_length_passed = None
            
        return max_3prime_complementation_length_passed
    
    def check_min_delta_g(self, min_delta_g=-5.0, folder='/dev/shm/addtag', o_oligo=None):
        o4 = None
        
        # min delta-G should be -3.0
        if ((min_delta_g != None) and (o_oligo != None)): # 'None' means skip this calculation
            # Calculate heterodimer delta-G
            #o4 = Structure.new_calculate_simple(folder, self.forward_primer.sequence, self.reverse_primer.sequence)
            o4 = o_oligo.find_structures(folder, self.forward_primer.sequence, self.reverse_primer.sequence)
            
            min_delta_g_passed = min_delta_g <= min(o4).delta_G
        else:
            min_delta_g_passed = None
        
        return min_delta_g_passed, o4
    
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
        labs = ['seq', 'amplicon_size', 'Tm', 'GC', 'min(dG)']
        vals = [
            (self.forward_primer.sequence, self.reverse_primer.sequence),
            self.get_amplicon_size(),
            tuple(round(x, 2) for x in self.get_tms()),
            self.get_gcs(),
            round(self.get_min_delta_G(), 2)
        ]
            
        return self.__class__.__name__ + '(' + ', '.join('='.join(map(str, x)) for x in zip(labs, vals)) + ')'

