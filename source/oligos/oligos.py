#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/oligos/oligos.py

# Import standard packages
import sys
import os
import random

# import non-standard package
import regex

# import AddTag-specific packages
if (__name__ == "__main__"):
    from unafold import Structure
    from primer import Primer, PrimerPair
else:
    from .unafold import Structure
    from .primer import Primer, PrimerPair

# Treat modules in PACKAGE_PARENT as in working directory
if (__name__ == "__main__"):
    # Relative path for package to import
    PACKAGE_PARENT = '..'
    # Obtain path of currently-running file
    SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
    # Convert to absolute path, and add to the PYTHONPATH
    sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))
    
    from nucleotides import rc
else:
    from ..nucleotides import rc

def random_sequence(length=100):
    return ''.join(random.choice(['A', 'C', 'G', 'T']) for k in range(length))

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

def get_3prime_homology_length(seq1, seq2, max_3prime_length=5):
    rc2 = rc(seq2)
    match_length_list = [0]
    
    for s1l in range(1, max_3prime_length+1):
        matches = regex.findall(seq1[-s1l:], rc2)
        if (len(matches) > 0):
            match_length_list.append(s1l)
        #print(s1l, matches)
    for s2l in range(1, max_3prime_length+1):
        matches = regex.findall(rc2[:s2l], seq1)
        if (len(matches) > 0):
            match_length_list.append(s2l)
        #print(s2l, matches)
    
    return max(match_length_list)

def check_potential_primer_pair(seq1, seq2, tm1, tm2, thorough=False, folder='/dev/shm', max_tm_difference=2.0, max_3prime_homology_length=3, min_delta_g=-3.0):
    max_tm_difference_passed = None
    max_3prime_homology_length_passed = None
    min_delta_g_passed = None
    
    # The difference in Tms should be as small as possible
    if (max_tm_difference != None):
        max_tm_difference_passed = abs(tm1-tm2) <= max_tm_difference
    
    # 3' ends of primers should NOT be complementary
    # As even a little dimerization will inhibit target annealing
    if (max_3prime_homology_length != None):
        max_3prime_homology_length_passed = get_3prime_homology_length(seq1, seq2, max_3prime_homology_length+1) <= max_3prime_homology_length
        
        # seq1 5'-ACAATACGAC-3'
        #               ||||
        #       seq2 3'-GCTGTTAAG-5' <-rev- 5'-GAATTGTCG-3    --rc-> CGACAATTC
        #                                    
    
    results = [max_tm_difference_passed, max_3prime_homology_length_passed, min_delta_g_passed]
    
    o = None
    
    # if previous tests all pass
    # Calculate heterodimer delta-G
    if (thorough or summarize(results)):
        if (min_delta_g != None):
            o = Structure.calculate_simple(folder, seq1, seq2)
            min_delta_g_passed = min_delta_g <= min(o).delta_G
    
    results = [max_tm_difference_passed, max_3prime_homology_length_passed, min_delta_g_passed]
    
    return results, o
    

def check_potential_primer(seq, thorough=False, folder='/dev/shm', length=(17,28), last5gc_count=(1,3), gc_clamp_length=(1,2), gc=(0.4,0.6), max_run_length=4, min_delta_g=-3.0, tm=(55,65)):
    """
    seq - ACGT sequence should be 5' to 3'
    """
    length_passed = None
    last5gc_count_passed = None
    gc_clamp_length_passed = None
    gc_passed = None
    max_run_length_passed = None
    min_delta_g_passed = None
    tm_passed = None
    
    # Check length of primer
    # primer length should be 17-28 nt long
    if (length != None):
        length_passed = length[0] <= len(seq) <= length[1]
    
    # Does the last 5 nt of the sequence have 1-3 C/G bases?
    # should avoid runs of 3-or-more Cs and Gs in 3' end
    if (last5gc_count != None):
        C_count = seq[-5:].count('C')
        G_count = seq[-5:].count('G')
        last5gc_count_passed = last5gc_count[0] <= C_count+G_count <= last5gc_count[1]
    
    # Does the 3' end of the sequence have the sequence 
    # 3' GC clamp
    if (gc_clamp_length != None):
        i = -1
        gcl = 0
        while (seq[i] in ['G', 'C']):
            i -= 1
            gcl += 1
        
        gc_clamp_length_passed = gc_clamp_length[0] <= gcl <= gc_clamp_length[1]
    
    gc_freq = None
    # %GC should be between 40-60
    if (gc != None):
        C_count = seq.count('C')
        G_count = seq.count('G')
        gc_freq = (C_count+G_count)/len(seq)
        gc_passed = gc[0] <= gc_freq <= gc[1]
    
    # Primers with long runs of a single base should generally be avoided as they can misprime
    if (max_run_length != None):
        matches = list(regex.finditer(r'(.)\1*', seq))
        max_run = 0
        for m in matches:
            max_run = max(max_run, len(m.group()))
        max_run_length_passed = max_run <= max_run_length
    
    results = [length_passed, last5gc_count_passed, gc_clamp_length_passed, gc_passed, max_run_length_passed, min_delta_g_passed, tm_passed]
    
    o1 = None
    o2 = None
    o3 = None
    
    if (thorough or summarize(results)):
        # min delta-G should be -3.0
        if (min_delta_g != None):
            # Calculate hairpin delta-G
            o1 = Structure.calculate_simple(folder, seq)
            
            # Calculate homodimer delta-G
            o2 = Structure.calculate_simple(folder, seq, seq)
            
            min_delta_g_passed = min_delta_g <= min(min(o1).delta_G, min(o2).delta_G)
        
        # Tm (against reverse_complement) should be between 50-70 or 55-65
        if (tm != None):
            # Calculate reverse-complement delta-G and Tm
            o3 = Structure.calculate_simple(folder, seq, rc(seq))
        
            tm_passed = tm[0] <= min(o3).melting_temperature <= tm[1]
        
    results = [length_passed, last5gc_count_passed, gc_clamp_length_passed, gc_passed, max_run_length_passed, min_delta_g_passed, tm_passed]
    
    return results, o1, o2, o3, gc_freq

def summarize(results):
    return all(x for x in results if (x != None))

class PrimerPair_old(object):
    def __init__(self, seqs, melting_temperatures, gcs, min_delta_G):
        self.seqs = seqs
        self.melting_temperatures = melting_temperatures
        self.gcs = gcs
        self.min_delta_G = min_delta_G
        self.locations = None
    def __repr__(self):
    
        labs = ['seq', 'Tm', 'GC', 'min(dG)']
        vals = [self.seqs, self.melting_temperatures, tuple(map(lambda x: round(x, 2), self.gcs)), self.min_delta_G]
            
        return 'PrimerPair(' + ', '.join('='.join(map(str, x)) for x in zip(labs, vals)) + ')'

def scan_sequence(seq, primer_size=(20,26), amplicon_size=(50,60)):
    # sliding window for F and R
    good_forward = []
    good_reverse = []
    for p_len in range(primer_size[0], primer_size[1]+1):
        print(seq)
        for pos in range(len(seq) - p_len+1):
            pf = seq[pos:pos+p_len]
            
            # Check pf (constraints) & (self, homodimer, rc) to see if it is good.
            pf_results, o_a, o_aa, o_ar, pf_gc = check_potential_primer(pf)
            
            if summarize(pf_results):
                good_forward.append(Primer(sequence=pf, position=pos, strand='+', o_hairpin=o_a, o_self_dimer=o_aa, o_reverse_complement=o_ar, gc=pf_gc, checks=pf_results))
            
            print(lr_justify(
                ' '*pos + pf,
                'fwd '+str(pf_results)
            ))
            
            pr = rc(pf)
            # check pr (constraints) & (self, homodimer, rc) to see if it is good.
            pr_results, o_b, o_bb, o_br, pr_gc = check_potential_primer(pr)
                        
            if summarize(pr_results):
                good_reverse.append(Primer(sequence=pr, position=pos, strand='-', o_hairpin=o_b, o_self_dimer=o_bb, o_reverse_complement=o_br, gc=pr_gc, checks=pr_results))
            
            print(lr_justify(
                ' '*pos + pr,
                'rev '+str(pr_results)
            ))
    
    # Filter pairs by amplicon size
    good_pairs = []
    for gf in good_forward:
        for gr in good_reverse:
            pp = PrimerPair(gf, gr, o_heterodimer=None, checks=None)
            
            if (amplicon_size[0] <= pp.get_amplicon_size() <= amplicon_size[1]):
                het_results, o_ab = check_potential_primer_pair(gf.sequence, gr.sequence, min(gf.o_reverse_complement).melting_temperature, min(gr.o_reverse_complement).melting_temperature)
                if summarize(het_results):
                    pp.o_heterodimer = o_ab
                    pp.checks = het_results
                    good_pairs.append(pp)
    
    return good_pairs


def scan_sequence_old(seq, primer_size=(20,26), amplicon_size=(50,60)):
    #rc_seq = rc(seq)
    
    #search_positions = (10, 90) # (start, end) positions within sequence to search for good primers
    #primer_size = (20, 26) # min, max
    #amplicon_size = (50, 60) # min, max
    
    selected_primers = []
    # sliding window for amplicons
    print(seq)
    for pos in range(len(seq) - amplicon_size[0]):
        for pf_len in range(primer_size[0], primer_size[1]):
            pf = seq[pos:pos+pf_len]
            
            # Check pf (constraints) & (self, homodimer, rc) to see if it is good.
            pf_results, o_a, o_aa, o_ar, pf_gc = check_potential_primer(pf)
            
            if summarize(pf_results):
                for amp_len in range(amplicon_size[0], amplicon_size[1]):
                    amp_seq = seq[pos:pos+amp_len]
                    
                    for pr_len in range(primer_size[0], primer_size[1]):
                        
                        pr = amp_seq[-pr_len:]
                        
                        # check pr (constraints) & (self, homodimer, rc) to see if it is good.
                        pr_results, o_b, o_bb, o_br, pr_gc = check_potential_primer(rc(pr))
                        
                        if summarize(pr_results):
                            # If both pf and pr are good, then calculate their (heterodimer)
                            het_results, o_ab = check_potential_primer_pair(pf, rc(pr), min(o_ar).melting_temperature, min(o_br).melting_temperature)
                            if summarize(het_results):
                                print(lr_justify(
                                    ' '*pos + pf + '-'*(len(amp_seq)-pf_len-pr_len) + pr + ' <------------',
                                    'Tm=('+str(min(o_ar).melting_temperature)+','+str(min(o_br).melting_temperature)+') min(dG)='+str(min(o_a+o_aa+o_b+o_bb+o_ab).delta_G)+' gc=('+str(pf_gc)+','+str(pr_gc)+')'
                                ))
                                selected_primers.append(PrimerPair((pf, rc(pr)), (min(o_ar).melting_temperature, min(o_br).melting_temperature), (pf_gc, pr_gc), min(o_a+o_aa+o_b+o_bb+o_ab).delta_G))
                            else:
                                print(lr_justify(
                                    ' '*pos + pf + '-'*(len(amp_seq)-pf_len-pr_len) + pr,
                                    'het '+str(het_results)
                                ))
                        else:
                            print(lr_justify(
                                ' '*pos + pf + '-'*(len(amp_seq)-pf_len-pr_len) + pr,
                                'rev '+str(pr_results)
                            ))
            else:
                print(lr_justify(
                    ' '*pos + pf,
                    'fwd '+str(pf_results)
                ))
    
    return selected_primers

def test():
    """Code to test the functions and classes"""
    
    # seq = 'TGGTCTGTCCCTAGCTATTTACGAGCCTGTATGTACGGAAGACACCGAGCTGGTCGCGCAGATCGAGTAAAATGCAGGGCGGATGCTCTCTAACAATACA' # 1 pair
    seq = 'TTCGTGTAGGATCACACCCGTTCCAAGATGTATAATCAGGAGACTCTTACGGTTACGAGGGACCCTCATCCAAGGACTCTAGGTGCAAAGTAACCGGTGG' # 2 pairs
    #seq = random_sequence(100)
    
    print(seq)
    selected_primers = scan_sequence(seq)
    
    for sp in selected_primers:
        fp = sp.forward_primer
        rp = sp.reverse_primer
        print(' '*fp.position + fp.sequence)
        print(' '*rp.position + rc(rp.sequence))
        print(sp)
    
    # Need to make sure each decent primer sequence is UNIQUE across the genome!
    

if (__name__ == "__main__"):
    test()
