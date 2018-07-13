#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/oligos/_unafold.py

# Import standard packages
import sys
import os
import subprocess

# Import non-standard packages
import regex

# import included AddTag-specific modules
if (__name__ == "__main__"):
    from oligo import Oligo, Primer, PrimerPair, lr_justify
else:
    from .oligo import Oligo, Primer, PrimerPair, lr_justify

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

debug = False

class UNAFold(Oligo):
    def __init__(self):
        super().__init__("UNAFold", "Markham, et al", 2008,
            citation="Markham, et al. UNAFold: Software for Nucleic Acid Folding and Hybridization. Bioinformatics: Structure, Function and Applications. Humana Press. p.3-31 (2008)."
        )
    
    def scan_sequence(self, seq, primer_size=(20,24), amplicon_size=(50,60)):
        # sliding window for F and R
        good_forward = []
        good_reverse = []
        for p_len in range(primer_size[0], primer_size[1]+1):
            if debug:
                print(seq)
            for pos in range(len(seq) - p_len+1):
                pf = seq[pos:pos+p_len]
                
                # Check pf (constraints) & (self, homodimer, rc) to see if it is good.
                pf_results, o_a, o_aa, o_ar, pf_gc = self.check_potential_primer(pf)
                
                if self.summarize(pf_results):
                    good_forward.append(Primer(sequence=pf, position=pos, strand='+', o_hairpin=o_a, o_self_dimer=o_aa, o_reverse_complement=o_ar, gc=pf_gc, checks=pf_results))
                
                if debug:
                    print(lr_justify(
                        ' '*pos + pf,
                        'fwd '+str(pf_results)
                    ))
                
                pr = rc(pf)
                # check pr (constraints) & (self, homodimer, rc) to see if it is good.
                pr_results, o_b, o_bb, o_br, pr_gc = self.check_potential_primer(pr)
                            
                if self.summarize(pr_results):
                    good_reverse.append(Primer(sequence=pr, position=pos, strand='-', o_hairpin=o_b, o_self_dimer=o_bb, o_reverse_complement=o_br, gc=pr_gc, checks=pr_results))
                
                if debug:
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
                    het_results, o_ab = self.check_potential_primer_pair(gf.sequence, gr.sequence, min(gf.o_reverse_complement).melting_temperature, min(gr.o_reverse_complement).melting_temperature)
                    if self.summarize(het_results):
                        pp.o_heterodimer = o_ab
                        pp.checks = het_results
                        good_pairs.append(pp)
        
        return good_pairs
    
    def summarize(self, results):
        return all(x for x in results if (x != None))
    
    def get_3prime_homology_length(self, seq1, seq2, max_3prime_length=5):
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

    def check_potential_primer_pair(self, seq1, seq2, tm1, tm2, thorough=False, folder='/dev/shm', max_tm_difference=2.0, max_3prime_homology_length=3, min_delta_g=-3.0):
        max_tm_difference_passed = None
        max_3prime_homology_length_passed = None
        min_delta_g_passed = None
        
        # The difference in Tms should be as small as possible
        if (max_tm_difference != None):
            max_tm_difference_passed = abs(tm1-tm2) <= max_tm_difference
        
        # 3' ends of primers should NOT be complementary
        # As even a little dimerization will inhibit target annealing
        if (max_3prime_homology_length != None):
            max_3prime_homology_length_passed = self.get_3prime_homology_length(seq1, seq2, max_3prime_homology_length+1) <= max_3prime_homology_length
            
            # seq1 5'-ACAATACGAC-3'
            #               ||||
            #       seq2 3'-GCTGTTAAG-5' <-rev- 5'-GAATTGTCG-3    --rc-> CGACAATTC
            #                                    
        
        results = [max_tm_difference_passed, max_3prime_homology_length_passed, min_delta_g_passed]
        
        o = None
        
        # if previous tests all pass
        # Calculate heterodimer delta-G
        if (thorough or self.summarize(results)):
            if (min_delta_g != None):
                o = Structure.calculate_simple(folder, seq1, seq2)
                min_delta_g_passed = min_delta_g <= min(o).delta_G
        
        results = [max_tm_difference_passed, max_3prime_homology_length_passed, min_delta_g_passed]
        
        return results, o
        

    def check_potential_primer(self, seq, thorough=False, folder='/dev/shm', length=(17,28), last5gc_count=(1,3), gc_clamp_length=(1,2), gc=(0.4,0.6), max_run_length=4, min_delta_g=-3.0, tm=(55,65)):
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
        
        if (thorough or self.summarize(results)):
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

class Structure(object):
    def __init__(self, seq1, seq2, delta_G, delta_H, delta_S, melting_temperature, sodium, magnesium, temperature, concentration):
        self.seq1 = seq1
        self.seq2 = seq2
        self.delta_G = delta_G
        self.delta_H = delta_H
        self.delta_S = delta_S
        self.melting_temperature = melting_temperature
        self.sodium = sodium
        self.magnesium = magnesium
        self.temperature = temperature
        self.concentration = concentration
    
    def __lt__(self, other):
        if (other != None):
            return ((self.delta_G, -self.melting_temperature) < (other.delta_G, -other.melting_temperature))
        else:
            return False
    # __gt__(), __le__(), __ne__(), __ge__(), __eq__()
    
    def __repr__(self):
        labs = ['dG', 'dH', 'dS', 'Tm']
        vals = [self.delta_G, self.delta_H, self.delta_S, self.melting_temperature]
            
        return 'Structure(' + ', '.join('='.join(map(str, x)) for x in zip(labs, vals)) + ')'
    
    @staticmethod
    def load_det_file(filename):
        data = {}
        with open(filename, 'r') as flo:
            for line in flo:
                # Structure 2: dG = -9.68  dH = -76.80  dS = -225.12  Tm = 42.3    
                m = regex.search(r'Structure\s+(\d+)', line)
                if m:
                    structure = int(m.group(1))
                    
                    s_data = {}
                    matches = regex.findall(r'(\w+) = ([-+]?\d*\.?\d*)', line)
                    for label, value in matches:
                        s_data[label] = float(value)
                    data[structure] = s_data
        return data
    
    @staticmethod
    def read_seq_file(seq_filename=None):
        seq = None
        if seq_filename:
            with open(seq_filename, 'r') as flo:
                seq = flo.read().rstrip()
        return seq
    
    @staticmethod
    def make_objects(det_filename, sodium, magnesium, temperature, concentration, seq1, seq2=None):
        data = Structure.load_det_file(det_filename)
        
        object_list = []
        for i, d in data.items():
            object_list.append(Structure(seq1, seq2, d['dG'], d['dH'], d['dS'], d['Tm'], sodium, magnesium, temperature, concentration))
        return object_list
    
    @classmethod
    def calculate_full(cls, folder_p, seq1, seq2, sodium=0.05, magnesium=0.0, temperature=25, concentration=0.00000025, output_basename=None):
        # Make temporary folder
        folder = os.path.join(folder_p, 'oligos')
        os.makedirs(folder, exist_ok=True)
        
        # Create sequence files for input
        with open(os.path.join(folder, 'A.seq'), 'w') as flo:
            print(seq1, file=flo)
        with open(os.path.join(folder, 'B.seq'), 'w') as flo:
            print(seq2, file=flo)
        with open(os.path.join(folder, 'rcA.seq'), 'w') as flo:
            print(rc(seq1), file=flo)
        with open(os.path.join(folder, 'rcB.seq'), 'w') as flo:
            print(rc(seq2), file=flo)
        
        # Default concetrations: 0.25 uM = 0.00025 mM = 0.00000025 M
        basic_command_list = [
            'UNAFold.pl',
            '--NA=DNA',
            '--temp='+str(temperature),
            '--sodium='+cls.float_to_str(sodium),
            '--magnesium='+cls.float_to_str(magnesium),
            '--Ct='+cls.float_to_str(concentration),
            '--max=100'
        ]
        if output_basename:
            outfile = os.path.join(folder, output_basename)
        else:
            outfile = os.devnull
        
        with open(outfile, 'w+') as flo:
            # Do hairpins
            command_list = basic_command_list + ['A.seq']
            cp = subprocess.run(command_list, shell=False, check=True, cwd=folder, stdout=flo, stderr=subprocess.STDOUT)
            command_list = basic_command_list + ['B.seq']
            cp = subprocess.run(command_list, shell=False, check=True, cwd=folder, stdout=flo, stderr=subprocess.STDOUT)
            
            # Do homodimers
            command_list = basic_command_list + ['A.seq', 'A.seq']
            cp = subprocess.run(command_list, shell=False, check=True, cwd=folder, stdout=flo, stderr=subprocess.STDOUT)
            command_list = basic_command_list + ['B.seq', 'B.seq']
            cp = subprocess.run(command_list, shell=False, check=True, cwd=folder, stdout=flo, stderr=subprocess.STDOUT)
            
            # Do heterodimer
            command_list = basic_command_list + ['A.seq', 'B.seq']
            cp = subprocess.run(command_list, shell=False, check=True, cwd=folder, stdout=flo, stderr=subprocess.STDOUT)
            
            # Do reverse complements
            command_list = basic_command_list + ['A.seq', 'rcA.seq']
            cp = subprocess.run(command_list, shell=False, check=True, cwd=folder, stdout=flo, stderr=subprocess.STDOUT)
            command_list = basic_command_list + ['B.seq', 'rcB.seq']
            cp = subprocess.run(command_list, shell=False, check=True, cwd=folder, stdout=flo, stderr=subprocess.STDOUT)
        
        # Make the objects
        #for det_filename in ['A.det', 'B.det', 'A-A.det', 'B-B.det', 'A-B.det']:
        #    pass
        a = cls.make_objects(os.path.join(folder, 'A.det'), sodium, magnesium, temperature, concentration, seq1)
        b = cls.make_objects(os.path.join(folder, 'B.det'), sodium, magnesium, temperature, concentration, seq2)
        aa = cls.make_objects(os.path.join(folder, 'A-A.det'), sodium, magnesium, temperature, concentration, seq1, seq1)
        bb = cls.make_objects(os.path.join(folder, 'B-B.det'), sodium, magnesium, temperature, concentration, seq2, seq2)
        ab = cls.make_objects(os.path.join(folder, 'A-B.det'), sodium, magnesium, temperature, concentration, seq1, seq2)
        ra = cls.make_objects(os.path.join(folder, 'A-rcA.det'), sodium, magnesium, temperature, concentration, seq1, rc(seq1))
        rb = cls.make_objects(os.path.join(folder, 'B-rcB.det'), sodium, magnesium, temperature, concentration, seq2, rc(seq2))
        
        return a, b, aa, bb, ab, ra, rb
    
    @classmethod
    def calculate_simple(cls, folder_p, seq1, seq2=None, sodium=0.05, magnesium=0.0, temperature=25, concentration=0.00000025, output_basename=None):
        """
        Writes files to temp folder 'folder_p'
        If 1 input sequence, then UNAFold is run on seq1 only
        If 2 input sequences, then UNAFold run on seq1+seq2
        """
        # Make temporary folder
        folder = os.path.join(folder_p, 'oligos')
        os.makedirs(folder, exist_ok=True)
        
        # Create sequence files for input
        with open(os.path.join(folder, 'A.seq'), 'w') as flo:
            print(seq1, file=flo)
        if seq2:
            with open(os.path.join(folder, 'B.seq'), 'w') as flo:
                print(seq2, file=flo)
        
        # Default concetrations: 0.25 uM = 0.00025 mM = 0.00000025 M
        basic_command_list = [
            'UNAFold.pl',
            '--NA=DNA',
            '--temp='+str(temperature),
            '--sodium='+cls.float_to_str(sodium),
            '--magnesium='+cls.float_to_str(magnesium),
            '--Ct='+cls.float_to_str(concentration),
            '--max=100'
        ]
        if output_basename:
            outfile = os.path.join(folder, output_basename)
        else:
            outfile = os.devnull
        
        with open(outfile, 'w+') as flo:
            if seq2:
                # Do heterodimer
                command_list = basic_command_list + ['A.seq', 'B.seq']
                cp = subprocess.run(command_list, shell=False, check=True, cwd=folder, stdout=flo, stderr=subprocess.STDOUT)
            else:
                # Do hairpins
                command_list = basic_command_list + ['A.seq']
                cp = subprocess.run(command_list, shell=False, check=True, cwd=folder, stdout=flo, stderr=subprocess.STDOUT)
        
        # Make the objects
        if seq2:
            return cls.make_objects(os.path.join(folder, 'A-B.det'), sodium, magnesium, temperature, concentration, seq1, seq2)
        else:
            return cls.make_objects(os.path.join(folder, 'A.det'), sodium, magnesium, temperature, concentration, seq1)
    
    @staticmethod
    def float_to_str(f):
        """
        Retrieved from https://stackoverflow.com/questions/38847690/convert-float-to-string-without-scientific-notation-and-false-precision/38847691#38847691
        """
        float_string = repr(f)
        if 'e' in float_string:  # detect scientific notation
            digits, exp = float_string.split('e')
            digits = digits.replace('.', '').replace('-', '')
            exp = int(exp)
            zero_padding = '0' * (abs(int(exp)) - 1)  # minus 1 for decimal point in the sci notation
            sign = '-' if (f < 0) else ''
            if (exp > 0):
                float_string = '{}{}{}.0'.format(sign, digits, zero_padding)
            else:
                float_string = '{}0.{}{}'.format(sign, zero_padding, digits)
        return float_string

def test():
    """Code to test the classes and functions in 'source/oligos/_unafold.py'"""
    
    C = UNAFold()
    print("===", C.name, "===")
    seq = 'TTCGTGTAGGATCACACCCGTTCCAAGATGTATAATCAGGAGACTCTTACGGTTACGAGGGACCCTCATCCAAGGACTCTAGGTGCAAAGTAACCGGTGG' # 2 pairs
    
    primer_pairs = C.scan_sequence(seq)
    for pp in primer_pairs:
        #print(pp, pp.forward_primer.strand, pp.reverse_primer.strand)
        print(pp)
    
    print(seq)
    for pp in primer_pairs:
        print(' '*pp.forward_primer.position + pp.get_formatted())

if (__name__ == "__main__"):
    test()