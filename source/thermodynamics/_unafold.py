#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/thermodynamics/_unafold.py

# Import standard packages
import sys
import os
import subprocess
import math
import logging

# Import non-standard packages
import regex

logger = logging.getLogger(__name__)

# import included AddTag-specific modules
if (__name__ == "__main__"):
    from oligo import Oligo, Structure
else:
    from .oligo import Oligo, Structure

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

class UNAFold(Oligo):
    def __init__(self):
        super().__init__("UNAFold", "Markham, et al", 2008,
            citation="Markham, et al. UNAFold: Software for Nucleic Acid Folding and Hybridization. Bioinformatics: Structure, Function and Applications. Humana Press. p.3-31 (2008)."
        )
    
    # def weighted_scan(self, seq, side, *args, primer_sizes =(20, 21, 19, 22, 23, 18, 24, 25, 26), tm_range=(52, 68), **kwargs): ##### unfinished #####
    #     """
    #     New function to make list of all potential primers, but only
    #     evaluate their 'weight' if they are good enough
    #     """
    #     good_primers = []
    #     if (side in ['left', 'forward']):
    #         for p_len in primer_sizes:
    #             for pos in range(len(seq)-p_len+1):
    #                 pf = seq[pos:pos+p_len]
    #                 pf_results, o_a, o_aa, o_ar, pf_gc = self.check_potential_primer(pf, tm=tm_range)
    #                 if self.summarize(pf_results):
    #                     good_primers.append(Primer(sequence=pf, position=pos, strand='+', o_hairpin=o_a, o_self_dimer=o_aa, o_reverse_complement=o_ar, gc=pf_gc, checks=pf_results))
    #     return good_primers
    
    def relaxing_scan(self, seq, side, *args, **kwargs): ##### unfinished #####
        """
        Repeatedly scans the region for primers. If no primers are found, then
        the constraints are relaxed, and the region is scanned again.
        Tries its best to find something.
        """
        range_primer_length = (18, 26)
        range_tm = (55, 65)
        min_delta_g = -3
        range_last5gc_count =(1,3)
        range_gc_clamp_length = (1,2)
        range_gc = (0.25,0.75)
        max_run_length = 4
        
        primer_set = []
        max_retries = 10
        retries = 0
        while ((len(primer_set) == 0) and (retries < max_retries)):
            primer_set = self.scan(seq, side, primer_size=range_primer_length, tm_range=range_tm)
            range_primer_length[1] += 2
            range_tm[0] -= 1
            range_tm[1] += 1
            range_gc[0] -= 0.05
            range_gc[1] += 0.05
            min_delta_g -= 1
        return primer_set
    
    # def scan_sequence(self, seq, primer_size=(20,24), amplicon_size=(50,60)):
    #     """
    #     Finds optimal primer pairs for input sequence
    #     """
    #     # sliding window for F and R
    #     good_forward = []
    #     good_reverse = []
    #     for p_len in range(primer_size[0], primer_size[1]+1):
    #         if debug:
    #             print(seq)
    #         for pos in range(len(seq) - p_len+1):
    #             pf = seq[pos:pos+p_len]
    #             
    #             # Check pf (constraints) & (self, homodimer, rc) to see if it is good.
    #             pf_results, o_a, o_aa, o_ar, pf_gc = self.check_potential_primer(pf)
    #             
    #             if self.summarize(pf_results):
    #                 good_forward.append(Primer(sequence=pf, position=pos, strand='+', o_hairpin=o_a, o_self_dimer=o_aa, o_reverse_complement=o_ar, gc=pf_gc, checks=pf_results))
    #             
    #             if debug:
    #                 print(lr_justify(
    #                     ' '*pos + pf,
    #                     'fwd '+str(pf_results)
    #                 ))
    #             
    #             pr = rc(pf)
    #             # check pr (constraints) & (self, homodimer, rc) to see if it is good.
    #             pr_results, o_b, o_bb, o_br, pr_gc = self.check_potential_primer(pr)
    #                         
    #             if self.summarize(pr_results):
    #                 good_reverse.append(Primer(sequence=pr, position=pos, strand='-', o_hairpin=o_b, o_self_dimer=o_bb, o_reverse_complement=o_br, gc=pr_gc, checks=pr_results))
    #             
    #             if debug:
    #                 print(lr_justify(
    #                     ' '*pos + pr,
    #                     'rev '+str(pr_results)
    #                 ))
    #     
    #     # Filter pairs by amplicon size
    #     good_pairs = []
    #     for gf in good_forward:
    #         for gr in good_reverse:
    #             pp = PrimerPair(gf, gr, o_heterodimer=None, checks=None)
    #             
    #             if (amplicon_size[0] <= pp.get_amplicon_size() <= amplicon_size[1]):
    #                 het_results, o_ab = self.check_potential_primer_pair(gf.sequence, gr.sequence, min(gf.o_reverse_complement).melting_temperature, min(gr.o_reverse_complement).melting_temperature)
    #                 if self.summarize(het_results):
    #                     pp.o_heterodimer = o_ab
    #                     pp.checks = het_results
    #                     good_pairs.append(pp)
    #     
    #     return good_pairs
    
    def find_structures(self, *args, **kwargs):
        """
        Should return the list of structures with delta-G values.
        """
        
        return UNAFoldStructure.new_calculate_simple(*args, **kwargs)

class UNAFoldStructure(Structure):
    
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
        data = UNAFoldStructure.load_det_file(det_filename)
        
        object_list = []
        for i, d in data.items():
            object_list.append(UNAFoldStructure(seq1, seq2, d['dG'], d['dH'], d['dS'], d['Tm'], sodium, magnesium, temperature, concentration))
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
    def new_calculate_simple(cls, folder_p, seq1, seq2=None, sodium=0.05, magnesium=0.0, temperature=25, concentration=0.00000025):
        # Make temporary folder
        folder = os.path.join(folder_p, 'oligos')
        os.makedirs(folder, exist_ok=True)
        
        # Run simplified UNAFold scripts to get thermodynamic properties
        deltaGs, deltaHs, deltaSs, Tms = cls.nnn_unafold(folder, seq1, seq2, sodium, magnesium, temperature, concentration)
        
        # Build a list of Structure objects
        object_list = []
        for i in range(len(deltaGs)):
            object_list.append(Structure(seq1, seq2, deltaGs[i], deltaHs[i], deltaSs[i], Tms[i], sodium, magnesium, temperature, concentration))
        
        # Return the list of Structure objects
        return object_list
    
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
                try:
                    command_list = basic_command_list + ['A.seq']
                    cp = subprocess.run(command_list, shell=False, check=True, cwd=folder, stdout=flo, stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError: # This error happens sometimes, i.e. when A.seq contains 'TTCTCCACTTCCATCACC'
                    # UNAFold.pl crashes, so we just write a crappy 'A.det' file
                    print("UNAFold.pl CRASH")
                    with open(os.path.join(folder, 'A.det'), 'w') as crash_flo:
                        print('Structure 1: dG = -999.0  dH = -999.0  dS = -999.0  Tm = 99.0', file=crash_flo)
        
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
    
    @staticmethod
    def parse_ct_file(filename):
        """
        Check each pair to determine if it's a homo- or heterodimer
        """
        homodimer_list = []
        with open(filename, 'r') as flo:
            structures = []
            structure = []
            for line in flo:
                sline = line.rstrip().split("\t")
                if (sline[1].startswith("dG")):
                    if (len(structure) > 0):
                        structures.append(structure)
                    structure = []
                    molecule = []
                else:
                    base, previous_index, next_index = sline[1], sline[2], sline[3]
                    molecule.append(base)
                    if (next_index == '0'):
                        structure.append(molecule)
                        molecule = []
            if (len(structure) > 0):
                structures.append(structure)
            
            for structure in structures:
                if (len(structure) == 2):
                    molecule0_iter = iter(structure[0])
                    molecule1_iter = iter(structure[1])
                    for nt1 in molecule1_iter:
                        nt0 = next(molecule0_iter)
                        if (nt0 != nt1):
                            homodimer_list.append(False)
                            break
                    else:     
                        homodimer_list.append(True)
                else: # Doing nothing would replicate what UNAFold does
                    homodimer_list.append(None)
        return homodimer_list
    
    @classmethod
    def nnn_unafold(cls, folder, seq1, seq2=None, sodium=0.05, magnesium=0.0, temperature=25, concentration=0.00000025):
        """
        Calculate deltaG, deltaH, deltaS, Tm.
        Faster than running 'UNAFold.pl'.
        
        Accepts 1 or 2 input sequences. Automatically runs either:
         * Hairpin     (1 input sequence: A=seq1, UNAFold run on A)
         * Homodimer   (2 identical input sequences: A=seq1=seq2, UNAFold run on A & A)
         * Heterodimer (2 input sequences: A=seq1 B=seq2, UNAFold run on A & B)
         
        Writes '*.det' file to temp 'folder'
        
        Returns four lists:
          ([deltaG, ...], [deltaH, ...], [deltaS, ...], [Tm, ...])
        """
        
        # Create sequence files for input
        with open(os.path.join(folder, 'A.seq'), 'w') as flo:
            print(seq1, file=flo)
        if seq2:
            if (seq2 != seq1):
                with open(os.path.join(folder, 'B.seq'), 'w') as flo:
                    print(seq2, file=flo)
        
        # 1 sequence
        # @command = ('hybrid-ss', @rules, @rules2, @rules3, '--tracebacks'=> $max);
        # Only used when --model=PG
        #hybrid-ss --NA DNA --tmin 25 --tmax 25 --sodium=0.05 --magnesium=0.0 --suffix DAT A.seq # <-- not actually tested yet
        
        #output_basename=None # add to function arguments
        #if output_basename:
        #    outfile = os.path.join(folder, output_basename)
        #else:
        #    outfile = os.devnull
        
        parameters_1 = [
            '--NA=DNA',
            '--tmin='+str(temperature),
            '--tmax='+str(temperature),
            '--tinc=1',
            '--sodium='+UNAFoldStructure.float_to_str(sodium),
            '--magnesium='+UNAFoldStructure.float_to_str(magnesium),
            '--maxloop='+str(30),
            '--mfold='+','.join(map(str, [5, -1, 100])),
        ]
        if (seq2 == None):
            prefix = 'A'
            # 1 sequence
            # Used by default, or when --model=EM
            # Creates files: A.ann, A.ct, A.dG, A.plot, A.run
            # Command: hybrid-ss-min --NA=DNA --tmin=25 --tmax=25 --tinc=1 --sodium=0.05 --magnesium=0.0 --maxloop=30 --mfold=5,-1,100 A.seq
            command_list = ['hybrid-ss-min'] + parameters_1 + ['A.seq']
        elif (seq1 == seq2):
            prefix = 'A-A'
            command_list = ['hybrid-min'] + parameters_1 + ['A.seq', 'A.seq']
        else:
            prefix = 'A-B'
            command_list = ['hybrid-min'] + parameters_1 + ['A.seq', 'B.seq']
        
        try:
            with open(os.devnull, 'w+') as flo: # Prevent printing to STDOUT
                cp = subprocess.run(command_list, shell=False, check=True, cwd=folder, stdout=flo, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError:
            # Some sequences, such as 'GAGAAGGAGAAGGAGAAG' paired with itself, will cause a segmentation fault
            return [math.inf], [math.inf], [math.nan], [math.nan]
        
        # If model=EM
        # Apparently, this is not needed
        # Creates file: A.h-num
        # Command: h-num.pl A
    #    with open(outfile, 'w+') as flo:
    #        command_list = ['h-num.pl', 'A']
    #        cp = subprocess.run(command_list, shell=False, check=True, cwd=folder, stdout=flo, stderr=subprocess.STDOUT)
        
        
        # Apparently, this is not needed
        # Creates file: A.ss-count
        # Command: ss-count.pl A.ct > A.ss-count
    #    with open('A.ss-count', 'w') as flo:
    #        command_list = ['ss-count.pl', 'A.ct']
    #        cp = subprocess.run(command_list, shell=False, check=True, cwd=folder, stdout=flo, stderr=None)
        
        # This calculates delta-H
        # if there is no --suffix option used '--suffix=DHD' overrides '--NA=DNA --sodium=0.05 --magnesium=0.0 --temperature=25'
        # Command: ct-energy --suffix=DHD A.ct > A.deltaH
        command_list = ['ct-energy', '--suffix=DHD', prefix+'.ct']
        cp = subprocess.run(command_list, shell=False, check=True, cwd=folder, stdout=subprocess.PIPE, stderr=None)
        deltaH_list = []
        for line in cp.stdout.decode().splitlines():
            deltaH_list.append(float(line))
        
        
        # This gives more precise delta-G calculations
        # Command: ct-energy --NA=DNA --sodium=0.05 --magnesium=0.0 --temperature=25 A.ct > A.deltaG
        parameters_2 = [
            '--NA=DNA',
            '--temperature='+str(temperature),
            '--sodium='+UNAFoldStructure.float_to_str(sodium),
            '--magnesium='+UNAFoldStructure.float_to_str(magnesium),
        ]
        command_list = ['ct-energy'] + parameters_2 + [prefix+'.ct']
        cp = subprocess.run(command_list, shell=False, check=True, cwd=folder, stdout=subprocess.PIPE, stderr=None)
        deltaG_list = []
        for line in cp.stdout.decode().splitlines():
            deltaG_list.append(float(line))
        
        # This gives the complicated information within the A.det file...
        #ct-energy --NA=DNA --sodium=0.05 --magnesium=0.0 --temperature=25 --verbose A.ct | ct-energy-det.pl --mode text
        
        number_structures = len(deltaH_list)
        homodimer_list = cls.parse_ct_file(os.path.join(folder, prefix+'.ct')) # Check each pair to determine if it's a homo- or heterodimer
        #deltaH_list = [-32.4, -49.8]
        #deltaG_list = [-3.43067, -2.56959]
        deltaS_list = []
        Tm_list = []
        factor_list = []
        
        #temperature = 25 # argument
        #concentration = 0.00000025 # argument
        R = 0.0019872 # Constant
        
        try:
            for i in range(number_structures):
                #if (len(homodimer_list) > i):
                #    homodimer = homodimer_list[i]
                #else:
                #    homodimer = None
                homodimer = homodimer_list[i]
                deltaH = deltaH_list[i]
                deltaG = deltaG_list[i]
                
                deltaS = 1000.0 * (deltaH - deltaG) / (273.15 + temperature)
                if (homodimer != None):
                    if (homodimer == True):
                        factor = 1
                    else:
                        factor = 4
                    Tm = 1000.0 * deltaH/(deltaS + 1000.0 * R * math.log(concentration/factor)) - 273.15 # Natural log (base e)
                else:
                    Tm = 1000.0 * deltaH/deltaS - 273.15
                    factor = None
                
                deltaS_list.append(deltaS)
                Tm_list.append(Tm)
                factor_list.append(factor)
        except IndexError:
            # Problematic sequences, such as 'TTCTCCACTTCCATCACC' will cause an error
            # So we artificially give them unreasonable results so they will be discarded
            for i in range(number_structures):
                deltaS_list.append(math.nan)
                Tm_list.append(math.nan)
                factor_list.append(None)
        
        # Write the det file
        with open(os.path.join(folder, prefix+'.det'), 'w') as flo:
            for i in range(number_structures):
                print('Structure {}: dG = {}  dH = {}  dS = {}  Tm = {}'.format(i+1, deltaG_list[i], deltaH_list[i], deltaS_list[i], Tm_list[i]), file=flo)
        
            
        
        # Create Probability dot plot 
        #my @command = ('hybrid-plot-ng', '--temperature' => $temp)
        #system(@command)
        
        # Create Energy dot plot
        #my @command = ('boxplot_ng', '-d', -c => 4);
        #system(@command)
        
        # Create structure plots
        #system($sirgraph, @flags, -ss => "${prefix}_$fold")
        #system($sirgraph, @flags, -p => "${prefix}_$fold")
        #system($sirgraph, @flags, $img, "${prefix}_$fold")
        #system('ps2pdfwr', "${prefix}_$fold.ps")
        
        # ????
        #system('ct2rnaml', $prefix)
        
        # These are all of the float data type
        return deltaG_list, deltaH_list, deltaS_list, Tm_list

def test():
    """Code to test the classes and functions in 'source/oligos/_unafold.py'"""
    
    C = UNAFold()
    print("===", C.name, "===")
    
    a = 'GAAATCGCTTAGCGCGAACTCAGACCAT'
    b = 'CCTAGCTATTTAATAAATC'
    c = 'TTCTCCACTTCCATCACCGT'
    
    print('Hairpin: {}'.format(repr(a)))
    for s in C.find_structures('.', a, None):
        print('', s)
    print('Homodimer: {} {}'.format(repr(a), repr(a)))
    for s in C.find_structures('.', a, a):
        print('', s)
    print('Heterodimer: {} {}'.format(repr(a), repr(b)))
    for s in C.find_structures('.', a, b):
        print('', s)
    print('Reverse-complements: {} {}'.format(repr(a), repr(rc(a))))
    for s in C.find_structures('.', a, rc(a)):
        print('', s)

def old_test():
    C = UNAFold()
    
    # Test out the Structure methods
    print("===", C.name, "===")
    
    print('single sequence/hairpin:')
    print(UNAFoldStructure.nnn_unafold('.', 'GAAATCGCTTAGCGCGAACTCAGACCAT'))
    print('homodimer:')
    print(UNAFoldStructure.nnn_unafold('.', 'GAAATCGCTTAGCGCGAACTCAGACCAT', 'GAAATCGCTTAGCGCGAACTCAGACCAT'))
    print('heterodimer:')
    print(UNAFoldStructure.nnn_unafold('.', 'GAAATCGCTTAGCGCGAACTCAGACCAT', 'CCTAGCTATTTAATAAATC'))
    
    # Try out the problem sequence
    print('problem single:')
    print(UNAFoldStructure.nnn_unafold('.', 'TTCTCCACTTCCATCACCGT'))
    print('problem homodimer:')
    print(UNAFoldStructure.nnn_unafold('.', 'TTCTCCACTTCCATCACCGT', 'TTCTCCACTTCCATCACCGT'))
    
#    p_results, o_a, o_aa, o_ar, p_gc = C.check_potential_primer('TTCTCCACTTCCATCACC') # problematic primer sequence
#    print(p_results)
#    print(o_a)
#    print(o_aa)
#    print(o_ar)
#    print(p_gc)
    
#    print("===", C.name, "===")
#    seq = 'TTCGTGTAGGATCACACCCGTTCCAAGATGTATAATCAGGAGACTCTTACGGTTACGAGGGACCCTCATCCAAGGACTCTAGGTGCAAAGTAACCGGTGG' # 2 pairs
#    
#    left_primers = C.scan(seq, 'left', primer_size=(18,22))
#    print('left_primer =', left_primers)
#    right_primers = C.scan(seq, 'right', primer_size=(18,22))
#    print('right_primer =', right_primers)
#    primer_pairs = C.pair(left_primers, right_primers, amplicon_size=(50,60))
#    print('primer_pairs =', primer_pairs)
    ##############
    #primer_pairs = C.scan_sequence(seq)
    #for pp in primer_pairs:
    #    #print(pp, pp.forward_primer.strand, pp.reverse_primer.strand)
    #    print(pp)
    #
    #print(seq)
    #for pp in primer_pairs:
    #    print(' '*pp.forward_primer.position + pp.get_formatted())

if (__name__ == "__main__"):
    test()