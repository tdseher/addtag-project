#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/oligos/unafold.py

# Import standard packages
import sys
import math
import os
import subprocess
import inspect

# import non-standard package
import regex

#################################
#s = inspect.stack()
#print(len(s))
#for i in range(len(s)):
#    print(i, s[i][1])
#print(s[0][1]) # Current file
#print(s[1][1]) # file that called this file?   <frozen importlib._bootstrap>
#print(s[-1][1]) # oligos.py
#print(os.path.expanduser(__file__))
#################################

# Treat modules in PACKAGE_PARENT as in working directory
if ((__name__ == "__main__") or (os.path.basename(inspect.stack()[-1][1]) == 'oligos.py')): # or __name__ == 'unafold'):
    # Relative path for package to import
    PACKAGE_PARENT = '..'
    # Obtain path of currently-running file
    SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
    # Convert to absolute path, and add to the PYTHONPATH
    sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))
    
    from nucleotides import rc
else:
    from ..nucleotides import rc

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
    """Code to test the functions and classes"""
    
    a, b, aa, bb, ab, ra, rb = Structure.calculate_full('/home/thaddeus/primertest', 'GAGGCTGTTGAGGGTACTAATG', 'GAATCGATATCGACGCGGCG')
    print(sorted(a)) # [Structure(dG=0.3, dH=-11.1, dS=-38.24, Tm=17.1), Structure(dG=0.38, dH=-11.4, dS=-39.51, Tm=15.4), Structure(dG=0.58, dH=-18.3, dS=-63.33, Tm=15.8), Structure(dG=1.23, dH=-20.8, dS=-73.88, Tm=8.4)]
    print(sorted(b)) # [Structure(dG=-1.53, dH=-25.5, dS=-80.39, Tm=44.0)]
    print(sorted(aa)) # [Structure(dG=-3.51, dH=-29.1, dS=-85.82, Tm=-22.3), Structure(dG=-3.09, dH=-42.0, dS=-130.5, Tm=-11.8)]
    print(sorted(bb)) # [Structure(dG=-9.86, dH=-60.2, dS=-168.83, Tm=29.3), Structure(dG=-9.68, dH=-76.8, dS=-225.12, Tm=27.6)]
    print(sorted(ab)) # [Structure(dG=-2.28, dH=-29.9, dS=-92.65, Tm=-35.1), Structure(dG=-1.69, dH=-22.7, dS=-70.46, Tm=-53.7), Structure(dG=-1.43, dH=-50.9, dS=-165.93, Tm=-17.2)]
    print(sorted(ra)) # [Structure(dG=-24.87, dH=-169.3, dS=-484.42, Tm=54.1)]
    print(sorted(rb)) # [Structure(dG=-26.52, dH=-166.5, dS=-469.51, Tm=58.2)]
    
    a = Structure.calculate_simple('/home/thaddeus/primertest', 'GAGGCTGTTGAGGGTACTAATG')
    print(a)

if (__name__ == "__main__"):
    test()
