#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/kimsong.py


# Algorithm subclass
# Name of program: CINDEL, DeepABE/DeepCpf1
# Program URL: http://big.hanyang.ac.kr/cindel/
# Paper URL: https://www.nature.com/articles/nmeth.4104
# Source URL: https://github.com/MyungjaeSong/Paired-Library
# Citation: Published: 19 December 2016
# In vivo high-throughput profiling of CRISPR-Cpf1 activity
# Hui K Kim, Myungjae Song, Jinu Lee, A Vipin Menon, Soobin Jung, Young-Mook Kang, Jae W Choi, Euijeon Woo, Hyun C Koh, Jin-Wu Nam & Hyongbum Kim
# Nature Methods volume 14, pages 153-159 (2017)
#
# CINDEL predicts the likelihood of getting a Cpf1-induced indel at the target locus
# Thus it serves as a decent proxy for an "on-target" score


# Import standard packages
import sys
import os
import math
import fractions
import random
import logging
logger = logging.getLogger(__name__)

# Import non-standard packages
import regex
import subprocess

# Treat modules in PACKAGE_PARENT as in working directory
if (__name__ == "__main__"):
    # Relative path for package to import
    PACKAGE_PARENT = '..'
    # Obtain path of currently-running file
    SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
    # Convert to absolute path, and add to the PYTHONPATH
    sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))
    
    from algorithm import SingleSequenceAlgorithm, PairedSequenceAlgorithm, BatchedSingleSequenceAlgorithm
    from nucleotides import rc, disambiguate_iupac
    from utils import which
else:
    from .algorithm import SingleSequenceAlgorithm, PairedSequenceAlgorithm, BatchedSingleSequenceAlgorithm
    from ..nucleotides import rc, disambiguate_iupac
    from ..utils import which

class Cindel(BatchedSingleSequenceAlgorithm):
    def __init__(self):
        super().__init__("CINDEL", "Kim, Song, et al", 2016,
            citation="Kim, Song, et al. In vivo high-throughput profiling of CRISPR-Cpf1 activity. Nature Methods 14, 153-159 (2017).",
            off_target=False,
            on_target=True,
            prefilter=False,
            postfilter=False,
            minimum=1.0,
            maximum=100.0,
            default=None
        )
    
    def weight(self, x):
        """Penalize any score less than 60"""
        return 1.0/(1+1.17**(50-x))
    
    def calculate(self, batch, *args, disambiguate=False, disambiguate_samples=512, **kwargs):
        # Old 'queries2' format: (index, calculate=True/False, (seq1, seq2, ...))
        # New 'queries2' format: (index, calculate=True/False, sequence, default_score=0.0)
        queries2 = []
        
        if sys.platform.startswith('win'):
            # Command line limit:   32768 characters = 2**15, /32=1024
            batch_size = 1000
        else:
            # Command line limit: 2097152 characters = 2**21, /32=65536
            batch_size = 10000
        
        # unpack the input sequences
        for i, query in enumerate(batch):
            sequence, target, pam, upstream, downstream = query
            built_query = self.build_query(target, pam, upstream, downstream)
            built_query = built_query.upper()
            
            # Check if a non-cannonical nt is found            
            #m = regex.search('[^ACGT]', built_query)
            m = regex.findall('[^ACGT]', built_query)
            if (len(m) > 0):
                # If so, then we either disambiguate or we skip
                # Disambiguate
                if (disambiguate and (len(m) <= 10)):
                    disamb_query = disambiguate_iupac(built_query)
                    disamb_sample = random.sample(disamb_query, min(len(disamb_query), disambiguate_samples))
                    for ds in disamb_sample:
                        queries2.append([i, True, ds, 0.0])
                
                #if should_skip:
                else:
                    queries2.append([i, False, built_query, 0.0])
                
            else:
                queries2.append([i, True, built_query, 0.0])
        
        # Obtain path of currently-running file
        SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
        WRAPPER_PATH = os.path.join(SCRIPT_DIR, "azimuth_wrapper.py")
        # Should use args.python2_path
        
        PYTHON2 = 'python'
        if sys.platform.startswith('win'):
            if which('py.exe'):
                PYTHON2 = 'py.exe'
            else:
                ppl = which('python.exe', full=True)
                for pp in ppl:
                    if 'python27' in pp.lower():
                        PYTHON2 = pp
                        break
        #print('PYHTON2 =', PYTHON2)
        
        program_path_list = [PYTHON2, WRAPPER_PATH]
        if (PYTHON2 == 'py.exe'):
            program_path_list.insert(1, '-2')
        
        #print("queries2", queries2)
        
        batch_list = []
        batch_count = 0
        for bi, q in enumerate(queries2):
            if (batch_count % batch_size == 0):
                batch_list.append([])
            if q[1]:
                batch_list[-1].append(q[2]) # Append the sequence
                batch_count += 1
        
        batch_scores = []
        for current_batch in batch_list:
            if (len(current_batch) > 0):
            #for batch_start in range(0, len(queries2), batch_size):
            #    current_batch = queries2[batch_start:batch_start+batch_size]
            #    command_list = ['python', WRAPPER_PATH] + [ x[2] for x in current_batch if x[1] ]
                #command_list = [PYTHON2, WRAPPER_PATH] + current_batch
                command_list = program_path_list + current_batch
                
                #print('command_list =', command_list[:10], '...'+str(len(command_list)) if (len(command_list) > 10) else '')
                #with open(error_file, 'w+') as flo:
                #    cp = subprocess.run(command_list, shell=False, check=True, stdout=flo, stderr=subprocess.STDOUT)
                #print(command_list)
                cp = subprocess.run(command_list, shell=False, check=True, stdout=subprocess.PIPE)
                output = cp.stdout.decode()
                
                for line in output.splitlines()[1:]: #skip first line: line.startswith('No model file specified')
                    i_seq, i_score = line.split(" ")
                    batch_scores.append(100*float(i_score))
                
        # Deal with skipped queries
        queries2_iter = iter(queries2)
        #batch_scores_iter = iter(batch_scores)
        #for s in batch_scores_iter:
        for s in batch_scores:
            q = next(queries2_iter)
            #s = next(batch_scores_iter)
            
            while (q[1] == False):
                logger.info("Skipping Azimuth Calculation: {}".format(q))
                q = next(queries2_iter)
            
            q[3] = s
        
        # Average scores of queries from same origin
        o_scores2 = []
        q_prev = None
        q_sum = 0
        q_count = 0
        for q in queries2:
            if q_prev:
                if (q_prev[0] == q[0]):
                    q_sum += q[3]
                    q_count += 1
                else:
                    o_scores2.append(q_sum/q_count)
                    q_sum = q[3]
                    q_count = 1
            else:
                q_sum = q[3]
                q_count = 1
            q_prev = q
        
        if (q_count > 0):
            o_scores2.append(q_sum/q_count)
        
        #print("==ASSERT==")
        #print("len(o_scores2) =", len(o_scores2), "len(batch) =", len(batch))
        #for a in o_scores2:
        #    print(a)
        #for a in batch:
        #    print(a)
        assert len(o_scores2) == len(batch)
        
        # Old outline for how the data is stored
        # i  process seqs        scores         score
        # 0  True    ('a', 'a')  (10.0, 11.0)   10.5
        # 1  True    ('a',)      (8.0)           8.0
        # 2  False   ('a')       (0.0)           0.0
        
        return o_scores2
    
    def build_query(self, seq, pam, upstream='', downstream='', missing='N'):
        # us     seq                pam  ds
        # ACAG CTGATCTCCAGATATGACCA|TGG GTT
        # CAGC TGATCTCCAGATATGACCAT|GGG TTT
        # CCAG AAGTTTGAGCCACAAACCCA|TGG TCA
        
        new = [missing] * 30
        for i, nt in enumerate(pam[:6]):
            new[24+i] = nt

        for i, nt in enumerate(seq[:-25:-1]):
            new[24-1-i] = nt

        for i, nt in enumerate(upstream[::-1]):
            ind = 24-len(seq)-1-i
            if (ind < 0):
                break
            new[ind] = nt

        for i, nt in enumerate(downstream):
            ind = 24+len(pam)+i
            if (ind >= len(new)):
                break
            new[ind] = nt
        query = ''.join(new)
        return query
    
    def filter_query(self, seq, score="min"):
        """
        If the query contains an invalid character, then replace the invalid
        characters with all combinations of A, C, G, and T, and choose the
        min, max, mean, score.
        """
        m = regex.search('[^ACGT]', seq)
        if m:
            pass
        pass

def test():
    """Code to test the functions and classes"""
    
    a = ('',   'CGATGGCTAGGATCGATTGA', 'TGG', '', '')
    b = ('',   'RYMKWSACGTbDHVNACGTA', 'TGG', '', '')
    c = ('',     'ATGSCTCGGATCGATTGA', 'AGG', '', '')
    d = ('',   'GCGATGCGCAGCTAGGCCGG', 'CGG', '', '')
    e = ('',   'CGAAGGCTCGGACCGATTGA', 'GGG', '', '')
    f = ('',   'CGCTGGCTAGGATCGATTGA', 'AGG', '', '')
    g = ('',   'AAAATTAACTATAGGTAAAG', 'TGG', '', '')
    h = ('',   'AACATCAACTCTAGCTAACG', 'CGG', '', '')
    i = ('',   'AACATCAACTCTACCTAACG', 'CGG', 'CCGA', 'AACA')
    j = ('',   'GTTAGCGGTATGTATATGTG', 'TGG', 'GGGA', 'CTCA')
    k = ('', 'CTCAACATGGTATGTATATGTG', 'TGG', 'TCGA', 'TTCA')
    l = ('',   'GGCATGCGCCATCGCCGGAC', 'NNN', 'NNNN', 'NNN')
    m = ('',   'GGCATGCGCCATCGCCGGAN', 'NNN', 'NNNN', 'NNN')
    n = ('',   'GAAAATTGGCATAACCACCA', 'AGG', 'ACAAAATC', 'TCATTGC')
    
    print("=== Azimuth ===")
    C = Azimuth()
    print(C.calculate([l])) # [0.0]
    print(C.calculate([l], disambiguate=True)) # [41.94695886685016]
    print(C.calculate([m], disambiguate=True)) # [0]
    print(C.calculate([g])) # [0.0]
    print(C.calculate([h])) # [0.0]
    print(C.calculate([g, h, i])) # [0.0, 0.0, 68.10224199640001]
    print(C.calculate([g, h, i], disambiguate=True)) # [47.92785688539728, 68.47839267067128, 68.10224199640001]
    print(C.calculate([n])) # [70.237901082]

if (__name__ == "__main__"):
    test()
