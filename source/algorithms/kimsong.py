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

# Keras must be set to use the 'theano' backend
# DeepCpf1 will work on either Python 2.7 or 3.5/63.


# Import standard packages
import sys
import os
import math
import fractions
import random
import logging
import importlib.util

# Import non-standard packages
import regex
import subprocess

logger = logging.getLogger(__name__)

# Treat modules in PACKAGE_PARENT as in working directory
if (__name__ == "__main__"):
    # Relative path for package to import
    PACKAGE_PARENT = '..'
    # Obtain path of currently-running file
    SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
    # Convert to absolute path, and add to the PYTHONPATH
    sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))
    
    from algorithm import BatchedSingleSequenceAlgorithm
    from nucleotides import rc, disambiguate_iupac, random_sequence
    from utils import which
else:
    from .algorithm import BatchedSingleSequenceAlgorithm
    from ..nucleotides import rc, disambiguate_iupac
    from ..utils import which

class DeepCpf1(BatchedSingleSequenceAlgorithm):
    
    logger = logger.getChild(__qualname__)
    
    METHOD_EXCLUSIVE = 0
    METHOD_UNION = 1
    CHROMATIN_IGNORED = 0
    CHROMATIN_INACCESSIBLE = 1
    CHROMATIN_ACCESSIBLE = 2
        
    def __init__(self):
        super().__init__(
            name="DeepCpf1",
            authors=['Kim, Hui K.', 'Song, Myungjae', 'Lee, Jinu', 'Menon, A. Vipin', 'Jung, Soobin', 'Kang, Young-Mook', 'Choi, Jae W.', 'Woo, Euijeon', 'Koh, Hyun C.', 'Nam, Jin-Wu', 'Kim, Hyongbum'],
            title='In vivo high-throughput profiling of CRISPR-Cpf1 activity',
            journal='Nature Methods',
            issuing='14:153-159',
            year=2016,
            doi='https://doi.org/10.1038/nmeth.4104',
            #citation="Kim, Song, et al. In vivo high-throughput profiling of CRISPR-Cpf1 activity. Nature Methods 14, 153-159 (2017).",
            off_target=False,
            on_target=True,
            prefilter=False,
            postfilter=False,
            minimum=1.0,
            maximum=150.0,
            default=None,
            weight_str='DeepCpf1:60+1.08' # TODO: Double-check this weight
        )
        
        self.loaded = False
        
        self.seq_model = None
        self.ca_model = None
    
    def is_available(self):
        """
        Determines if the prerequisites for DeepCpf1 have been met.
        :return: True or False
        """
        # TODO: Finish this
        #       Add 'theano'
        #       Have it look for the 'kimsing_ca_scores.h5'
        #       Have it look for the 'kimsong_seq_scores.h5'
        spec = importlib.util.find_spec('keras')
        if spec:
            return True
        else:
            return False
    
    def load(self):
        global zeros
        global ones
        
        from numpy import zeros, ones
        
        global backend
        global Model
        global Input
        global Multiply
        global Dense
        global Dropout
        global Flatten
        global Convolution1D
        global AveragePooling1D
        
        # Force Keras to use 'theano' backend
        # by setting environmental variable before import
        os.environ['KERAS_BACKEND'] = 'theano'
        
        from keras import backend as backend
        from keras.models import Model
        from keras.layers import Input
        from keras.layers.merge import Multiply
        from keras.layers.core import Dense, Dropout, Flatten
        from keras.layers.convolutional import Convolution1D, AveragePooling1D
        
        # Raise exception if Theano is missing
        if backend.backend() != 'theano':
            raise Exception("DeepCpf1 requires 'Keras' backend to be set to 'theano'.")
        
        self.seq_model = self.build_seq_model()
        self.ca_model = self.build_ca_model()
    
    # def weight(self, x):
    #     """Penalize any score less than 60""" # Is this appropriate? Should it be 40?
    #     return 1.0/(1+1.08**(60-x))
    
    def build_seq_model(self):
        self.logger.info("Building models")
        Seq_deepCpf1_Input_SEQ = Input(shape=(34,4))
        Seq_deepCpf1_C1 = Convolution1D(80, 5, activation='relu')(Seq_deepCpf1_Input_SEQ)
        Seq_deepCpf1_P1 = AveragePooling1D(2)(Seq_deepCpf1_C1)
        Seq_deepCpf1_F = Flatten()(Seq_deepCpf1_P1)
        Seq_deepCpf1_DO1= Dropout(0.3)(Seq_deepCpf1_F)
        Seq_deepCpf1_D1 = Dense(80, activation='relu')(Seq_deepCpf1_DO1)
        Seq_deepCpf1_DO2= Dropout(0.3)(Seq_deepCpf1_D1)
        Seq_deepCpf1_D2 = Dense(40, activation='relu')(Seq_deepCpf1_DO2)
        Seq_deepCpf1_DO3= Dropout(0.3)(Seq_deepCpf1_D2)
        Seq_deepCpf1_D3 = Dense(40, activation='relu')(Seq_deepCpf1_DO3)
        Seq_deepCpf1_DO4= Dropout(0.3)(Seq_deepCpf1_D3)
        Seq_deepCpf1_Output = Dense(1, activation='linear')(Seq_deepCpf1_DO4)
        Seq_deepCpf1 = Model(inputs=[Seq_deepCpf1_Input_SEQ], outputs=[Seq_deepCpf1_Output])
        
        # Load scores from the data file into the model
        try:
            self.logger.info("Loading scores for the SEQ models")
            Seq_deepCpf1.load_weights(os.path.join(os.path.dirname(__file__), 'kimsong_seq_scores.h5')) # Renamed 'Seq_deepCpf1_weights.h5' file
        except FileNotFoundError:
            raise Exception("Could not find DeepCpf1 data files containing scores")
        
        return Seq_deepCpf1
        
    def build_ca_model(self):
        DeepCpf1_Input_SEQ = Input(shape=(34,4))
        DeepCpf1_C1 = Convolution1D(80, 5, activation='relu')(DeepCpf1_Input_SEQ)
        DeepCpf1_P1 = AveragePooling1D(2)(DeepCpf1_C1)
        DeepCpf1_F = Flatten()(DeepCpf1_P1)
        DeepCpf1_DO1= Dropout(0.3)(DeepCpf1_F)
        DeepCpf1_D1 = Dense(80, activation='relu')(DeepCpf1_DO1)
        DeepCpf1_DO2= Dropout(0.3)(DeepCpf1_D1)
        DeepCpf1_D2 = Dense(40, activation='relu')(DeepCpf1_DO2)
        DeepCpf1_DO3= Dropout(0.3)(DeepCpf1_D2)
        DeepCpf1_D3_SEQ = Dense(40, activation='relu')(DeepCpf1_DO3)
        
        DeepCpf1_Input_CA = Input(shape=(1,))
        DeepCpf1_D3_CA = Dense(40, activation='relu')(DeepCpf1_Input_CA)
        DeepCpf1_M = Multiply()([DeepCpf1_D3_SEQ, DeepCpf1_D3_CA])
        
        DeepCpf1_DO4= Dropout(0.3)(DeepCpf1_M)
        DeepCpf1_Output = Dense(1, activation='linear')(DeepCpf1_DO4)
        DeepCpf1 = Model(inputs=[DeepCpf1_Input_SEQ, DeepCpf1_Input_CA], outputs=[DeepCpf1_Output])
        
        # Load scores from the data file into the model
        try:
            self.logger.info("Loading scores for the CA model")
            DeepCpf1.load_weights(os.path.join(os.path.dirname(__file__), 'kimsong_ca_scores.h5')) # Renamed 'DeepCpf1_weights.h5' file
        except FileNotFoundError:
            raise Exception("Could not find DeepCpf1 data files containing scores")
        
        return DeepCpf1
    
    def calculate(self, batch, *args, disambiguate=False, disambiguate_samples=512, method=None, chromatin=None, **kwargs):
        
        if (method == None):
            method = self.METHOD_EXCLUSIVE
        
        if (chromatin == None):
            chromatin = self.CHROMATIN_IGNORED
        
        if not self.loaded:
            self.load()
            self.loaded = True
        
        # This should build a temporary file that gets read-in by DeepCpf1.py, DeepCpf1.py should write an output file.
        # Tab-delimited
        # INT        34-nt target (4 bp + PAM + 23 bp protospacer + 3 bp)     <EMPTY>     0 or 1
        # 1	TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT		1
        # 2	GTTATTTGAGCAATGCCACTTAATAAACATGTAA		0
        
        
        
        
        # COLUMN    I/O  DESCRIPTION
        # 0         I/O  integer that just counts up for each input/output sequence
        # 1         I/O  34 nt target sequence: US (4 nt) + PAM (TTTN) + SPACER (23 nt) + DS (3 nt)
        # 2         I/O  EMPTY
        # 3         I/O  (0 or 1) indicating whether endonuclease can cleave this site (binary chromatin information)
        #                1 means the site is hypersensitive to endonuclease (the chromatin is accessible)
        #                0 means the site is non-sensitive to endonuclease (the chromatin is inaccessible)
        # 4           O  The on-target score when ignoring chromatin accessibility state
        # 5           O  The on-target score taking chromatin accessibility into account
        
        
        
        
        # 'queries2' format: (index, calculate=True/False, sequence, default_score=0.0)
        queries2 = []
        
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
        
        # Internal equivalent of PREPROCESS
        data_n = len(queries2)
        SEQ = zeros((data_n, 34, 4), dtype=int)
        
        # q[2] contains the 34-nt built_query sequence
        
        if (method == self.METHOD_EXCLUSIVE):
            for i, q in enumerate(queries2):
                for j, qq in enumerate(q[2]):
                    if (qq == 'A'):
                        SEQ[i, j, 0] = 1
                    elif (qq == 'C'):
                        SEQ[i, j, 1] = 1
                    elif (qq == 'G'):
                        SEQ[i, j, 2] = 1
                    elif (qq == 'T'):
                        SEQ[i, j, 3] = 1
            
        elif (method == self.METHOD_UNION):
            iupac = {
                'A': ['A'],
                'C': ['C'],
                'G': ['G'],
                'T': ['T'],
                'R': ['A', 'G'],
                'Y': ['C', 'T'],
                'M': ['A', 'C'],
                'K': ['G', 'T'],
                'W': ['A', 'T'],
                'S': ['C', 'G'],
                'B': ['C', 'G', 'T'],
                'D': ['A', 'G', 'T'],
                'H': ['A', 'C', 'T'],
                'V': ['A', 'C', 'G'],
                'N': ['A', 'C', 'G', 'T'],
            }
            
            for i, q in enumerate(queries2):
                for j, qq in enumerate(q[2]):
                    if 'A' in iupac[qq]:
                        SEQ[i, j, 0] = 1
                    if 'C' in iupac[qq]:
                        SEQ[i, j, 1] = 1
                    if 'G' in iupac[qq]:
                        SEQ[i, j, 2] = 1
                    if 'T' in iupac[qq]:
                        SEQ[i, j, 3] = 1
        
        self.logger.info("Predicting on test data")
        if (chromatin == self.CHROMATIN_IGNORED):
            score_seq = self.seq_model.predict([SEQ], batch_size=50, verbose=0)
        elif (chromatin == self.CHROMATIN_INACCESSIBLE):
            CA0 = zeros((data_n, 1), dtype=int)
            #score_ca_0 = self.ca_model.predict([SEQ, CA0], batch_size=50, verbose=0) * 3
            score_seq = self.ca_model.predict([SEQ, CA0], batch_size=50, verbose=0) * 3
        elif (chromatin == self.CHROMATIN_ACCESSIBLE):
            CA1 = ones((data_n, 1), dtype=int)*100
            #score_ca_1 = self.ca_model.predict([SEQ, CA1], batch_size=50, verbose=0) * 3
            score_seq = self.ca_model.predict([SEQ, CA1], batch_size=50, verbose=0) * 3
        
        scores = []
        
        if disambiguate:
            for i, q in enumerate(queries2):
                if q[1]:
                    q[3] = float(score_seq[i]) # float(score_ca_0[i]) float(score_ca_1[i])
                #else:
                #    self.logger.info("Skipping DeepCpf1 calculation: {}".format(q))
            
            q_avg_list = []
            current_i = None
            current_list = None
            
            def loop_process():
                if (current_i != None):
                    q_avg = current_list[0][:]
                    #q_avg[2] = batch[q_avg[0]]
                    q_avg[3] = sum(x[3] for x in current_list)/len(current_list)
                    q_avg_list.append(q_avg)
            
            for q in queries2:
                if (q[0] == current_i):
                    current_list.append(q)
                else:
                    loop_process()
                    
                    current_i = q[0]
                    current_list = [q]
            loop_process()
            
            for q in q_avg_list:
                #print(q)
                scores.append(q[3])
        else:
            for i, q in enumerate(queries2):
                q[3] = float(score_seq[i])
                #print(q, float(score_seq[i]), float(score_ca_0[i]), float(score_ca_1[i]))
                #print(q)
                scores.append(q[3])
        
        assert len(scores) == len(batch)
        
        return scores
    
    def build_query(self, target, pam, upstream='', downstream='', missing='N'):
        '''
        34 nt target sequence: US (4 nt) + PAM (TTTN) + SPACER (23 nt) + DS (3 nt)
        '''
        
        # us   pam  seq                     ds
        # ACAG TTTA CTGATCTCCAGATATGACCATGG GTT
        # CAGC TTTC TGATCTCCAGATATGACCATGGG TTT
        # CCAG TTTG AAGTTTGAGCCACAAACCCATGG TCA
        
        new = [missing] * 34
        
        _pam = pam[-8:]
        for i, nt in enumerate(_pam):
            new[8-len(_pam)+i] = nt
        
        _us = upstream[max(0, len(upstream)-(8-len(_pam))):]
        for i, nt in enumerate(_us):
            new[8-len(_pam)-len(_us)+i] = nt
        
        _target = target[:26]
        for i, nt in enumerate(_target):
            new[8+i] = nt
        
        _ds = downstream[:34-8-len(_target)]
        for i, nt in enumerate(_ds):
            new[8+len(_target)+i] = nt
        
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
    # V = (sequence, target, pam, upstream, downstream)
    a = ('', 'CGATGGCTAGGATCGATTGACCA', 'TTTA', 'CCGG', 'GTT')
    b = ('', 'RYMKWSACGTbDHVNACGTAAAC', 'TTTC', 'AATA', 'AAA')
    c = ('', 'ATGSCTCGGATCGATTGAGGATA', 'TTTT', 'ACCA', 'TAT')
    d = ('', 'GCGATGCGCAGCTAGGCCGGTTG', 'TTTG', 'GAGG', 'TTG')
    e = ('', 'CGAAGGCTCGGACCGATTGATGG', 'TTTA', 'CCTT', 'GGT')
    f = ('', 'CGCTGGCTAGGATCGATTGAAAT', 'TTAA', 'TTTA', 'GGC')
    g = ('', 'AAAATTAACTATAGGTAAAGCCA', 'TTTN', 'TATT', 'CCG')
    h = ('', 'AACATCAACTCTAGCTAACGCCA', 'TATT', 'CCAA', 'CCA')
    i = ('', 'AACATCAACTCTACCTAACGCCA', 'TTAT', 'CCGA', 'ACA')
    j = ('', 'GTTAGCGGTAATGTG', 'TTTN', 'GGGA', 'CCA')
    k = ('', 'CTCAAACCATTCATGGTATGTATATGTG', 'TTTG', 'TCGA', 'TTA')
    l = ('', 'GGCATGCGCCATCGCCGGACCAA', 'NNNN', 'NNNN', 'NNN')
    m = ('', 'GGCATGCGCCATCGCCGGANGGA', 'NNNN', 'NNNN', 'NNN')
    n = ('', 'GAAAATTGGCATAACCACCACAA', 'TTTC', 'ACAAAATC', 'TCATTGC')
    o = ('', 'AATGGAGTCGTGAGCGCAAGAAC', 'TTTG', 'TGAC', 'GCT') # should be 55.699310, 8.181222, 46.077484
    p = ('', 'AGCAATGCCACTTAATAAACATG', 'TTTG', 'GTTA', 'TAA') # should be 53.469837, 7.932923, 44.095497
    
    print("=== DeepCpf1 ===")
    C = DeepCpf1()
    for I, V in enumerate([a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p]):
        print(I, [V], C.calculate([V]))
        print(I, [V], C.calculate([V], disambiguate=True))
    
    for METHOD in [DeepCpf1.METHOD_EXCLUSIVE, DeepCpf1.METHOD_UNION]:
        for CHROMATIN in [DeepCpf1.CHROMATIN_IGNORED, DeepCpf1.CHROMATIN_INACCESSIBLE, DeepCpf1.CHROMATIN_ACCESSIBLE]:
            print('MEHTOD={}, CHROMATIN={}'.format(METHOD, CHROMATIN))
            Q = [g, h, i]
            print(Q, C.calculate(Q, method=METHOD, chromatin=CHROMATIN))
            print(Q, C.calculate(Q, method=METHOD, chromatin=CHROMATIN, disambiguate=True))
    
    
    
    import time
    
    start = time.time()
    seqs = []
    for I in range(100000):
        seqs.append(('', random_sequence(23), random_sequence(4), random_sequence(4), random_sequence(3)))
    scores = C.calculate(seqs, disambiguate=True)
    
    print('min={}, max={}, mean={}, time={}'.format(min(scores), max(scores), sum(scores)/len(scores), time.time()-start))
    # Scores range from 1 to 143 with a mean of 60
    
    def figure(scores, C, x_min=0, x_max=100, title=''):
        import matplotlib.pyplot as plt
        
        fig, ax1 = plt.subplots(figsize=(8, 4))
        num_bins = 50
        n, bins, patches = ax1.hist(scores, num_bins, color='orange', histtype='step', label='Histogram') # The histogram
        ax2 = ax1.twinx()
        n, bins, patches = ax2.hist(scores, num_bins, density=True, histtype='step', cumulative=True, label='Cumulative') # plot the cumulative histogram
        
        x = list(range(x_min, x_max+1))
        y = []
        for i in x:
            y.append(C.weight(i))
        ax2.plot(x, y, 'k--', linewidth=1.5, label='Weight')
        
        ax2.grid(True)
        
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper left')
        
        ax2.set_title(title)
        ax1.set_xlabel('Score')
        ax1.set_ylabel('Count')
        ax2.set_ylabel('Frequency')
        fig.tight_layout()
        plt.show()
    
    figure(scores, C, x_max=140, title='Random DeepCpf1 scores')

if (__name__ == "__main__"):
    test()
