#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/doench.py

# Import standard packages
import sys
import os
import math
import fractions
import random
import logging
import subprocess
import importlib.util

logger = logging.getLogger(__name__)

# Import non-standard packages
import regex

# Treat modules in PACKAGE_PARENT as in working directory
if (__name__ == "__main__"):
    # Relative path for package to import
    PACKAGE_PARENT = '..'
    # Obtain path of currently-running file
    SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
    # Convert to absolute path, and add to the PYTHONPATH
    sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))
    
    from algorithm import SingleSequenceAlgorithm, PairedSequenceAlgorithm, BatchedSingleSequenceAlgorithm
    from nucleotides import rc, disambiguate_iupac, random_sequence
    from utils import which
else:
    from .algorithm import SingleSequenceAlgorithm, PairedSequenceAlgorithm, BatchedSingleSequenceAlgorithm
    from ..nucleotides import rc, disambiguate_iupac
    from ..utils import which

class Doench2014(SingleSequenceAlgorithm):
    def __init__(self):
        super().__init__(
            name="Doench-2014",
            authors=['Doench, John G.', 'Hartenian, Ella', 'Graham, Daniel B.', 'Tothova, Zuzana', 'Hegde, Mudra', 'Smith, Ian', 'Sullender, Meagan', 'Ebert, Benjamin L.', 'Xavier, Ramnik J.', 'Root, David E.'],
            title='Rational design of highly active sgRNAs for CRISPR-Cas9–mediated gene inactivation',
            journal='Nature Biotechnology',
            issuing='32(12):1262-1267',
            year=2014,
            doi='https://doi.org/10.1038/nbt.3026',
            #citation="Doench, et al. Rational design of highly active sgRNAs for CRISPR-Cas9–mediated gene inactivation. Nature Biotechnology 32, 1262–1267 (2014)",
            off_target=True,
            on_target=True,
            prefilter=False,
            postfilter=False,
            minimum=1.0,
            maximum=100.0,
            default=None,
            rgn_list=('Cas9',),
            weight_str=None
        )
    
    def calculate(self, intended, *args, **kwargs):
        sequence, target, pam, upstream, downstream = intended
        
        return self.score(target, pam, upstream, downstream)
    
    def score(self, seq, pam, upstream='', downstream=''):
        """
        Function to calculate the sgRNA on-target efficacy score,
        as defined in Doench et al 2014
          seq = the gRNA sequence
          pam = the PAM sequence
          upstream = 0-10 bp immediately upstream of the gRNA
          downstream = 0-10 bp immediately downstream of the PAM
        Should return a score between 0.0 & 100.0, with higher numbers being
        better.
        
        **CANNOT HANDLE IUPAC AMBIGUITIES YET**
        """
        # Metric for "efficiency"
        #  Doench
        #   Range: 0-100. Linear regression model trained on 880 guides transfected
        #   into human MOLM13/NB4/TF1 cells (three genes) and mouse cells
        #   (six genes). Delivery: lentivirus. The Fusi score can be considered an
        #   updated version this score, as their training data overlaps a lot.
        #   See Doench et al.: http://www.nature.com/nbt/journal/v32/n12/full/nbt.3026.html
        
        # Code retrieved from:
        #  https://github.com/maximilianh/crisporWebsite/doenchScore.py
        # Calculates the sgRNA on-target efficacy score from the article
        # "Rational design of highly active sgRNAs for CRISPR-Cas9-mediated gene inactivation"
        # by J Doench et al. 2014
        # The authors' web tool is available at http://www.broadinstitute.org/rnai/public/analysis-tools/sgrna-design
        # Thanks to Cameron Mac Pherson at Pasteur Paris for fixing my original version. Maximilian Haeussler 2014
        
        # We anchor the scoring algorithm at the PAM sequence
        # Inefficient code
        new = ['-'] * 30
        for i, nt in enumerate(pam):
            new[24+i] = nt

        for i, nt in enumerate(seq[::-1]):
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
        new = ''.join(new)
        
        score = 0.59763615 # also called  'intercept'
        
        #guideSeq = seq[4:24]
        guideSeq = new[4:24]
        gcCount = guideSeq.count("G") + guideSeq.count("C")
        if gcCount <= 10:
            gcWeight = -0.2026259 # gcLow
        else:
            gcWeight = -0.1665878 # gcHigh
        score += abs(10 - gcCount) * gcWeight
        
        # ...shouldn't it be pos-1, so the scores are 0-indexed???
        for pos, model_seq, weight in OLD_SCORES:
            #if (seq[pos:pos + len(model_seq)] == model_seq):
            if (new[pos:pos + len(model_seq)] == model_seq):
                score += weight
        
        return (1.0/(1.0+math.exp(-score))) * 100

class Doench2016(PairedSequenceAlgorithm):
    def __init__(self):
        super().__init__(
            name="CFD",
            authors=['Doench, John G.', 'Fusi, Nicolo', 'Sullender, Meagan', 'Hegde, Mudra', 'Vaimberg, Emma W.', 'Donovan, Katherine F.', 'Smith, Ian', 'Tothova, Zuzana', 'Wilen, Craig', 'Orchard, Robert', 'Virgin, Herbert W.', 'Listgarten, Jennifer', 'Root, David E.'],
            title='Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9',
            journal='Nature Biotechnology',
            issuing='34(2):184-191',
            year=2016,
            doi='https://doi.org/10.1038/nbt.3437',
            #citation="Doench, Fusi, et al. Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. Nature Biotechnology 34, 184–191 (2016).",
            off_target=True,
            prefilter=False,
            postfilter=False,
            minimum=1.0,
            maximum=100.0,
            default=100.0,
            rgn_list=('Cas9'),
            weight_str='CFD:40+1.2' # Penalize any score less than 50
        )
    
    # def weight(self, x):
    #     """Severely penalize any score less than 50"""
    #     return 1.0/(1+1.2**(40-x))
    
    def calculate(self, intended, potential, *args, **kwargs):
        on_sequence, on_target, on_pam, on_upstream, on_downstream = intended
        off_sequence, off_target, off_pam, off_upstream, off_downstream = potential
        
        return self.score(on_target, off_target, on_pam)
    
    #def calculate(self, seq1, seq2, pam, *args, **kwargs):
    def score(self, seq1, seq2, pam, max_length=20):
        """
        Calculate the CFD score for a given gRNA and off-target sequence.
        
        Input oligonucleotides must be aligned with no gaps, and contain
        ONLY A, C, G, & T residues:
          seq1 = 20 bp genome sequence
          seq2 = 20 bp candidate gRNA sequence
           pam = 3 bp PAM sequence downstream of both seq1 & seq2
        
        Output:
          Cutting Frequency Determination (CFD) score
        
        The on-target score is from Doench, Fusi et al. (2016)
        doi:10.1038/nbt.3437 and measures activity at the target location.
        Higher scores give higher confidence that the guide will be active.
        Scores range from 0-100, and should be used to rank guides relative to
        each other.
        
        This score only penalizes mismatches in gRNA sequence, then modulates
        the score based on the PAM.
        
        The Cutting Frequency Determination (CFD) score is calculated by using the
        percent activity values provided in a matrix of penalties based on mismatches
        of each possible type at each position within the guide RNA sequence. 
        
        For example, if the interaction between the sgRNA and DNA has a single
        rG:dA ("rna G aligning with dna A") mismatch in position 6, then that
        interaction receives a score of 0.67. If there are two or more mismatches,
        then individual mismatch values are multiplied together. For example, 
        an rG:dA mismatch at position 7 coupled with an rC:dT mismatch at position
        10 receives a CFD score of 0.57 x 0.87 = 0.50.
        
        The on-target score represents the cleavage efficiency of Cas9.
        You can think of the score as the probability a given gRNA will be in
        top 20% of cleavage activity. Note that the scoring system is not linear,
        and only 5% of gRNAs receive a score of 60 or higher.
        
        see: http://portals.broadinstitute.org/gpp/public/software/sgrna-scoring-help
        """
        # Also called the Fusi score
        # Metric for "efficiency"
        #  Fusi
        #   Range: 0-100. Boosted Regression Tree model, trained on data produced
        #   by Doench et al (881 guides, MOLM13/NB4/TF1 cells + unpublished
        #   additional data). Delivery: lentivirus. See Fusi et al. 2015:
        #   http://biorxiv.org/content/early/2015/06/26/021568
        #   implemented in Azimuth?
        #    https://github.com/MicrosoftResearch/Azimuth
        #    https://www.microsoft.com/en-us/research/project/azimuth/

        
        # Return a score of 0 if sequences are of different length
        #if (len(seq1) != len(seq2)):
        #    return 0.0
        
        # Only calculate Doench score with cannonical nucleotides
        m1 = regex.search(r'[^ATCGatcg]', seq1)
        m2 = regex.search(r'[^ATCGatcg]', seq2)
        mp = regex.search(r'[^ATCGatcg]', pam[-2:])
        if ((m1 != None) or (m2 != None) or (mp != None)):
            # Otherwise return a score of 0
            return 0.0
        
        # Start with a perfect score
        score = 1
        
        # Ensure inputs are RNA sequences
        seq2 = seq2.replace('T','U')
        seq1 = seq1.replace('T','U')
        
        # Make algorithm function for spacers greater than 20 nt
        # by trimming off the left-most bases
        seq2 = seq2[-max_length:]
        seq1 = seq1[-max_length:]
        
        # Modify algorithm to allow for less than 20 nt gRNAs
        # loop should start at far-right (of longer sequence) and move left
        shorter, longer = sorted([seq1, seq2], key=len)
        #for i in range(-1, -len(shorter)-1, -1):
        for i in range(-len(shorter), 0):
            if (seq1[i] != seq2[i]):
                key = 'r'+seq1[i]+':d'+rc(seq2[i], kind="rna")+','+str(20+i+1)
                score *= MISMATCH_SCORES[key]
        
        # Reduce the score multiplicatively if there is a mismatch
        #for i in range(len(seq2)):
        #    if (seq1[i] != seq2[i]):
        #        key = 'r'+seq1[i]+':d'+rc(seq2[i], kind="rna")+','+str(i+1)
        #        score *= MISMATCH_SCORES[key]
        
        # Modulate the aggregate mismatch score by the PAM score
        # Exclude the most upstream residue in PAM motif (usually 'N'),
        # as it does not significantly affect score calculation
        score *= PAM_SCORES[pam[-2:]]
        
        return score * 100

class Azimuth(BatchedSingleSequenceAlgorithm):
    logger = logger.getChild(__qualname__)
    
    def __init__(self):
        super().__init__(
            name="Azimuth",
            authors=['Doench, John G.', 'Fusi, Nicolo', 'Sullender, Meagan', 'Hegde, Mudra', 'Vaimberg, Emma W.', 'Donovan, Katherine F.', 'Smith, Ian', 'Tothova, Zuzana', 'Wilen, Craig', 'Orchard, Robert', 'Virgin, Herbert W.', 'Listgarten, Jennifer', 'Root, David E.'],
            title='Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9',
            journal='Nature Biotechnology',
            issuing='34(2):184-191',
            year=2016,
            doi='https://doi.org/10.1038/nbt.3437',
            #citation="Doench, Fusi, et al. Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. Nature Biotechnology 34, 184-191 (2016).",
            off_target=False,
            on_target=True,
            prefilter=False,
            postfilter=False,
            minimum=1.0,
            maximum=100.0,
            default=None,
            rgn_list=('Cas9',),
            weight_str='Azimuth:52+1.17' # Penalize any score less than 60
        )
    
    # def weight(self, x):
    #     """Penalize any score less than 60"""
    #     return 1.0/(1+1.17**(52-x))
    
    # Not currently used
    def fix_flank(self, upstream='', downstream=''):
        # if upstream is missing any values, then it must be filled with 'Ns'
        #   '' becomes 'NNNN'
        # 'AC' becomes 'NNAC'
        # same with downstream
        
        # upstream needs 4 nt sequence
        us = upstream
        while(len(us) < 4):
            us = 'N' + us
        
        # downstream needs 3 nt sequence
        ds = downstream
        while(len(ds) < 3):
            ds = ds + 'N'
        
        return us, ds
    
    def calculate(self, batch, *args, disambiguate=False, disambiguate_samples=512, **kwargs):
        # Old 'queries2' format: (index, calculate=True/False, (seq1, seq2, ...))
        # New 'queries2' format: (index, calculate=True/False, sequence, default_score=0.0)
        queries2 = []
        
        #if ("disambiguate" in kwargs):
        #    disambiguate = kwargs["disambiguate"]
        #else:
        #    disambiguate = False
        #if ("batch_size" in kwargs):
        #    batch_size = kwargs["batch_size"]
        #else:
        #    batch_size = 50000
        #if ("disambiguate_samples" in kwargs):
        #    disambiguate_samples = max(1, int(kwargs["disambiguate_samples"]))
        #else:
        #    disambiguate_samples = 512
        
#        if sys.platform.startswith('win'):
#            # Command line limit:   32768 characters = 2**15, /32=1024
#            batch_size = 1000
#        else:
#            # Command line limit: 2097152 characters = 2**21, /32=65536
#            batch_size = 10000
        
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
        
        spec = importlib.util.find_spec('azimuth')
        if spec:
            PYTHON = sys.executable
            command_list = [PYTHON, WRAPPER_PATH]
        else:
            PYTHON = 'python'
            if sys.platform.startswith('win'):
                if which('py.exe'):
                    PYTHON = 'py.exe'
                else:
                    ppl = which('python.exe', full=True)
                    for pp in ppl:
                        if 'python27' in pp.lower():
                            PYTHON = pp
                            break
            #print('PYHTON2 =', PYTHON2)
        
            command_list = [PYTHON]
            if (PYTHON == 'py.exe'):
                command_list.append('-2.7')
            command_list.append(WRAPPER_PATH)
        
        #print("queries2", queries2)
        
        batch_list = []
        for q in queries2:
            if q[1]:
                batch_list.append(q[2])
        
        cp = subprocess.run(command_list, input=bytes(' '.join(batch_list), 'utf-8'), shell=False, stdout=subprocess.PIPE)
        out_lines = cp.stdout.decode().splitlines()
        
        batch_scores = []
        for line in out_lines:
            if line.startswith('No model file specified'):
                pass
            else:
                i_seq, i_score = line.split(" ")
                batch_scores.append(100*float(i_score))
        
        
        
        
        
#        batch_list = []
#        batch_count = 0
#        for bi, q in enumerate(queries2):
#            if (batch_count % batch_size == 0):
#                batch_list.append([])
#            if q[1]:
#                batch_list[-1].append(q[2]) # Append the sequence
#                batch_count += 1
        
#        batch_scores = []
#        for current_batch in batch_list:
#            if (len(current_batch) > 0):
#            #for batch_start in range(0, len(queries2), batch_size):
#            #    current_batch = queries2[batch_start:batch_start+batch_size]
#            #    command_list = ['python', WRAPPER_PATH] + [ x[2] for x in current_batch if x[1] ]
#                #command_list = [PYTHON2, WRAPPER_PATH] + current_batch
#                command_list = program_path_list + current_batch
#                
#                #print('command_list =', command_list[:10], '...'+str(len(command_list)) if (len(command_list) > 10) else '')
#                #with open(error_file, 'w+') as flo:
#                #    cp = subprocess.run(command_list, shell=False, check=True, stdout=flo, stderr=subprocess.STDOUT)
#                #print(command_list)
#                cp = subprocess.run(command_list, shell=False, check=True, stdout=subprocess.PIPE)
#                output = cp.stdout.decode()
#                
#                for line in output.splitlines()[1:]: #skip first line: line.startswith('No model file specified')
#                    i_seq, i_score = line.split(" ")
#                    batch_scores.append(100*float(i_score))
                
        # Deal with skipped queries
        queries2_iter = iter(queries2)
        #batch_scores_iter = iter(batch_scores)
        #for s in batch_scores_iter:
        for s in batch_scores:
            q = next(queries2_iter)
            #s = next(batch_scores_iter)
            
            while (q[1] == False):
                self.logger.info("Skipping Azimuth Calculation: {}".format(q))
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

class Elevation(): # Subclass needed
    # Paper:
    #  https://www.nature.com/articles/s41551-017-0178-6
    # Source code:
    #  https://github.com/Microsoft/Elevation
    # Website:
    #  https://www.microsoft.com/en-us/research/project/crispr/
    # Web interface:
    #  https://crispr.ml/

    def __init__(self):
        super().__init__(
            name="Elevation",
            authors=['Listgarten, Jennifer', 'Weinstein, Michael', 'Kleinstiver, Benjamin P.', 'Sousa, Alexander A.', 'Joung, J. Keith', 'Crawford, Jake', 'Gao, Kevin', 'Hoang, Luong', 'Elibol, Melih', 'Doench, John G.', 'Fusi, Nicolo'],
            title='Prediction of off-target activities for the end-to-end design of CRISPR guide RNAs',
            journal='Nature Biomedical Engineering',
            issuing='2(1):38-47',
            year=2018,
            doi='https://doi.org/10.1038/s41551-017-0178-6',
            #citation=("Listgarten, et al. Prediction of off-target activities for the end-to-end design of "
            #          "CRISPR guide RNAs. Nature Biomedical Engineering volume 2, pages 38–47 (2018)"),
            off_target=True,
            on_target=True,
            prefilter=False,
            postfilter=False,
            minimum=1.0,
            maximum=100.0,
            default=None,
            weight_str='Elevation:50+1.17' # Uninformed, placeholder weight
        )

def load_mismatch_scores(file_path, sep="\t", approximation=False):
    """Load the mismatch scores defined by Doench et al (2016)"""
    # Data extracted into tab-delimited text as 'doench_mismatch_scores.txt'
    
    mismatch_scores = {}
    with open(file_path, 'r') as flo:
        for line in flo:
            line = line.rstrip()
            if (len(line) > 0):
                m = regex.match(r'^\s*#', line)
                if not m:
                    sline = line.split(sep)
                    k = 'r' + sline[0] + ':d' + sline[1] + ',' + sline[2]
                    if approximation:
                        mismatch_scores[k] = float(fractions.Fraction(sline[4]))
                    else:
                        mismatch_scores[k] = float(sline[3])
    return mismatch_scores

def load_pam_scores(file_path, sep="\t", approximation=False):
    """
    Load the PAM scors defined by Doench et al (2016)
    These represent the 'NGG' scores, excluding the 'N'
    """
    pam_scores = {}
    with open(file_path, 'r') as flo:
        for line in flo:
            line = line.rstrip()
            if (len(line) > 0):
                m = regex.match(r'^\s*#', line)
                if not m:
                    sline = line.split(sep)
                    if approximation:
                        pam_scores[sline[0]] = float(fractions.Fraction(sline[2]))
                    else:
                        pam_scores[sline[0]] = float(sline[1])
    return pam_scores

def load_old_scores(file_path, sep="\t"):
    """
    Function to open and parse the tab-delimited 'params' file for
    Doench et al (2014), and return a list of tuples"""
    params = []
    with open(file_path, 'r') as flo:
        for line in flo:
            line = line.rstrip()
            if (len(line) > 0):
                m = regex.match(r'^\s*#', line)
                if not m:
                    sline = line.split(sep)
                    params.append((int(sline[0]), sline[1], float(sline[2])))
    return params

# Load scores from the data files when this module is imported
try:
    MISMATCH_SCORES = load_mismatch_scores(os.path.join(os.path.dirname(__file__), 'doench_mismatch_scores.txt'))
    PAM_SCORES = load_pam_scores(os.path.join(os.path.dirname(__file__), 'doench_pam_scores.txt'))
    OLD_SCORES = load_old_scores(os.path.join(os.path.dirname(__file__), 'doench_params.txt'))
except FileNotFoundError:
    raise Exception("Could not find file with mismatch scores or PAM scores")

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
    
    print("=== Doench2014 ===")
    C = Doench2014()
    print(C.calculate(a)) # 49.18420140306892
    print(C.calculate(b)) # 16.668179672077244
    print(C.calculate(c)) # 32.30335056113513
    print(C.calculate(d)) # 43.173287661422385
    print(C.calculate(e)) # 62.564461061593754
    print(C.calculate(f)) # 22.59673843434246
    
    print("=== Doench2016 ===")
    C = Doench2016()
    print(C.calculate(a, a)) # 100.0
    print(C.calculate(a, b)) # 0.0
    print(C.calculate(a, c)) # 0.0
    print(C.calculate(a, d)) # 0.008843903254636097
    print(C.calculate(a, e)) # 21.482277090941643
    print(C.calculate(a, f)) # 42.8571429
    print(C.calculate(j, k)) # 12.26204765820043
    
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
    
    
    import time
    
    start = time.time()
    seqs = []
    for I in range(10000):
        seqs.append(('', random_sequence(20), random_sequence(3), random_sequence(4), random_sequence(3)))
    scores = C.calculate(seqs, disambiguate=True)
    
    print('min={}, max={}, mean={}, time={}'.format(min(scores), max(scores), sum(scores)/len(scores), time.time()-start))
    # Scores range from -4 to 92 with a mean of 52
    
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
    
    figure(scores, C, title='Random Azimuth scores')

if (__name__ == "__main__"):
    test()

#if (__name__ == "__main__"):
#   print(Second().calculate("ACGT"))
