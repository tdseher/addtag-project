#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/doench.py

# Import standard packages
import sys
import os
import math
import fractions

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
    from nucleotides import rc
else:
    from .algorithm import SingleSequenceAlgorithm, PairedSequenceAlgorithm, BatchedSingleSequenceAlgorithm
    from ..nucleotides import rc

class Doench2014(SingleSequenceAlgorithm):
    def __init__(self):
        super().__init__("Doench-2014", "Doench, et al", 2014,
            citation="Doench, et al. Rational design of highly active sgRNAs for CRISPR-Cas9–mediated gene inactivation. Nature Biotechnology 32, 1262–1267 (2014)",
            off_target=True,
            prefilter=False,
            postfilter=True,
            minimum=1.0,
            maximum=100.0,
            default=None
        )
    
    def calculate(self, intended, *args, **kwargs):
        sequence, target, pam, upstream, downstream = intended
        
        return self.score(target, pam, upstream, downstream)
    
    def score(self, seq, pam, upstream='', downstream=''):
        """Function to calculate the sgRNA on-target efficacy score,
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
        super().__init__("CFD", "Doench, Fusi, et al", 2016,
            citation="Doench, Fusi, et al. Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. Nature Biotechnology 34, 184–191 (2016).",
            off_target=True,
            prefilter=False,
            postfilter=True,
            minimum=1.0,
            maximum=100.0,
            default=100.0
        )
    
    def calculate(self, intended, potential, *args, **kwargs):
        on_sequence, on_target, on_pam, on_upstream, on_downstream = intended
        off_sequence, off_target, off_pam, off_upstream, off_downstream = potential
        
        return self.score(on_target, off_target, on_pam)
    
    #def calculate(self, seq1, seq2, pam, *args, **kwargs):
    def score(self, seq1, seq2, pam):
        """Calculate the CFD score for a given gRNA and off-target sequence.
        
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
        m1 = regex.search('[^ATCGatcg]', seq1)
        m2 = regex.search('[^ATCGatcg]', seq2)
        mp = regex.search('[^ATCGatcg]', pam[-2:])
        if ((m1 != None) or (m2 != None) or (mp != None)):
            # Otherwise return a score of 0
            return 0.0
        
        # Start with a perfect score
        score = 1
        
        # Ensure inputs are RNA sequences
        seq2 = seq2.replace('T','U')
        seq1 = seq1.replace('T','U')
        
        # modify algorithm to allow for less than 20 nt gRNAs
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
    def __init__(self):
        super().__init__("Azimuth", "Doench, Fusi, et al", 2016,
            citation="Doench, Fusi, et al. Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9. Nature Biotechnology 34, 184-191 (2016).",
            off_target=False,
            prefilter=False,
            postfilter=True,
            minimum=1.0,
            maximum=100.0,
            default=None
        )
    
    def calculate(self, batch, *args, **kwargs):
        queries = []
        skips = {}
        
        # unpack the input sequences
        for i, query in enumerate(batch):
            sequence, target, pam, upstream, downstream = query
            built_query = self.build_query(target, pam, upstream, downstream)
            m = regex.search('[^ACGT]', built_query)
            if m:
                skips[i] = built_query
            else:
                queries.append(built_query)
        
        # Obtain path of currently-running file
        SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
        WRAPPER_PATH = os.path.join(SCRIPT_DIR, "azimuth_wrapper.py")
        # Should use args.python2_path
        
        o_scores = []
        line_count = 0
        if (len(queries) > 0):
            command_list = ['python', WRAPPER_PATH] + queries
            
            #with open(error_file, 'w+') as flo:
            #    cp = subprocess.run(command_list, shell=False, check=True, stdout=flo, stderr=subprocess.STDOUT)
            cp = subprocess.run(command_list, shell=False, check=True, stdout=subprocess.PIPE)
            output = cp.stdout.decode()
        
            for line in output.splitlines():
                if not line.startswith('No model file specified'):
                    while line_count in skips:
                        o_scores.append(0.0)
                        line_count += 1
                    i_seq, i_score = line.split(" ")
                    o_scores.append(100*float(i_score))
                    line_count += 1
        for i in range(len(batch)-line_count):
            o_scores.append(0.0)
        
        return o_scores
    
    def build_query(self, seq, pam, upstream='', downstream=''):
        # us     seq                pam  ds
        # ACAG CTGATCTCCAGATATGACCA|TGG GTT
        # CAGC TGATCTCCAGATATGACCAT|GGG TTT
        # CCAG AAGTTTGAGCCACAAACCCA|TGG TCA
        
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
    """Load the PAM scors defined by Doench et al (2016)
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
    """Function to open and parse the tab-delimited 'params' file for
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
    
    a = ('', 'CGATGGCTAGGATCGATTGA', 'TGG', '', '')
    b = ('', 'RYMKWSACGTbDHVNACGTA', 'TGG', '', '')
    c = ('',   'ATGSCTCGGATCGATTGA', 'AGG', '', '')
    d = ('', 'GCGATGCGCAGCTAGGCCGG', 'CGG', '', '')
    e = ('', 'CGAAGGCTCGGACCGATTGA', 'GGG', '', '')
    f = ('', 'CGCTGGCTAGGATCGATTGA', 'AGG', '', '')
    g = ('', 'AAAATTAACTATAGGTAAAG', 'TGG', '', '')
    h = ('', 'AACATCAACTCTAGCTAACG', 'CGG', '', '')
    i = ('', 'AACATCAACTCTACCTAACG', 'CGG', 'CCGA', 'AACA')
    
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
    
    print("=== Azimuth ===")
    C = Azimuth()
    print(C.calculate([g])) # [0.0]
    print(C.calculate([h])) # [0.0]
    print(C.calculate([g, h, i])) # [0.0, 0.0, 68.10224199640001]

if (__name__ == "__main__"):
    test()

#if (__name__ == "__main__"):
#   print(Second().calculate("ACGT"))