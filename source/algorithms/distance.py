#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/distance.py

# Import standard packages
import os

# Import non-standard packages
import regex

if (__name__ == "__main__"):
    from algorithm import SingleSequenceAlgorithm, PairedSequenceAlgorithm
else:
    from .algorithm import SingleSequenceAlgorithm, PairedSequenceAlgorithm

class GlobalAlignment(object):
    __slots__ = [
        'seq1',
        'seq2',
        'cumulative_scores',
        'align1',
        'align2',
        'bitscore',
        'matches',
        'substitutions',
        'insertions',
        'deletions',
        'errors',
        'gap_opens',
    ]
    previous_alignment = None
    
    def score_extremes(self, scoring, gaps=False):
        """
        Returns the minimum and maximum bitscore possible,
        either ignoring or including potential gaps.
        This function assumes self.seq1 is the motif with IUPAC ambiguities.
        'NGG' for the PAM, for instance.
        """
        iupac = [
            'A', 'C', 'G', 'T',
            'R', 'Y', 'M', 'K', 'W', 'S',
            'B', 'D', 'H', 'V',
            'N',
        ] # note the lack of '-'
        
        minimums = []
        maximums = []
        
        if gaps:
            for a in range(len(self.align1)):
                scores = []
                if (self.align1[a] == '~'):
                    for i in iupac:
                        scores.append(scoring[('-', i, 'open')])
                elif (self.align1[a] == '-'):
                    for i in iupac:
                        scores.append(scoring[('-', i, 'extend')])
                elif (self.align2[a] == '~'):
                    for i in iupac:
                        scores.append(scoring[(i, '-', 'open')])
                elif (self.align2[a] == '-'):
                    for i in iupac:
                        scores.append(scoring[(i, '-', 'extend')])
                else:
                    for i in iupac:
                        scores.append(scoring[(self.align1[a], i)])
                minimums.append(min(scores))
                maximums.append(max(scores))
        else:
            for a in self.seq1:
                scores = []
                for i in iupac:
                    scores.append(scoring[(a, i)])
                minimums.append(min(scores))
                maximums.append(max(scores))
        
        return sum(minimums), sum(maximums)
    
    @classmethod
    def align(cls, seq1, seq2, scoring):
        if cls.previous_alignment:
            if ((seq1 != cls.previous_alignment.seq1) or (seq2 != cls.previous_alignment.seq2)):
                cls.previous_alignment = cls(seq1, seq2, scoring)
        else:
            cls.previous_alignment = cls(seq1, seq2, scoring)
        return cls.previous_alignment
    
    def __init__(self, seq1, seq2, scoring, print_matrices=False):
        """
        Align two sequences globally. Utilizes affine gap penalties.
        Calculates the following for single optimal alignment:
         * aligned sequences,
         * cumulative scores, final bitscore,
         * matches, substitutions, insertions, deletions,
         * total errors,
         * gap opens
        
        Arguments:
          seq1
          seq2
          scoring
        
        Optionally calculates dynamic programming scoring and traceback
        matrices.
        
        Gap opens are designated with '~' and gap extensions with '-'.
        Case is ignored. IUPAC ambiguities are counted as substitutions.
        """
        # Convert to upper-case
        seq1 = seq1.upper()
        seq2 = seq2.upper()
        
        # Get length of input sequences
        len1, len2 = len(seq1), len(seq2)
        
        # Generate DP table and traceback path pointer matrix
        score = self.zeros((len1+1, len2+1))      # the DP table
        choices = [['']*(len2+1) for i in range(len1+1)]
        
        eseq1 = '-' + seq1
        eseq2 = '-' + seq2
        
        # Set the first row and column of the DP table
        #score[0][0] = 0 # Implicit
        for i1 in range(1, len1 + 1):
            if (i1 == 1):
                gap_type = 'open'
            else:
                gap_type = 'extend'
            score[i1][0] = score[i1-1][0] + scoring[(eseq1[i1], eseq2[0], gap_type)]
            if (gap_type == 'open'):
                choices[i1][0] = 'D'
            else:
                choices[i1][0] = 'd'
        for i2 in range(1, len2 + 1):
            if (i2 == 1):
                gap_type = 'open'
            else:
                gap_type = 'extend'
            score[0][i2] = score[0][i2-1] + scoring[(eseq1[0], eseq2[i2], gap_type)]
            if (gap_type == 'open'):
                choices[0][i2] = 'I'
            else:
                choices[0][i2] = 'i'
        
        # Calculate the rest of the DP table
        for i1 in range(1, len1 + 1): # Cycle through rows
            for i2 in range(1, len2 + 1): # Cycle through columns
                # If the best arrow pointing into this cell was a deletion
                if (choices[i1-1][i2] in ['D', 'd']):
                    deletion_gap_type = 'extend'
                else:
                    deletion_gap_type = 'open'
                # If the best arrow pointing into this cell was an insertions
                if (choices[i1][i2-1] in ['I', 'i']):
                    insertion_gap_type = 'extend'
                else:
                    insertion_gap_type = 'open'
                
                match = score[i1-1][i2-1] + scoring[(seq1[i1-1], seq2[i2-1])]
                delete = score[i1-1][i2] + scoring[(seq1[i1-1], '-', deletion_gap_type)]
                insert = score[i1][i2-1] + scoring[('-', seq2[i2-1], insertion_gap_type)]
                
                v = {"M":match, "D":delete, "I":insert}
                best = sorted(v, key=lambda x: v[x], reverse=True)[0]
                
                score[i1][i2] = v[best]
                if (best == 'M'):
                    if (seq1[i1-1] != seq2[i2-1]):
                        best = 'S'
                elif (best == 'I'):
                    if (insertion_gap_type == 'extend'):
                        best = 'i'
                elif (best == 'D'):
                    if (deletion_gap_type == 'extend'):
                        best = 'd'
                choices[i1][i2] = best
        
        if print_matrices:
            # Print the score matrix
            print("Score matrix")
            new_score = []
            temp = [list(eseq2)] + score
            for i, char in enumerate(' '+eseq1):
                new_score.append([char] + temp[i])
            
            for row in new_score:
                sline = []
                for val in row:
                    if isinstance(val, str):
                        sline.append('{:>7}'.format(val))
                    else:
                        sline.append('{:7.2f}'.format(val))
                print(' '.join(sline))
            
            # Print the traceback matrix
            print("Traceback matrix")
            new_choices = []
            temp = [list(eseq2)] + choices
            for i, char in enumerate(' '+eseq1):
                new_choices.append([char] + temp[i])
            
            for r in new_choices:
                print(' '.join(map(lambda x: '{:>2}'.format(x), r)))
        
        # Perform traceback
        count = 0
        i1, i2 = len1, len2
        bitscore = score[i1][i2]
        cumulative_traceback_scores = []
        a1 = ''
        a2 = ''
        matches = 0
        substitutions = 0
        insertions = 0
        deletions = 0
        gap_opens = 0
        while((i1 != 0) or (i2 != 0)):
            if (count > len1+1+len2+1):
                raise Exception("Alignment failed to converge on valid solution")
            cumulative_traceback_scores.append(score[i1][i2])
            T = choices[i1][i2]
            if (T == 'M'):
                a1 += seq1[i1-1]
                a2 += seq2[i2-1]
                i1 -= 1
                i2 -= 1
                matches += 1
            elif (T == 'S'):
                a1 += seq1[i1-1]
                a2 += seq2[i2-1]
                i1 -= 1
                i2 -= 1
                substitutions += 1
            elif (T == 'I'):
                a1 += '~'
                a2 += seq2[i2-1]
                i2 -= 1
                insertions += 1
                gap_opens += 1
            elif (T == 'i'):
                a1 += '-'
                a2 += seq2[i2-1]
                i2 -= 1
                insertions += 1
            elif (T == 'D'):
                a1 += seq1[i1-1]
                a2 += '~'
                i1 -= 1
                deletions += 1
                gap_opens += 1
            elif (T == 'd'):
                a1 += seq1[i1-1]
                a2 += '-'
                i1 -= 1
                deletions += 1
            count += 1
        errors = substitutions + insertions + deletions
        
        # Store the final values as instance attributes
        self.seq1 = seq1
        self.seq2 = seq2
        self.cumulative_scores = cumulative_traceback_scores[::-1]
        self.align1 = a1[::-1]
        self.align2 = a2[::-1]
        self.bitscore = bitscore
        self.matches = matches
        self.substitutions = substitutions
        self.insertions = insertions
        self.deletions = deletions
        self.errors = errors
        self.gap_opens = gap_opens
        
    def __repr__(self):
        return self.__class__.__name__ + '(' + ' '.join([
        'M='+str(self.matches),
        'S='+str(self.substitutions),
        'I='+str(self.insertions),
        'D='+str(self.deletions),
        'E='+str(self.errors),
        'b='+str(self.bitscore),
        # 'c='+str(self.cumulative_scores),
        ]) + ')'
    
    def zeros(self, shape):
        retval = []
        for x in range(shape[0]):
            retval.append([])
            for y in range(shape[1]):
                retval[-1].append(0)
        return retval

def build_score(match=1, transition=-1, transversion=-1.2, gap_open=-4, gap_extend=-2):
    characters = ['-', 'A', 'C', 'G', 'T']
    scoring = {}
    for c1 in characters:
        for c2 in characters:
            if ((c1 == '-') or (c2 == '-')):
                scoring[(c1, c2, 'open')] = gap_open
                scoring[(c1, c2, 'extend')] = gap_extend
            elif (c1 == c2):
                scoring[(c1, c2)] = match
            elif (((c1 in ['A', 'G']) and (c2 in ['A', 'G'])) or ((c1 in ['C', 'T']) and (c2 in ['C', 'T']))):
                scoring[(c1, c2)] = transition
            else:
                scoring[(c1, c2)] = transversion
    return scoring
    
def build_iupac_score(single_scoring):
    """
    Extrapolate input scoring matrix to handle IUPAC ambiguity codes.
    Input argument should be the scoring dict, which can be generated
    with the 'build_score()' function.
    """
    #single_scoring = build_score()
    iupac_scoring = {}
    
    iupac = {
        '-': ['-'],
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
    
    for c1 in iupac.keys():
        for c2 in iupac.keys():
            scores = []
            go_scores = []
            ge_scores = []
            for i1 in iupac[c1]:
                for i2 in iupac[c2]:
                    if ((i1 == '-') or (i2 == '-')):
                        go_scores.append(single_scoring[(i1, i2, 'open')])
                        ge_scores.append(single_scoring[(i1, i2, 'extend')])
                    else:
                        scores.append(single_scoring[(i1, i2)])
            
            if (len(scores) > 0):
                iupac_scoring[(c1, c2)] = sum(scores)/len(scores)
            
            go_concat = scores + go_scores
            if (len(go_scores) > 0):
                iupac_scoring[(c1, c2, 'open')] = sum(go_concat)/len(go_concat)
            
            ge_concat = scores + ge_scores
            if (len(ge_scores) > 0):
                iupac_scoring[(c1, c2, 'extend')] = sum(ge_concat)/len(ge_concat)
    return iupac_scoring

# # Scoring matrix of log-likelihood ratios with weights for
# # transitions/transversions, gap penalties, etc
# scoring = {
#     # Matches
#     ('A', 'A'): 1,
#     ('C', 'C'): 1,
#     ('G', 'G'): 1,
#     ('T', 'T'): 1,
#     # Transitions (Purines: A<->G, Pyrimidines: C<->T)
#     ('A', 'G'): -1,
#     ('G', 'A'): -1,
#     ('C', 'T'): -1,
#     ('T', 'C'): -1,
#     # Transversions
#     ('A', 'T'): -1.2, # Weak
#     ('T', 'A'): -1.2, # Weak
#     ('C', 'G'): -1.2, # Strong
#     ('G', 'C'): -1.2, # Strong
#     ('A', 'C'): -1.2, # Amino
#     ('C', 'A'): -1.2, # Amino
#     ('G', 'T'): -1.2, # Keto
#     ('T', 'G'): -1.2, # Keto
#     # Insertions
#     ('-', 'A', 'open'): -4,
#     ('-', 'C', 'open'): -4,
#     ('-', 'G', 'open'): -4,
#     ('-', 'T', 'open'): -4,
#     ('-', 'A', 'extend'): -2,
#     ('-', 'C', 'extend'): -2,
#     ('-', 'G', 'extend'): -2,
#     ('-', 'T', 'extend'): -2,
#     # Deletions
#     ('A', '-', 'open'): -4,
#     ('C', '-', 'open'): -4,
#     ('G', '-', 'open'): -4,
#     ('T', '-', 'open'): -4,
#     ('A', '-', 'extend'): -2,
#     ('C', '-', 'extend'): -2,
#     ('G', '-', 'extend'): -2,
#     ('T', '-', 'extend'): -2,
#     # Gaps should never be aligned to each other
#     ('-', '-', 'open'): -4,
#     ('-', '-', 'extend'): -2,
# }
# 
# scoring = build_iupac_score(scoring)

def load_scores(file_path, sep="\t"):
    """Load the PAM scors defined by Doench et al (2016)
    These represent the 'NGG' scores, excluding the 'N'
    """
    scoring = {}
    with open(file_path, 'r') as flo:
        for line in flo:
            line = line.rstrip()
            if (len(line) > 0):
                m = regex.match(r'^\s*#', line)
                if not m:
                    sline = line.split(sep)
                    if (sline[2] == ''):
                        scoring[(sline[0], sline[1])] = float(sline[3])
                    else:
                        scoring[(sline[0], sline[1], sline[2])] = float(sline[3])
    return scoring

class Substitutions(PairedSequenceAlgorithm):
    def __init__(self):
        super().__init__("Substitutions", "Seher", 2017,
            citation="AddTag",
            off_target=False,
            prefilter=False,
            postfilter=True,
            minimum=0.0,
            maximum=5.0,
            default=0.0
        )
    
    def calculate(self, intended, potential, *args, **kwargs):
        on_sequence, on_target, on_pam, on_upstream, on_downstream = intended
        off_sequence, off_target, off_pam, off_upstream, off_downstream = potential
        
        return self.score(on_target, off_target)
    
    def score(self, seq1, seq2):
        """"""
        return GlobalAlignment.align(seq1, seq2, SCORES).substitutions

class Insertions(PairedSequenceAlgorithm):
    def __init__(self):
        super().__init__("Insertions", "Seher", 2017,
            citation="AddTag",
            off_target=False,
            prefilter=False,
            postfilter=True,
            minimum=0.0,
            maximum=2.0,
            default=0.0
        )
    
    def calculate(self, intended, potential, *args, **kwargs):
        on_sequence, on_target, on_pam, on_upstream, on_downstream = intended
        off_sequence, off_target, off_pam, off_upstream, off_downstream = potential
        
        return self.score(on_target, off_target)
    
    def score(self, seq1, seq2):
        """"""
        return GlobalAlignment.align(seq1, seq2, SCORES).insertions

class Deletions(PairedSequenceAlgorithm):
    def __init__(self):
        super().__init__("Deletions", "Seher", 2017,
            citation="AddTag",
            off_target=False,
            prefilter=False,
            postfilter=True,
            minimum=0.0,
            maximum=2.0,
            default=0.0
        )
    
    def calculate(self, intended, potential, *args, **kwargs):
        on_sequence, on_target, on_pam, on_upstream, on_downstream = intended
        off_sequence, off_target, off_pam, off_upstream, off_downstream = potential
        
        return self.score(on_target, off_target)
    
    def score(self, seq1, seq2):
        """"""
        return GlobalAlignment.align(seq1, seq2, SCORES).deletions

class Errors(PairedSequenceAlgorithm):
    def __init__(self):
        super().__init__("Errors", "Seher", 2017,
            citation="AddTag",
            off_target=False,
            prefilter=False,
            postfilter=True,
            minimum=0.0,
            maximum=5.0,
            default=0.0
        )
    
    def calculate(self, intended, potential, *args, **kwargs):
        on_sequence, on_target, on_pam, on_upstream, on_downstream = intended
        off_sequence, off_target, off_pam, off_upstream, off_downstream = potential
        
        return self.score(on_target, off_target)
    
    def score(self, seq1, seq2):
        """"""
        return GlobalAlignment.align(seq1, seq2, SCORES).errors

class PamIdentity(SingleSequenceAlgorithm):
    def __init__(self):
        """
        Quantifies the similarity between the PAM motif and the actual PAM
        of the target. Normalizes based on gapless minimum and maximum
        possible bitscores.
        """
        super().__init__("PAM-Identity", "Seher", 2017,
            citation="AddTag",
            off_target=False,
            prefilter=False, # Filter spacers before aligning
            postfilter=True, # Filter alignments before calculating scores
            minimum=75.0,
            maximum=100.0,
            default=None
        )
    
    def calculate(self, potential, *args, **kwargs):
        off_sequence, off_target, off_pam, off_upstream, off_downstream = potential
        
        return self.score(off_pam, kwargs['motif'])
    
    def score(self, pam, motif):
        """"""
        # Separate spacer and PAM motifs
        if ('>' in motif):
            spacer_motif, pam_motif = motif.split('>')
        elif ('<' in motif):
            pam_motif, spacer_motif = motif.split('<')
        #else:
        #    raise Exception()
        
        a = GlobalAlignment.align(pam_motif, pam, SCORES)
        min_score, max_score = a.score_extremes(SCORES)
        
        return (a.bitscore - min_score)/(max_score - min_score)*100

# Load scores from the data files when this module is imported
try:
    SCORES = load_scores(os.path.join(os.path.dirname(__file__), 'distance_scores.txt'))
    SCORES = build_iupac_score(SCORES)
except FileNotFoundError:
    raise Exception("Could not find file with scores")

def test():
    """Code to test the functions and classes"""
    
    a = ('',  'ACGTYGGCCATTKAGGCAg', 'TGG', '', '')
    b = ('', 'ACATAGGTCCACTTAGGCaG', 'TGG', '', '')
    c = ('', 'GCATTGCCACTTTGGCAGAG', 'AGG', '', '')
    
    print("=== Substitutions ===")
    C = Substitutions()
    print(C.calculate(a, b))
    print(C.calculate(a, c))
    
    print("=== Insertions ===")
    C = Insertions()
    print(C.calculate(a, b))
    print(C.calculate(a, c))
    
    print("=== Deletions ===")
    C = Deletions()
    print(C.calculate(a, b))
    print(C.calculate(a, c))
    
    print("=== Errors ===")
    C = Errors()
    print(C.calculate(a, b))
    print(C.calculate(a, c))
    
    print("=== PAM-Identity ===")
    C = PamIdentity()
    print(C.calculate(('', '', 'GGG', '', ''), motif='N{{20}}>NGG'))
    print(C.calculate(('', '', 'GAG', '', ''), motif='N{{20}}>NGG'))
    print(C.calculate(('', '', 'GAG', '', ''), motif='N{{20}}>NRG'))
    print(C.calculate(('', '', 'CCACAAA', '', ''), motif='N{{20}}>NHRBMAW'))
    print(C.calculate(('', '', 'ATGCCAT', '', ''), motif='N{{20}}>NHRBMAW'))
    print(C.calculate(('', '', 'CGACAAA', '', ''), motif='N{{20}}>NHRBMAW'))
    print(C.calculate(('', '', 'CGACAAC', '', ''), motif='N{{20}}>NHRBMAW'))
    print(C.calculate(('', '', 'CCACAGA', '', ''), motif='N{{20}}>NHRBMAW'))
    print(C.calculate(('', '', 'CCCCAAG', '', ''), motif='N{{20}}>NHRBMAW'))
    
    print("=== Test ===")
    print(GlobalAlignment.align('GGG', 'GGG', SCORES))
    print(GlobalAlignment.align('GGG', 'CCC', SCORES))
    print(GlobalAlignment.align('NRG', 'AGG', SCORES))
    print(GlobalAlignment.align('NRG', 'TCG', SCORES))
    print(GlobalAlignment.align('NRG', 'NRG', SCORES))
    print(GlobalAlignment.align('NRG', 'NNN', SCORES))
    print(GlobalAlignment.align('NRG', 'RRR', SCORES))
    print(GlobalAlignment.align('N', 'N', SCORES), 'N', 'N')
    print(GlobalAlignment.align('N', 'A', SCORES), 'N', 'A') #-0.6
    print(GlobalAlignment.align('N', 'B', SCORES), 'N', 'B')
    print(GlobalAlignment.align('R', 'N', SCORES), 'R', 'N')
    print(GlobalAlignment.align('R', 'R', SCORES), 'R', 'R')
    print(GlobalAlignment.align('R', 'G', SCORES), 'R', 'G')
    print(GlobalAlignment.align('R', 'T', SCORES), 'R', 'T')
    
    print('------')
    a = GlobalAlignment.align('NARRG', 'TAAG', SCORES)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(SCORES))
    print(a.score_extremes(SCORES, gaps=True), 'gaps')
    
    print('------')
    a = GlobalAlignment.align('NARG', 'TTAAG', SCORES)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(SCORES))
    print(a.score_extremes(SCORES, gaps=True), 'gaps')
    
    print('------')
    a = GlobalAlignment.align('NGG', 'GGG', SCORES)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(SCORES))
    print(a.score_extremes(SCORES, gaps=True), 'gaps')
    
    print('------')
    a = GlobalAlignment.align('NGG', 'ARG', SCORES)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(SCORES))
    print(a.score_extremes(SCORES, gaps=True), 'gaps')
    
    print('------')
    a = GlobalAlignment.align('NRG', 'AGG', SCORES)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(SCORES))
    print(a.score_extremes(SCORES, gaps=True), 'gaps')
    
    print('------')
    a = GlobalAlignment.align('NGG', 'ACG', SCORES)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(SCORES))
    print(a.score_extremes(SCORES, gaps=True), 'gaps')
    
    print('------')
    a = GlobalAlignment.align('NHRBMAW', 'CCACAAA', SCORES)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(SCORES))
    print(a.score_extremes(SCORES, gaps=True), 'gaps')
    
    print('------')
    a = GlobalAlignment.align('NHRBMAW', 'CGACAAA', SCORES)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(SCORES))
    print(a.score_extremes(SCORES, gaps=True), 'gaps')

if (__name__ == "__main__"):
    test()
