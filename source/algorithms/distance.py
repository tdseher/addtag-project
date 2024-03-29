#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/distance.py

# Import standard packages
import os

# Import non-standard packages
import regex

if (__name__ == "__main__"):
    from algorithm import SingleSequenceAlgorithm, PairedSequenceAlgorithm
    from distance_matrix import load_scores, build_iupac_score, build_score
else:
    from .algorithm import SingleSequenceAlgorithm, PairedSequenceAlgorithm
    from .distance_matrix import load_scores, build_iupac_score, build_score

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
        #iupac = [
        #    'A', 'C', 'G', 'T',
        #    'R', 'Y', 'M', 'K', 'W', 'S',
        #    'B', 'D', 'H', 'V',
        #    'N',
        #] # note the lack of '-'
        iupac = list(set([item for sublist in scoring for item in sublist if not item in ('open', 'extend', '-')]))
        
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
    def align(cls, seq1, seq2, scoring, print_matrices=False):
        if cls.previous_alignment:
            if ((seq1 != cls.previous_alignment.seq1) or (seq2 != cls.previous_alignment.seq2)):
                cls.previous_alignment = cls(seq1, seq2, scoring, print_matrices=print_matrices)
        else:
            cls.previous_alignment = cls(seq1, seq2, scoring, print_matrices=print_matrices)
        return cls.previous_alignment
    
    def get_highest_score(self, m, d, i):
        # MATCH = 0b1
        # O_DEL = 0b100
        # O_INS = 0b10000
        if ((m > d) and (m > i)):
            return 0b1, m
        elif ((d > m) and (d > i)):
            return 0b100, d
        elif ((i > m) and (i > d)):
            return 0b10000, i
        elif ((m == d) and (m > i)):
            return 0b1, m
        elif ((m == i) and (m > d)):
            return 0b1, m
        elif ((d == i) and (d > m)):
            return 0b100, d
        else:
            return 0b1, m
    
    def __init__(self, seq1, seq2, scoring, print_matrices=False):
        """
        Align two sequences globally using the Needleman-Wunsch algorithm.
        Utilizes affine gap penalties. Calculates the following for single
        optimal alignment:
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
        # Binary conversion tables for 'choices' table
        NONE = 0b0 # Nothing
        MATCH = 0b1 # A match
        SUB = 0b10 # A substitution
        O_DEL = 0b100 # Start of a new deletion
        C_DEL = 0b1000 # Continue an existing deletion
        O_INS = 0b10000 # Start a new insertion
        C_INS = 0b100000 # Continue an existing insertion
        # char2bin = {
        #     '': NONE,
        #     'M': MATCH,
        #     'S': SUB,
        #     'D': O_DEL,
        #     'd': C_DEL,
        #     'I': O_INS,
        #     'i': C_INS,
        # }
        
        # Convert to upper-case
        seq1 = seq1.upper()
        seq2 = seq2.upper()
        
        # Get length of input sequences
        len1, len2 = len(seq1), len(seq2)
        
        # Generate DP table and traceback path pointer matrix
        score = self.zeros((len1+1, len2+1))      # the DP table
        choices = [[NONE]*(len2+1) for i in range(len1+1)]
        
        eseq1 = '-' + seq1
        eseq2 = '-' + seq2
        
        # Set the first row and column of the DP table
        #score[0][0] = 0 # Implicit
        for i1 in range(1, len1 + 1):
            if (i1 == 1):
                gap_type = 'open'
                choices[i1][0] = O_DEL
            else:
                gap_type = 'extend'
                choices[i1][0] = C_DEL
            score[i1][0] = score[i1-1][0] + scoring[(eseq1[i1], eseq2[0], gap_type)]
        
        for i2 in range(1, len2 + 1):
            if (i2 == 1):
                gap_type = 'open'
                choices[0][i2] = O_INS
            else:
                gap_type = 'extend'
                choices[0][i2] = C_INS
            score[0][i2] = score[0][i2-1] + scoring[(eseq1[0], eseq2[i2], gap_type)]
        
        # Calculate the rest of the DP table
        for i1 in range(1, len1 + 1): # Cycle through rows
            for i2 in range(1, len2 + 1): # Cycle through columns
                # If the best arrow pointing into this cell was a deletion
                if (choices[i1-1][i2] in [O_DEL, C_DEL]):
                    deletion_gap_type = 'extend'
                else:
                    deletion_gap_type = 'open'
                # If the best arrow pointing into this cell was an insertions
                if (choices[i1][i2-1] in [O_INS, C_INS]):
                    insertion_gap_type = 'extend'
                else:
                    insertion_gap_type = 'open'
                
                match = score[i1-1][i2-1] + scoring[(seq1[i1-1], seq2[i2-1])]
                delete = score[i1-1][i2] + scoring[(seq1[i1-1], '-', deletion_gap_type)]
                insert = score[i1][i2-1] + scoring[('-', seq2[i2-1], insertion_gap_type)]
                
                
                
                #v = {"M":match, "D":delete, "I":insert}
                #best = sorted(v, key=lambda x: v[x], reverse=True)[0]
                #score[i1][i2] = v[best]
                
                best, score[i1][i2] = self.get_highest_score(match, delete, insert)
                
                if (best == MATCH):
                    if (seq1[i1-1] != seq2[i2-1]):
                        best = SUB
                elif (best == O_INS):
                    if (insertion_gap_type == 'extend'):
                        best = C_INS
                elif (best == O_DEL):
                    if (deletion_gap_type == 'extend'):
                        best = C_DEL
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
            
            bin2char = {
                NONE: '',
                MATCH: 'M',
                SUB: 'S',
                O_DEL: 'D',
                C_DEL: 'd',
                O_INS: 'I',
                C_INS: 'i',
            }
            
            choices2 = []
            for r in choices:
                xx = []
                for x in r:
                    xx.append(bin2char[x])
                choices2.append(xx)
            
            temp = [list(eseq2)] + choices2
            for i, char in enumerate(' '+eseq1):
                new_choices.append([char] + temp[i])
            
            for r in new_choices:
                print(' '.join(map(lambda x: '{:>2}'.format(x), r)))
            
            # TODO: Defer printing out the 'Traceback matrix' until the traceback is performed.
            #       Then include '-', '\', and '|' for insertion, match/sub, and deletion in the
            #       table to show the path.
        
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
            if (T == MATCH):
                a1 += seq1[i1-1]
                a2 += seq2[i2-1]
                i1 -= 1
                i2 -= 1
                matches += 1
            elif (T == SUB):
                a1 += seq1[i1-1]
                a2 += seq2[i2-1]
                i1 -= 1
                i2 -= 1
                substitutions += 1
            elif (T == O_INS):
                a1 += '~'
                a2 += seq2[i2-1]
                i2 -= 1
                insertions += 1
                gap_opens += 1
            elif (T == C_INS):
                a1 += '-'
                a2 += seq2[i2-1]
                i2 -= 1
                insertions += 1
            elif (T == O_DEL):
                a1 += seq1[i1-1]
                a2 += '~'
                i1 -= 1
                deletions += 1
                gap_opens += 1
            elif (T == C_DEL):
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

class Substitutions(PairedSequenceAlgorithm):
    def __init__(self):
        super().__init__(
            name="Substitutions",
            authors=["Needleman, Saul B.", "Wunsch, Christian D."],
            title='A general method applicable to the search for similarities in the amino acid sequence of two proteins',
            journal='Journal of Molecular Biology',
            issuing='48(3):443-453',
            year=1970,
            doi='https://doi.org/10.1016/0022-2836(70)90057-4',
            off_target=False,
            prefilter=False,
            postfilter=True,
            minimum=0.0,
            maximum=5.0,
            default=0.0,
            weight_str=None
        )
    
    def is_available(self):
        """
        Determines if the prerequisites for the Algorithm have been met.
        :return: True or False
        """
        return True
    
    def calculate(self, intended, potential, *args, **kwargs):
        on_sequence, on_side, on_target, on_pam, on_upstream, on_downstream = intended
        off_sequence, off_side, off_target, off_pam, off_upstream, off_downstream = potential
        
        return self.score(on_target, off_target)
    
    def score(self, seq1, seq2):
        """"""
        return GlobalAlignment.align(seq1, seq2, SCORES).substitutions

class Insertions(PairedSequenceAlgorithm):
    def __init__(self):
        super().__init__(
            name="Insertions",
            authors=["Needleman, Saul B.", "Wunsch, Christian D."],
            title='A general method applicable to the search for similarities in the amino acid sequence of two proteins',
            journal='Journal of Molecular Biology',
            issuing='48(3):443-453',
            year=1970,
            doi='https://doi.org/10.1016/0022-2836(70)90057-4',
            off_target=False,
            prefilter=False,
            postfilter=True,
            minimum=0.0,
            maximum=2.0,
            default=0.0,
            weight_str=None
        )
    
    def is_available(self):
        """
        Determines if the prerequisites for the Algorithm have been met.
        :return: True or False
        """
        return True
    
    def calculate(self, intended, potential, *args, **kwargs):
        on_sequence, on_side, on_target, on_pam, on_upstream, on_downstream = intended
        off_sequence, off_side, off_target, off_pam, off_upstream, off_downstream = potential
        
        return self.score(on_target, off_target)
    
    def score(self, seq1, seq2):
        """"""
        return GlobalAlignment.align(seq1, seq2, SCORES).insertions

class Deletions(PairedSequenceAlgorithm):
    def __init__(self):
        super().__init__(
            name="Deletions",
            authors=["Needleman, Saul B.", "Wunsch, Christian D."],
            title='A general method applicable to the search for similarities in the amino acid sequence of two proteins',
            journal='Journal of Molecular Biology',
            issuing='48(3):443-453',
            year=1970,
            doi='https://doi.org/10.1016/0022-2836(70)90057-4',
            off_target=False,
            prefilter=False,
            postfilter=True,
            minimum=0.0,
            maximum=2.0,
            default=0.0,
            weight_str=None
        )
    
    def is_available(self):
        """
        Determines if the prerequisites for the Algorithm have been met.
        :return: True or False
        """
        return True
    
    def calculate(self, intended, potential, *args, **kwargs):
        on_sequence, on_side, on_target, on_pam, on_upstream, on_downstream = intended
        off_sequence, off_side, off_target, off_pam, off_upstream, off_downstream = potential
        
        return self.score(on_target, off_target)
    
    def score(self, seq1, seq2):
        """"""
        return GlobalAlignment.align(seq1, seq2, SCORES).deletions

class Errors(PairedSequenceAlgorithm):
    def __init__(self):
        super().__init__(
            name="Errors",
            authors=["Needleman, Saul B.", "Wunsch, Christian D."],
            title='A general method applicable to the search for similarities in the amino acid sequence of two proteins',
            journal='Journal of Molecular Biology',
            issuing='48(3):443-453',
            year=1970,
            doi='https://doi.org/10.1016/0022-2836(70)90057-4',
            off_target=False,
            prefilter=False,
            postfilter=True,
            minimum=0.0,
            maximum=5.0,
            default=0.0,
            weight_str=None
        )
    
    def is_available(self):
        """
        Determines if the prerequisites for the Algorithm have been met.
        :return: True or False
        """
        return True
    
    def calculate(self, intended, potential, *args, **kwargs):
        on_sequence, on_side, on_target, on_pam, on_upstream, on_downstream = intended
        off_sequence, off_side, off_target, off_pam, off_upstream, off_downstream = potential
        
        return self.score(on_target, off_target)
    
    def score(self, seq1, seq2):
        """"""
        return GlobalAlignment.align(seq1, seq2, SCORES).errors

class PamIdentity(SingleSequenceAlgorithm):
    def __init__(self):
        """
        Quantifies the similarity between the PAM motif and the actual PAM
        of the target. Normalizes based on gapless minimum and maximum
        possible bitscores. Returns the maximum for all pams that follow the
        motif.
        """
        super().__init__(
            name="PAM-Identity",
            authors=["Seher, Thaddeus D."],
            title='[AddTag paper]',
            journal='',
            issuing='',
            year=2017,
            doi='',
            off_target=False,
            prefilter=False, # Filter spacers before aligning
            postfilter=True, # Filter alignments before calculating scores
            minimum=75.0,
            maximum=100.0,
            default=None,
            weight_str=None
        )
    
    def is_available(self):
        """
        Determines if the prerequisites for the Algorithm have been met.
        :return: True or False
        """
        return True
    
    def calculate(self, potential, *args, **kwargs):
        off_sequence, off_side, off_target, off_pam, off_upstream, off_downstream = potential
        
        return self.score(off_pam, kwargs['parsed_motif'])
    
    def score(self, pam, parsed_motif):
        """"""
        # Store each score in this list
        scores = []
        
        # Separate spacer and PAM motifs
        spacers, pams, side = parsed_motif
        for p in pams:
            a = GlobalAlignment.align(p, pam, SCORES)
            min_score, max_score = a.score_extremes(SCORES)
            scores.append((a.bitscore - min_score)/(max_score - min_score)*100)
        
        # Return the max score
        return max(scores)

# Build test substitution matrix
# SCORES = build_iupac_score(build_score())

# Load substitution matrix scores from the data files when this module is imported
try:
    SCORES = load_scores(os.path.join(os.path.dirname(__file__), 'distance_scores.txt'))
    SCORES = build_iupac_score(SCORES)
except FileNotFoundError:
    raise Exception("Could not find file with scores")

def test():
    """Code to test the functions and classes"""
    
    a = ('', '>',  'ACGTYGGCCATTKAGGCAg', 'TGG', '', '')
    b = ('', '>', 'ACATAGGTCCACTTAGGCaG', 'TGG', '', '')
    c = ('', '>', 'GCATTGCCACTTTGGCAGAG', 'AGG', '', '')
    
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
#    print(C.calculate(('', '', 'GGG', '', ''), motif='N{{20}}>NGG'))
#    print(C.calculate(('', '', 'GAG', '', ''), motif='N{{20}}>NGG'))
#    print(C.calculate(('', '', 'GAG', '', ''), motif='N{{20}}>NRG'))
#    print(C.calculate(('', '', 'CCACAAA', '', ''), motif='N{{20}}>NHRBMAW'))
#    print(C.calculate(('', '', 'ATGCCAT', '', ''), motif='N{{20}}>NHRBMAW'))
#    print(C.calculate(('', '', 'CGACAAA', '', ''), motif='N{{20}}>NHRBMAW'))
#    print(C.calculate(('', '', 'CGACAAC', '', ''), motif='N{{20}}>NHRBMAW'))
#    print(C.calculate(('', '', 'CCACAGA', '', ''), motif='N{{20}}>NHRBMAW'))
#    print(C.calculate(('', '', 'CCCCAAG', '', ''), motif='N{{20}}>NHRBMAW'))
#    print(C.calculate(('', '', 'AAATGG', '', ''), motif='N{{21,22}}>AAANGG'))
#    print(C.calculate(('', '', 'AAATGG', '', ''), motif='N{{21,22}}>A{{3,4}}NGG'))
    print(C.calculate(('', '>', '', 'GGG', '', ''), parsed_motif=(['NNNNNNNNNNNNNNNNNNNN'], ['NGG'], '>')))
    print(C.calculate(('', '>', '', 'GAG', '', ''), parsed_motif=(['NNNNNNNNNNNNNNNNNNNN'], ['NGG'], '>')))
    print(C.calculate(('', '>', '', 'GAG', '', ''), parsed_motif=(['NNNNNNNNNNNNNNNNNNNN'], ['NRG'], '>')))
    print(C.calculate(('', '>', '', 'CCACAAA', '', ''), parsed_motif=(['NNNNNNNNNNNNNNNNNNNN'], ['NHRBMAW'], '>')))
    print(C.calculate(('', '>', '', 'ATGCCAT', '', ''), parsed_motif=(['NNNNNNNNNNNNNNNNNNNN'], ['NHRBMAW'], '>')))
    print(C.calculate(('', '>', '', 'CGACAAA', '', ''), parsed_motif=(['NNNNNNNNNNNNNNNNNNNN'], ['NHRBMAW'], '>')))
    print(C.calculate(('', '>', '', 'CGACAAC', '', ''), parsed_motif=(['NNNNNNNNNNNNNNNNNNNN'], ['NHRBMAW'], '>')))
    print(C.calculate(('', '>', '', 'CCACAGA', '', ''), parsed_motif=(['NNNNNNNNNNNNNNNNNNNN'], ['NHRBMAW'], '>')))
    print(C.calculate(('', '>', '', 'CCCCAAG', '', ''), parsed_motif=(['NNNNNNNNNNNNNNNNNNNN'], ['NHRBMAW'], '>')))
    print(C.calculate(('', '>', '', 'AAATGG', '', ''), parsed_motif=(['NNNNNNNNNNNNNNNNNNNNN', 'NNNNNNNNNNNNNNNNNNNNNN'], ['AAANGG'], '>')))
    print(C.calculate(('', '>', '', 'AAATGG', '', ''), parsed_motif=(['NNNNNNNNNNNNNNNNNNNNN', 'NNNNNNNNNNNNNNNNNNNNNN'], ['AAANGG', 'AAAANGG'], '>')))
    print(C.calculate(('', '>', '', 'AAAATGG', '', ''), parsed_motif=(['NNNNNNNNNNNNNNNNNNNNN', 'NNNNNNNNNNNNNNNNNNNNNN'], ['AAANGG', 'AAAANGG'], '>')))
    print(C.calculate(('', '>', '', 'AATATGG', '', ''), parsed_motif=(['NNNNNNNNNNNNNNNNNNNNN', 'NNNNNNNNNNNNNNNNNNNNNN'], ['AAANGG', 'AAAANGG'], '>')))
    print(C.calculate(('', '>', '', 'CAAATGG', '', ''), parsed_motif=(['NNNNNNNNNNNNNNNNNNNNN', 'NNNNNNNNNNNNNNNNNNNNNN'], ['AAANGG', 'AAAANGG'], '>')))
    print(C.calculate(('', '>', '', 'AAAATCG', '', ''), parsed_motif=(['NNNNNNNNNNNNNNNNNNNNN', 'NNNNNNNNNNNNNNNNNNNNNN'], ['AAANGG', 'AAAANGG'], '>')))
    print(C.calculate(('', '>', '', 'AAAATTG', '', ''), parsed_motif=(['NNNNNNNNNNNNNNNNNNNNN', 'NNNNNNNNNNNNNNNNNNNNNN'], ['AAANGG', 'AAAANGG'], '>')))
    print(C.calculate(('', '>', '', 'AAAATAG', '', ''), parsed_motif=(['NNNNNNNNNNNNNNNNNNNNN', 'NNNNNNNNNNNNNNNNNNNNNN'], ['AAANGG', 'AAAANGG'], '>')))
    
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
    
    print('=======================================')
    a = GlobalAlignment.align('NARRG', 'TAAG', SCORES)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(SCORES))
    print(a.score_extremes(SCORES, gaps=True), 'gaps')
    
    print('=======================================')
    a = GlobalAlignment.align('NARG', 'TTAAG', SCORES)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(SCORES))
    print(a.score_extremes(SCORES, gaps=True), 'gaps')
    
    print('=======================================')
    a = GlobalAlignment.align('NGG', 'GGG', SCORES)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(SCORES))
    print(a.score_extremes(SCORES, gaps=True), 'gaps')
    
    print('=======================================')
    a = GlobalAlignment.align('NGG', 'ARG', SCORES)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(SCORES))
    print(a.score_extremes(SCORES, gaps=True), 'gaps')
    
    print('=======================================')
    a = GlobalAlignment.align('NRG', 'AGG', SCORES)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(SCORES))
    print(a.score_extremes(SCORES, gaps=True), 'gaps')
    
    print('=======================================')
    a = GlobalAlignment.align('NGG', 'ACG', SCORES)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(SCORES))
    print(a.score_extremes(SCORES, gaps=True), 'gaps')
    
    print('=======================================')
    a = GlobalAlignment.align('NHRBMAW', 'CCACAAA', SCORES, print_matrices=False)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(SCORES))
    print(a.score_extremes(SCORES, gaps=True), 'gaps')
    
    print('=======================================')
    a = GlobalAlignment.align('NHRBMAW', 'CGACAAA', SCORES, print_matrices=False)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(SCORES))
    print(a.score_extremes(SCORES, gaps=True), 'gaps')
    
    print('=======================================')
    a = GlobalAlignment.align('', 'CGACAAA', SCORES)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(SCORES))
    print(a.score_extremes(SCORES, gaps=True), 'gaps')
    
    print('=======================================')
    a = GlobalAlignment.align('AACGTACGTACAATAGTTTTACGATAACCGATAGCGATACCCATTAGACTATA', 'AAATAACGTAACTACAATAGTTCTACGATAACCGATATTGCGTTACCCAATAGACGATA', SCORES, print_matrices=False)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(SCORES))
    print(a.score_extremes(SCORES, gaps=True), 'gaps')
    
    print('=======================================')
    a = GlobalAlignment.align('ACAATTACC', 'AATCCG', SCORES)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(SCORES))
    print(a.score_extremes(SCORES, gaps=True), 'gaps')
    
    print('=======================================')
    a = GlobalAlignment.align('ACGGTTGC', 'AGCGTC', SCORES)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(SCORES))
    print(a.score_extremes(SCORES, gaps=True), 'gaps')
    
    print('=======================================')
    scores = build_score(match=10, transition=-30, transversion=-30, gap_open=-40, gap_extend=-1)
    a = GlobalAlignment.align('AAATTTGC', 'CGCCTTAC', scores)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(scores))
    print(a.score_extremes(scores, gaps=True), 'gaps')
    # AAATTTG~----C score=-98 <-- This is what it returns
    #      :|     |
    # ~----CGCCTTAC
    
    # ~-----AAATTTGC score=-70 <-- This is the optimal alignment
    #       |      |
    # CGCCTTA~-----C
    # TODO: Discover why this global alignment returns a non-optimal alignment, and fix it
    
    print('=======================================')
    # ACGGTTAACCA  matches    -> +9  final score -> +3
    # ||||  |||||  gap open   -> -4
    # ACGG~-AACCA  gap extend -> -2
    a = GlobalAlignment.align('ACGGTTAACCA', 'ACGGAACCA', SCORES)
    print(a.align1)
    print(a.align2)
    print(a.bitscore)
    print(a.score_extremes(SCORES))
    print(a.score_extremes(SCORES, gaps=True), 'gaps')

if (__name__ == "__main__"):
    test()
