#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/donors.py

# Import standard packages
import logging

# Import non-standard packages
import regex

# Import included AddTag-specific modules
from . import nucleotides





class Node(object):
    def __init__(self, data=None, parent=None, quant=(1,1)):
        self.data = data
        self.parent = parent
        self.quant = quant
    
    def __repr__(self):
        if (self.quant[0] == self.quant[1]):
            n = self.quant[0]
        else:
            n = self.quant
        return self.__class__.__name__ + '({}, n={})'.format(self.data, n)

class Seq(Node):
    def listify(self):
        return [self.data*n for n in range(self.quant[0], self.quant[1]+1)]

class And(Node):
    def __init__(self, data=None, parent=None, quant=(1,1)):
        super().__init__(data=data, parent=parent, quant=quant)
        if (data == None):
            self.data = []
    
    def add(self, d):
        self.data.append(d)
    
    def andlast(self):
        a = And(parent=self, data=[self.data[-1]])
        self.data[-1].parent = a
        self.data[-1] = a
    
    def listify(self):
        if ((self.data != None) and (len(self.data) > 0)):
            sequences = [''] # ['a', 'b', 'c']
            for d in self.data:
                d_list = d.listify()
                
                # Copy all existing sequences
                if (len(d_list) > 0):
                    new_sequences = [seq for seq in sequences for i in range(len(d_list))]
                else:
                    new_sequences = sequences
                
                # Join the new string to each duplicated sequence
                for i in range(len(sequences)):
                    for j, d in enumerate(d_list):
                        new_sequences[i*len(d_list)+j] += d
                
                sequences = new_sequences
            
            # Remove duplicate values
            return list(set([seq*n for seq in sequences for n in range(self.quant[0],self.quant[1]+1)]))
        else:
            return []
    
class Or(Node):
    def __init__(self, data=None, parent=None, quant=(1,1)):
        super().__init__(data=data, parent=parent, quant=quant)
        if (data == None):
            self.data = []
    
    def add(self, d):
        self.data.append(d)
    
    def andlast(self):
        a = And(parent=self, data=[self.data[-1]])
        self.data[-1].parent = a
        self.data[-1] = a
    
    def listify(self):
        #return [[d.listify() for d in self.data]*n for n in range(self.quant[0], self.quant[1]+1)]
        #return [seq*n for seq in sequences for n in range(self.quant[0],self.quant[1]+1)]
        sequences = []
        for d in self.data:
            d_list = d.listify()
            for l in d_list:
                for i in range(self.quant[0], self.quant[1]+1):
                    sequences.append(l*i)
        # Remove duplicate values
        return list(set(sequences))

def parse_submotif(submotif):
    n = And()
    root = n
    or_next = False
    
    # If the actual class is important, then they can be stored in the capturing groups with this regex: (\d+)|(\w+)|([()])|([{}])|(,)
    #token_list = regex.findall(r'\w+|[()]|[{}]|,', submotif)
    token_list = regex.findall(r'[acgtrymkwsbdhvnACGTRYMKWSBDHVN./|\\]+|[()]|\{(?:\d+|\d*,\d+)\}|,', submotif)
    #token_list = regex.findall(r'[acgtrymkwsbdhvnACGTRYMKWSBDHVN./|\\]+|\d+|[()]|[{}]|,', submotif)
    #token_list = regex.findall(r'\w+(?!\{)|\w|[()]|[{}]|,', 'AT(NRG{1,2},(NNG){2,3},GAA,GAT,CAA)') # Makes a separate token for 'NR' and 'G'
    for token in token_list:
        if (token == '('):
            if not or_next and (len(n.data) > 0):
                if isinstance(n.data[-1], Seq):
                    n.andlast()
                    n = n.data[-1]
                elif isinstance(n.data[-1], And):
                    n = n.data[-1]
            n.add(Or(parent=n))
            n = n.data[-1]
        elif (token == ')'):
            while (not isinstance(n, Or)):
                n = n.parent
            n = n.parent
        elif token.startswith('{'):
            quant = token[1:-1].split(',')
            if (len(quant) == 1):
                quant = (int(quant[0]), int(quant[0]))
            elif (len(quant) == 2):
                if (quant[0] == ''):
                    quant = (0, int(quant[1]))
                elif (quant[1] == ''):
                    raise Exception("Motif quantifier '" + token + "' contains no maximum value")
                else:
                    quant = tuple(int(x) for x in quant)
                if (quant[0] > quant[1]):
                    raise Exception("Motif quantifier '" + token + "' minimum and maximum lengths are invalid")
            
            if (isinstance(n.data[-1], Seq) and (len(n.data[-1].data) > 1)):
                if isinstance(n, And):
                    d_first, d_last = n.data[-1].data[:-1], n.data[-1].data[-1:]
                    n.data[-1].data =d_first
                    n.add(Seq(d_last, parent=n, quant=quant))
                else: # elif isinstance(n, Or):
                    a = And(parent=n)
                    a.add(Seq(data=n.data[-1].data[:-1], parent=a))
                    a.add(Seq(data=n.data[-1].data[-1:], parent=a, quant=quant))
                    n.data[-1] = a
            else:
                n.data[-1].quant = quant
        elif (token == ','):
            if isinstance(n, And):
                while (not isinstance(n, Or)):
                    n = n.parent
            #or_next = True
        else:
            if (or_next or (len(n.data) == 0)):
                n.add(Seq(token, parent=n))
            else:
                if isinstance(n.data[-1], Seq):
                    n.andlast()
                    n = n.data[-1]
                elif isinstance(n.data[-1], And):
                    n = n.data[-1]
                n.add(Seq(token, parent=n))
        if (token == ','):
            or_next = True
        else:
            or_next = False
        #print(token, root, or_next)
    
    #return n.listify()
    #return n
    return root.listify()

class Motif(object):
    motifs = []
    def __init__(self, motif_string):
        # Save the input string: 'N{3,5}>NRG'
        self.motif_string = motif_string
        
        # Create [spacers, pams, side] list: (['NNN', 'NNNN', 'NNNNN'], ['NRG'], '>')
        self.parsed_list, self.spacer_sense_cuts, self.spacer_antisense_cuts, self.pam_sense_cuts, self.pam_antisense_cuts = self.parse_motif(self.motif_string)
        
        # Build the regex string: '((?:[ACGT][ACGT][ACGT])|(?:[ACGT][ACGT][ACGT][ACGT])|(?:[ACGT][ACGT][ACGT][ACGT][ACGT]))((?:[ACGT][AG]G))'
        self.regex_string = self.build_motif_regex(self.parsed_list[0], self.parsed_list[1], self.parsed_list[2])
        
        # Compile the regex
        self.compiled_regex = self.compile_regex(self.regex_string, ignore_case=True)
        
        # Save object to publicly-accessible, class list
        self.motifs.append(self)
        
        # Example:
        # motif          TTTN<N{19}/.{4}\
        # genome      ...CCAAAGCCAACGACTTTAGCTAGCTAAAGGACCTATGCCCATTACATGCCGCCAA...
        # sequence                     TTTAGCTAGCTAAAGGACCTATGCCCA
        # upstream           AGCCAACGAC
        # pam                          TTTA
        # target                           GCTAGCTAAAGGACCTATG
        # cut sites                                          ><  ><
        # sense cuts                                         ><
        # antisense cuts                                         ><
        # double-strand cuts
        # downstream                                              TTACATGCCG
        
    def parse_motif(self, motif):
        """
        Takes string of SPACER>PAM and returns ([spacer, spacer, ...], [pam, pam ...], '>')
        
        Example inputs: 'G{,2}N{19,20}>NGG', 'N{20,21}>NNGRRT', 'TTTN<N{20,23}'
        """
        # TODO: Add a linked quantifier that uses brackets '[' and ']'. For example
        #       N[12,13]|N[8,7]>NGG    would match   NNNNNNNNNNNN|NNNNNNNN>NGG
        #                              and           NNNNNNNNNNNNN|NNNNNNN>NGG
        #       Alternatiely, just rely on whether-or-not the left quantifier is smaller than the right quantifier?
        #       N{12,14}|N{8,6}>NGG
        
        # Possible quantifier delimiters
        #   N┤6,7├
        #   N╣6,7╠
        #   N╡6,7╞
        #   N╢6,7╟
        #   N<6,7>
        #   N{6,7}
        #   N(6,7)
        #   N[6,7]
        #   Nq6,7p
        #   N‹6,7›
        #   N«6,7»
        #   N‘6,7’
        #   N“6,7”
        #   N←6,7→
        #   N﴾6,7﴿
        
        # Possible cut characters
        #   sense, antisense, both
        #   /\|   '\' can be an escape
        #   ^v|   'v' is reserved
        #   *.|   '.' is reserved
        #   !¡|   '¡' is hard to type
        #   ?¿|   '¿' is hard to type
        #   ↑↓↕   best, but hard to type
        #   ╛╕╡   difficult
        #   ╧╤╪   difficult
        
        
        # TODO: Consider implementing a way of parsing nickase RGNs (only cut on one strand)
        # TODO: Consider implementing a way of parsing paired nickase RGNs
        # Still need to implement:
        #  Identify gRNA pairs on opposite strands appropriate for "nickase" Cas9
        #  allow for PAM-out and PAM-in orientations
        #  --motifs, where '.' is any non-spacer nucleotide (from the "Dimeric CRISPR RNA-guided FokI nucleases for highly specific genome editing" paper)
        #    'CCN<N{20}.{13,18}N{20}>NGG'     # PAM-out
        #    'N{20}>NGG.{13,18}CCN<N{20}'     # PAM-in
        #
        # The FokI-dCas9 restriction site can cut anywhere between either PAM
        # site (for PAM-out). It typically cuts with a 3' overhang
        # (Targeting individual subunits of the FokI restriction endonuclease to specific DNA strands):
        #   5'-.../.... ...-3'
        #   3'-... ....\...-5'
        #
        # Maybe this? It isn't perfect, but it lets the user specify the cut sites.
        #   'CCN<N{20}.{5,9}/....\.{9,5}N{20}>NGG' # One quantifier counts up while the other counts down?
        
        gt_count = motif.count('>')
        lt_count = motif.count('<')
        lb_count = motif.count('{')
        rb_count = motif.count('}')
        # Make sure motif does not violate basic rules
        if (gt_count + lt_count < 1):
            raise Exception("Motif lacks distinction between spacer and PAM sequences ('>' or '<' character)")
        elif (gt_count + lt_count > 1):
            raise Exception("Motif has too many '>' or '<' characters")
        if (lb_count != rb_count):
            raise Exception("Motif braces '{' and '}' do not match")
        if (motif.count(' ') > 0):
            raise Exception("Motif contains invalid space ' ' characters")
        if (motif.count('{}') > 0):
            raise Exception("Motif contains invalid quantifier '{}'")
        if (motif.count('{,}') > 0):
            raise Exception("Motif contains invalid quantifier '{,}'")
        
        if (gt_count == 1):
            spacer_motif, pam_motif = motif.split('>')
            side = '>'
        elif (lt_count == 1):
            pam_motif, spacer_motif = motif.split('<')
            side = '<'
        
        spacer_sequences, spacer_sense_cuts, spacer_antisense_cuts = self.parse_motif_cuts(spacer_motif)
        pam_sequences, pam_sense_cuts, pam_antisense_cuts = self.parse_motif_cuts(pam_motif)
        
        return (spacer_sequences, pam_sequences, side), spacer_sense_cuts, spacer_antisense_cuts, pam_sense_cuts, pam_antisense_cuts
    
    # List of common motifs obtained from
    #  https://github.com/maximilianh/crisporWebsite/crispor.py
    
    # from (https://benchling.com/pub/cpf1):
    # Cpf1 is an RNA-guided nuclease, similar to Cas9. It recognizes a T-rich
    # PAM, TTTN, but on the 5' side of the guide. This makes it distinct from
    # Cas9, which uses an NGG PAM on the 3' side. The cut Cpf1 makes is
    # staggered. In AsCpf1 and LbCpf1, it occurs 19 bp after the PAM on the
    # targeted (+) strand and 23 bp on the other strand, as shown here:
    #                        Cpf1
    #   TTTC GAGAAGTCATCTAATAAGG|CCAC TGTTA
    #   AAAG CTCTTCAGTAGATTATTCC GGTG|ACAAT
    #   =PAM =========gRNA====== =
    #
    # Benchling suggests 20 nt guides can be used for now, as there is not real
    # suggestion for optimal guide length. Robust guide scores for Cpf1 are
    # still in development, but simple scoring based on the number of off-target
    # sites is available on Benchling.
    # 
    # Cpf1 requires only a crRNA for activity and does not need a tracrRNA to
    # also be present.
    #
    # Two Cpf1-family proteins, AsCpf1 (from Acidaminococcus)
    # and LbCpf1 (from Lachnospiraceae), have been shown to perform efficient
    # genome editing in human cells.
    #
    # Why use Cpf1 over Cas9?
    #  see https://benchling.com/pub/cpf1
    
    #####################
    ### New code here ###
    #####################
    
    def parse_motif_cuts(self, submotif):
        """
        Helper function that parses either the SPACER or PAM motif.
        Decodes quantifiers and returns a list
        """
        # Expand 'submotif' into a list of sequences and sort it
        sequences = sorted(parse_submotif(submotif))
        
        # Define empty list holding all cut sites
        sense_cuts = []
        antisense_cuts = []
        
        # Iterate through sequences, and record the location of any cut sites
        # Then strip out any cut sites
        for i, seq in enumerate(sequences):
            sense_cuts.append([])
            antisense_cuts.append([])
            found = 0
            
            nseq = seq
            for j, c in enumerate(seq):
                if (c  == '/'):
                    sense_cuts[-1].append(j-found)
                    nseq = nseq.replace('/', '', 1)
                    found += 1
                elif (c == '\\'):
                    antisense_cuts[-1].append(j-found)
                    nseq = nseq.replace('\\', '', 1)
                    found += 1
                elif (c == '|'):
                    sense_cuts[-1].append(j-found)
                    antisense_cuts[-1].append(j-found)
                    nseq = nseq.replace('|', '', 1)
                    found += 1
            sequences[i] = nseq
        
        return sequences, sense_cuts, antisense_cuts
    
    def compile_regex(self, re_pattern, ignore_case=False, enhance_match=False, best_match=False):
        """
        Returns compiled regex
        """
        # Choose the regex flags
        myflags = 0
        if ignore_case:
            myflags |= regex.IGNORECASE
        if enhance_match:
            myflags |= regex.ENHANCEMATCH
        if best_match:
            myflags |= regex.BESTMATCH
            
        c = regex.compile(re_pattern, flags=myflags) # regex.ENHANCEMATCH|regex.IGNORECASE
        return c
    
    def build_motif_regex(self, spacers, pams, side, anchored=False):
        """
        Returns regex pattern string for motif
        """
        spacer_pattern = '(' + '|'.join([self.build_regex_pattern(x, capture=False) for x in spacers]) + ')'
        pam_pattern = '(' + '|'.join([self.build_regex_pattern(x, capture=False) for x in pams]) + ')'
        if (side == '>'): # SPACER>PAM
            if anchored:
                re_pattern = '^'+ spacer_pattern + pam_pattern + '$'
            else:
                re_pattern = spacer_pattern + pam_pattern
        elif (side == '<'): # PAM<SPACER
            if anchored:
                re_pattern = '^'+ pam_pattern + spacer_pattern + '$'
            else:
                re_pattern = pam_pattern + spacer_pattern
        return re_pattern
    
    def build_regex_pattern(self, iupac_sequence, max_substitutions=0, max_insertions=0, max_deletions=0, max_errors=0, capture=True):
        """
        Build a regular expression pattern for the nucleotide search, taking IUPAC
        ambiguities into account.
        
        So far for DNA only
        Expects string as input (NOT a list)
        """
        # Format the fuzzy matching restrictions
        fuzzy = []
        if (max_substitutions > 0):
            fuzzy.append('s<={}'.format(max_substitutions))
        if (max_insertions > 0):
            fuzzy.append('i<={}'.format(max_insertions))
        if (max_deletions > 0):
            fuzzy.append('d<={}'.format(max_deletions))
        if (max_errors > 0):
            fuzzy.append('e<={}'.format(max_errors))
        if (len(fuzzy) > 0):
            fuzzy = '{' + ','.join(fuzzy) + '}'
        else:
            fuzzy = ''
        
        # Convert the IUPAC sequence to an equivalent regex
        iupac = {
            'a': 'a',
            'c': 'c',
            'g': 'g',
            't': 't',
            'r': '[ag]',
            'y': '[ct]',
            'm': '[ac]',
            'k': '[gt]',
            'w': '[at]',
            's': '[cg]',
            'b': '[cgt]',
            'd': '[agt]',
            'h': '[act]',
            'v': '[acg]',
            'n': '[acgt]',
            
            'A': 'A',
            'C': 'C',
            'G': 'G',
            'T': 'T',
            'R': '[AG]',
            'Y': '[CT]',
            'M': '[AC]',
            'K': '[GT]',
            'W': '[AT]',
            'S': '[CG]',
            'B': '[CGT]',
            'D': '[AGT]',
            'H': '[ACT]',
            'V': '[ACG]',
            'N': '[ACGT]',
            
            '.': '.', # Any character
            
            '|': r'\|', # Double-stranded cut
            '/': r'\/', # Sense cut
            '\\': r'\\', # anti-sense cut
        }
        sequence = ''.join(map(lambda x: iupac[x], iupac_sequence))
        if capture:
            pattern = '(' + sequence + ')' + fuzzy
        else:
            pattern = '(?:' + sequence + ')' + fuzzy
        #logger.info('Built regex string: {!r}'.format(pattern))
        return pattern
    
    def generate_sequence(self, compositions=None, complement=False, samples=100):
        """
        Generate a random sequence according to input sequence composition
        or with uniform composition of None.
        """
        
        if compositions:
            seqs = []
            scores = []
            for i in range(samples):
                s = bartag.random_motif_sequence(self.motif_string)
                seqs.append(s)
                scores.append(nucleotides.sequence_likelihood(s, compositions)/len(s))
            
            #sorted_seqs = sorted(seqs, key=lambda x: seqs[x], reverse=not complement)
            #return sorted_seqs[0]
            if complement:
                return min(zip(scores, seqs))[1]
            else:
                return max(zip(scores, seqs))[1]
        else:
            return bartag.random_motif_sequence(self.motif_string)

class OnTargetMotif(Motif):
    motifs = []

class OffTargetMotif(Motif):
    motifs = []
