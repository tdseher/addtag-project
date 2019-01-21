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
        
        spacer_sequences, spacer_sense_cuts, spacer_antisense_cuts = self.parse_motif_helper(spacer_motif)
        pam_sequences, pam_sense_cuts, pam_antisense_cuts = self.parse_motif_helper(pam_motif)
        
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
    # Two Cp1-family proteins, AsCpf1 (from Acidaminococcus)
    # and LbCpf1 (from Lachnospiraceae), have been shown to perform efficient
    # genome editing in human cells.
    #
    # Why use Cpf1 over Cas9?
    #  see https://benchling.com/pub/cpf1
    
    def parse_motif_helper(self, submotif):
        """
        Helper function that parses either the SPACER or PAM motif.
        Decodes quantifiers and returns a list
        """
        # Keep track of expanded sequences
        sequences = ['']
        
        # Keep track if a quantifier is being parsed
        quantifier = None
        
        # Iterate through the characters
        for i, c in enumerate(submotif):
            if (c == '{'): # Start the quantifier
                quantifier = ''
            elif (c == '}'): # End the quantifier
                quantifier_list = quantifier.split(',')
                if (len(quantifier_list) == 1):
                    min_length = int(quantifier)
                    max_length = min_length
                    
                elif (len(quantifier_list) == 2):
                    if (quantifier_list[0] == ''):
                        min_length = 0
                    else:
                        min_length = int(quantifier_list[0])
                    if (quantifier_list[1] == ''):
                        raise Exception("Motif quantifier '{" + quantifier + "}' contains no maximum value")
                    else:
                        max_length = int(quantifier_list[1])
                    if (min_length > max_length):
                        raise Exception("Motif quantifier '{" + quantifier + "}' minimum and maximum lengths are invalid")
                
                last_chars = [ x[-1] for x in sequences ]
                
                sequences = [ x[:-1] for x in sequences ]
                
                new_sequences = []
                for j, s in enumerate(sequences):
                    for length in range(min_length, max_length+1):
                        new_sequences.append(s + last_chars[j]*length)
                sequences = new_sequences
                quantifier = None
            elif (quantifier != None): # add current character to quantifier if it is open
                # We strip out any cut sites
                #if (c not in ['/', '\\', '|']):
                quantifier += c
            else: # add the current character to the expanded sequences
                for j in range(len(sequences)):
                    sequences[j] = sequences[j] + c
        
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
