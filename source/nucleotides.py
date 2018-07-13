#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/nucleotides.py

# List general Python imports
import sys
import random
import difflib
import logging

# import non-standard package
import regex

logger = logging.getLogger(__name__)

class SlidingWindow(object):
    """Iterator for a sliding window across a nucleotide sequence"""
    def __init__(self, sequence, window=20, start=0, stop=None, step=1):
        """Initialize a new instance of this iterator"""
        self.sequence = sequence
        self.window = window
        self.step = step
        self.start = start
        if (stop == None):
            self.stop = len(self.sequence)
        else:
            self.stop = stop
        self.position = self.start
    def __iter__(self):
        """Return this object as an iterator"""
        return self
    def __next__(self):
        """Return the next element if there is one, otherwise, raise a
        StopIteration exception
        """
        s = self.position
        e = self.position + self.window
        #if ((self.position >= self.stop) or (end > len(self.sequence))):
        if (e > self.stop):
            raise StopIteration
        seq = self.sequence[s:e]
        self.position += self.step
        return s, e, seq

def random_sequence(length=100):
    return ''.join(random.choice(['A', 'C', 'G', 'T']) for k in range(length))

def rc(seq, kind="dna"):
    """Returns the reverse-complement of a string containing DNA or RNA characters"""
    if (kind == "dna"):
        complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB') # exclude ws, WS
    elif (kind == "rna"):
        complements = str.maketrans('acgturymkbdhvACGTURYMKBDHV', 'ugcaayrkmvhdbTGCAAYRKMVHDB') # exclude ws, WS
    else:
        raise ValueError("'" + str(kind) + "' is an invalid argument for rc()")
    return seq.translate(complements)[::-1]

def ngg():
    '''
    Generator that returns NGG
    '''
    for m in ['A', 'C', 'G', 'T']:
        yield m + 'GG'

def kmers(k, y=''):
    '''
    Generator for recursive definition of k-mers in string form.
    k >= 0
    y for internal use
    '''
    if k==0:
        yield y
    else:
        for m in ['A', 'C', 'G', 'T']:
            yield from kmers(k-1, m+y)

def permute_genotypes(genes, y=''):
    """
    Generator for recursive definition of genotypes in string form.
    Example usage:
      >>> list(permute_genotypes(['aA', 'bB']))
      ['ab', 'aB', 'Ab', 'AB']
      >>> list(permute_genotypes(['ACGT']*2))
      ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    """
    if (len(genes) == 0):
        yield y
    else:
        for allele in genes[0]:
            # Code for Python 3
            yield from permute_genotypes(genes[1:], y+allele)
            
            # Code for Python 2 (not tested for efficiency)
            # _iter = iter(permute_genotypes(genes[1:], y+allele))
            # try:
            #     while True: #broken by StopIteration
            #         yield next(_iter)
            # except StopIteration as e:
            #     if e.args:
            #         yield e.args[0]

def lcs(string1, string2):
    """Find the longest common substring between two strings"""
    matcher = difflib.SequenceMatcher(None, string1, string2, True)
    match = matcher.find_longest_match(0, len(string1), 0, len(string2))
    # Match(a=0, b=15, size=9)
    return match
    # print(string1[match.a: match.a + match.size])
    # print(string2[match.b: match.b + match.size])

def ridentities(seq1, seq2):
    """
    Returns the number of invariant nucleotides from the 3' side
    of the alignment. Assumes alignments are right-justified.
    
    input:
      seq1 = 'GCTATCTGGGACCAGGCTGA'
      seq2 =  'CTATCTGCGACCAGGCTGA'
    output:
      11
    """
    identities = 0
    min_len = min(len(seq1), len(seq2))
    for i in range(-1,-min_len-1, -1):
        if (seq1[i] != seq2[i]):
            break
        else:
            identities += 1
    return identities

def lidentities(seq1, seq2):
    """
    Returns the number of invariant nucleotides from the 5' side
    of the alignment. Assumes alignments are left-justified.
    
    input:
      seq1 = 'GCTATCTGGGACCAGGCTGA'
      seq2 = 'GCTATCTGCGACCAGGCTG'
    output:
      8
    """
    identities = 0
    min_len = min(len(seq1), len(seq2))
    for i in range(min_len):
        if (seq1[i] != seq2[i]):
            break
        else:
            identities += 1
    return identities

def filter_polyt(sequences, max_allowed=4):
    """Uses filter to remove sequences with consecutive Ts (in genome) or
    consecutive Us (in gRNA)
    
    subsequence   consecutive Ts   polymerase termination frequency
         TTTTTT   6                mostly
          TTTTT   5                sometimes
           TTTT   4                rarely
    """
    # consider using a generator instead
    #def filter_polyt(sequences, max_allowed):
    #    for s in sequences:
    #        if not ('T'*(max_allowed+1) in s):
    #            yield s
    
    return list(filter(lambda x: 'T'*(max_allowed+1) not in x, sequences))

# So far for DNA only
def build_regex_pattern(iupac_sequence, max_substitutions=0, max_insertions=0, max_deletions=0, max_errors=0, capture=True):
    """
    Build a regular expression pattern for the nucleotide search, taking IUPAC
    ambiguities into account.
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
    }
    sequence = ''.join(map(lambda x: iupac[x], iupac_sequence))
    if capture:
        pattern = '(' + sequence + ')' + fuzzy
    else:
        pattern = '(?:' + sequence + ')' + fuzzy
    logger.info('Built regex string: {!r}'.format(pattern))
    return pattern

# So far for DNA only
def build_regex(iupac_sequence, case_sensitive=False, max_substitutions=0, max_insertions=0, max_deletions=0, max_errors=0):
    """Build a regular expression for the nucleotide search, taking IUPAC
    ambiguities into account.
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
    
    # Choose the regex flags
    myflags = regex.ENHANCEMATCH | regex.IGNORECASE
    if case_sensitive:
        myflags = regex.ENHANCEMATCH
    #myflags = regex.IGNORECASE
    #if case_sensitive:
    #    myflags = 0
    
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
    }
    sequence = ''.join(map(lambda x: iupac[x], iupac_sequence))
    pattern = '(' + sequence + ')' + fuzzy
    logger.info('Compiled regex: {!r}'.format(pattern))
    compiled_regex = regex.compile(pattern, flags=myflags)
    return compiled_regex

def disambiguate_iupac(iupac_sequence, kind="dna"):
    """converts a string containing IUPAC nucleotide sequence to a list
    of non-iupac sequences.
    """
    if (kind == 'dna'):
        iupac = {
            'a': ['a'],
            'c': ['c'],
            'g': ['g'],
            't': ['t'],
            'r': ['a', 'g'],
            'y': ['c', 't'],
            'm': ['a', 'c'],
            'k': ['g', 't'],
            'w': ['a', 't'],
            's': ['c', 'g'],
            'b': ['c', 'g', 't'],
            'd': ['a', 'g', 't'],
            'h': ['a', 'c', 't'],
            'v': ['a', 'c', 'g'],
            'n': ['a', 'c', 'g', 't'],
            
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
    elif (kind == 'rna'):
        iupac = {
            'a': ['a'],
            'c': ['c'],
            'g': ['g'],
            't': ['u'],
            'u': ['u'],
            'r': ['a', 'g'],
            'y': ['c', 'u'],
            'm': ['a', 'c'],
            'k': ['g', 'u'],
            'w': ['a', 'u'],
            's': ['c', 'g'],
            'b': ['c', 'g', 'u'],
            'd': ['a', 'g', 'u'],
            'h': ['a', 'c', 'u'],
            'v': ['a', 'c', 'g'],
            'n': ['a', 'c', 'g', 'u'],
            
            'A': ['A'],
            'C': ['C'],
            'G': ['G'],
            'T': ['U'],
            'U': ['U'],
            'R': ['A', 'G'],
            'Y': ['C', 'U'],
            'M': ['A', 'C'],
            'K': ['G', 'U'],
            'W': ['A', 'U'],
            'S': ['C', 'G'],
            'B': ['C', 'G', 'U'],
            'D': ['A', 'G', 'U'],
            'H': ['A', 'C', 'U'],
            'V': ['A', 'C', 'G'],
            'N': ['A', 'C', 'G', 'U'],
        }
    sequences = ['']
    for char in iupac_sequence:
        # If the current character is an ambiguous base,
        # then copy all existing sequences
        for i in range(len(sequences)):
            for j in range(len(iupac[char]) - 1):
                sequences.append(sequences[i])
        # Sort the sequences... not very efficient
        sequences = sorted(sequences)
        # Add to each sequence a new, disambiguated nt
        for i in range(len(sequences)):
            sequences[i] += iupac[char][i % len(iupac[char])]
    
    # Throw an error if there is a non-canonical nucleotide
    assert all(map(lambda seq: all(map(lambda x: x in 'AaCcGgTtUu', seq)), sequences)) == True
    
    # Full list of non-canonical nucleotides/nucleobases can be found in:
    #  Chawla et al (2015) (https://doi.org/10.1093/nar/gkv606)
    #  .     low quality N, exotic base, or missing bp
    #  -     gap
    #  x, X  rarely used alternative for N, or exotic base
    
    return sequences

def count_errors(seq1, seq2):
    """Counts number of substitutions, insertions, and deletions"""
    # There is a bug in this function that makes it inaccurate. I will need
    # to replace it with a custom solution
    m = regex.match('(?:'+seq1+'){e}', seq2, flags=regex.ENHANCEMATCH|regex.IGNORECASE)
    # substitutions, insertions, deletions
    return m.fuzzy_counts

def split_target_sequence2(seq, pams):
    """
    Searches for all PAMs in sequence, and splits it. Assumes PAM at 3' end.
    Splits at the shortest PAM...
    Returns gRNA, PAM
    """
    # Build a regex to only match strings with PAM sites specified in args.pams
    re_pattern = '^(.*)(?:' + '|'.join(map(lambda x: build_regex_pattern(x), pams)) + ')$'
    #m = regex.match('^(.*)([ACGT]GG)$', nt)
    logger.info(re_pattern)
    m = regex.match(re_pattern, seq)
    if m:
        return m.group(1), m.group(2)
    else:
        # Force finding the PAM...
        #return seq, ''
        l = max(map(len, pams))
        return seq[:-l], seq[-l:]

def split_target_sequence(seq, pams, force=False, side='>'):
    """
    Searches for all PAMs in sequence, and splits it accordingly.
    Side specifies where to expect the PAM.
    Returns the longest matched PAM sequence.
    Returns: gRNA, PAM
    """
    # Build a regex to only match strings with PAM sites specified in args.parsed_motifs
    # by default, regex finds the longest pattern
    if (side == '>'): # SPACER>PAM
        re_pattern = '|'.join(map(lambda x: build_regex_pattern(x)+'$', pams))
        m = regex.search(re_pattern, seq, flags=regex.ENHANCEMATCH|regex.IGNORECASE)
        if m:
            return seq[:m.start()], seq[m.start():] # gRNA, PAM
        elif force: # Force finding the PAM
            l = max(map(len, pams))
            return seq[:-l], seq[-l:] # gRNA, PAM
        else:
            return None
        
    elif (side == '<'): # PAM<SPACER
        re_pattern = '|'.join(map(lambda x: '^'+build_regex_pattern(x), pams))
        m = regex.search(re_pattern, seq, flags=regex.ENHANCEMATCH|regex.IGNORECASE)
        if m:
            return seq[:m.end()], seq[m.end():] # gRNA, PAM
        elif force: # Force finding the PAM
            l = max(map(len, pams))
            return seq[:l], seq[l:] # gRNA, PAM
        else:
            return None
    
    else:
        return None

def motif_search2(sequence, side, compiled_regex):
    """
    Return list of all (target, start, end, spacer, PAM) found within sequence
    """
    
    return_matches = []
    
    matches = compiled_regex.finditer(sequence, overlapped=True)
    
    if (side == '>'): # SPACER>PAM
        for m in matches:
            seq = m.group()
            start, end = m.span()
            spacer, pam = m.groups()
            return_matches.append((seq, start, end, spacer, pam))
    
    elif (side == '<'): # PAM<SPACER
        for m in matches:
            seq = m.group()
            start, end = m.span()
            spacer, pam = m.groups()[::-1] # reverse the order
            return_matches.append((seq, start, end, spacer, pam))
    
    return return_matches

def motif_search(sequence, spacers, pams, side):
    """
    Return list of all (target, start, end, spacer, PAM) found within sequence
    """
    spacer_pattern = '(' + '|'.join([build_regex_pattern(x, capture=False) for x in spacers]) + ')'
    pam_pattern = '(' + '|'.join([build_regex_pattern(x, capture=False) for x in pams]) + ')'
    
    return_matches = []
    
    if (side == '>'): # SPACER>PAM
        re_pattern = spacer_pattern + pam_pattern
        matches = regex.finditer(re_pattern, sequence, flags=regex.ENHANCEMATCH|regex.IGNORECASE, overlapped=True)
        for m in matches:
            seq = m.group()
            start, end = m.span()
            spacer, pam = m.groups()
            return_matches.append((seq, start, end, spacer, pam))
    
    elif (side == '<'): # PAM<SPACER
        re_pattern = pam_pattern + spacer_pattern
        matches = regex.finditer(re_pattern, sequence, flags=regex.ENHANCEMATCH|regex.IGNORECASE, overlapped=True)
        for m in matches:
            seq = m.group()
            start, end = m.span()
            spacer, pam = m.groups()[::-1] # reverse the order
            return_matches.append((seq, start, end, spacer, pam))
    
    return return_matches

def motif_conformation2(sequence, side, compiled_regex):
    """
    Checks if the sequence matches at least one of the possible SPACER>PAM
    combinations.
    
    Returns (gRNA, PAM) if valid, otherwise None
    """
    m = compiled_regex.match(sequence)
    if m:
        if (side == '>'): # SPACER>PAM
            return m.groups() # (gRNA, PAM)
        elif (side == '<'): # PAM<SPACER
            return m.groups()[::-1] # (gRNA, PAM)
    else:
        return None

def motif_conformation(sequence, spacers, pams, side):
    """
    Checks if the sequence matches at least one of the possible SPACER>PAM
    combinations.
    
    Arguments:
     sequence - A nucleotide sequence that may or may not match the SPACER>PAM motif
      spacers - list of acceptable spacer motifs
         pams - list of acceptable pam motifs
         side - '>' or '<'
    
    Returns (gRNA, PAM) if valid, otherwise None
    """
    spacer_pattern = '(' + '|'.join([build_regex_pattern(x, capture=False) for x in spacers]) + ')'
    pam_pattern = '(' + '|'.join([build_regex_pattern(x, capture=False) for x in pams]) + ')'
    if (side == '>'): # SPACER>PAM
        re_pattern = '^'+ spacer_pattern + pam_pattern + '$'
        m = regex.match(re_pattern, sequence, flags=regex.ENHANCEMATCH|regex.IGNORECASE)
        if m:
            return m.groups() # (gRNA, PAM)
        else:
            return None
    elif (side == '<'): # PAM<SPACER
        re_pattern = '^'+ pam_pattern + spacer_pattern + '$'
        m = regex.match(re_pattern, sequence, flags=regex.ENHANCEMATCH|regex.IGNORECASE)
        if m:
            return m.groups()[::-1] # (gRNA, PAM)
        else:
            return None

def build_dDNA_regex_comp2(dDNA, minimum_length=2):
    """incomplete"""
    seqs = []
    l = len(dDNA)
    
    left_seq = ''
    left_count = 0
    left_iter = iter(dDNA)
    for li in range(minimum_length):
        try:
            left_count += 1
            left_seq += next(left_iter)
        except StopIteration:
            pass
    
    while(left_count <= l):
        right_seq = ''
        right_count = 0
        right_iter = iter(dDNA[::-1])
        for ri in range(minimum_length):
            try:
                right_count += 1
                right_seq = next(right_iter) + right_seq
            except StopIteration:
                pass
        
        while(right_count < l-left_count):
            print('left_count=', left_count, 'left_seq=', left_seq, 'right_count=', right_count, 'right_seq=', right_seq)
            try:
                right_count += 1
                right_seq = next(right_iter) + right_seq
                seqs.append((left_seq, right_seq))
            except StopIteration:
                pass
        
        try:
            left_count += 1
            left_seq += next(left_iter)
            
        except StopIteration:
            #left_count += 1
            seqs.append((dDNA, ''))
        
    return seqs

def build_dDNA_regex_comp1(dDNA):
    """Starts at length=0 for both left and right"""
    seqs = []
    l = len(dDNA)
    
    left_seq = ''
    left_count = 0
    left_iter = iter(dDNA)
    
    while(left_count <= l):
        right_seq = ''
        right_count = 0
        right_iter = iter(dDNA[::-1])
        
        while(right_count < l-left_count):
            #print('left_count=', left_count, 'left_seq=', left_seq, 'right_count=', right_count, 'right_seq=', right_seq)
            try:
                right_count += 1
                right_seq = next(right_iter) + right_seq
                seqs.append((left_seq, right_seq))
            except StopIteration:
                pass
        
        try:
            left_count += 1
            left_seq += next(left_iter)
            
        except StopIteration:
            #left_count += 1
            seqs.append((dDNA, ''))
        
    return seqs
        
def split_dDNA_flanks(dDNA):
    """Starts at length=1 at for both left and right"""
    seqs = []
    l = len(dDNA)
    
    left_seq = ''
    for i, left_nt in enumerate(dDNA):
        left_seq += left_nt
        
        right_seq = ''
        for j, right_nt in enumerate(dDNA[::-1]):
            if (l-j-1 > i):
                right_seq = right_nt + right_seq
                seqs.append((left_seq, right_seq))
    return seqs

def filter_dDNA_flanks(lists, minimum_length=40):
    seqs = []
    for left, right in lists:
        if ((len(left) >= minimum_length) and (len(right) >= minimum_length)):
            seqs.append((left, right))
    return seqs

def build_dDNA_regex(dDNA):
    splits = filter_dDNA_flanks(split_dDNA_flanks(dDNA))
    flags=regex.ENHANCEMATCH | regex.IGNORECASE
    
    pattern = '|'.join([left+'(.{,2000})'+right for left, right in splits]) # '(.*)' is too slow
    
    compiled_regex = regex.compile(pattern, flags=flags)
    return compiled_regex

def test():
    """Code to test the classes and functions in 'source/nucleotides.py'"""
    
    print("=== SlidingWindow ===")
    
    print("=== rc ===")
    
    print("=== filter_polyt ===")
    
    print("=== build_regex_pattern ===")
    
    print("=== build_regex ===")
    
    print("=== find_target_matches ===")
    
    print("=== disambiguate_iupac ===")
    
    print("=== count_errors ===")
    
    print("=== split_target_sequence ===")
    

if (__name__ == '__main__'):
    test()
