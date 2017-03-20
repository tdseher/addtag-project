#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/nucleotides.py

# List general Python imports
import sys
import difflib

# import non-standard package
import regex

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

def rc(seq, kind="dna"):
    """Returns the reverse-complement of a string containing DNA or RNA characters"""
    if (kind == "dna"):
        complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB') # exclude ws, WS
    elif (kind == "rna"):
        complements = str.maketrans('acgturymkbdhvACGTURYMKBDHV', 'ugcaayrkmvhdbTGCAAYRKMVHDB') # exclude ws, WS
    else:
        raise ValueError("'" + str(kind) + "' is an invalid argument for rc()")
    return seq.translate(complements)[::-1]

def lcs(string1, string2):
    """Find the longest common substring between two strings"""
    matcher = difflib.SequenceMatcher(None, string1, string2, True)
    match = matcher.find_longest_match(0, len(string1), 0, len(string2))
    # Match(a=0, b=15, size=9)
    return match
    # print(string1[match.a: match.a + match.size])
    # print(string2[match.b: match.b + match.size])

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
def build_regex_pattern(iupac_sequence, max_substitutions=0, max_insertions=0, max_deletions=0, max_errors=0):
    """Build a regular expression pattern for the nucleotide search, taking IUPAC
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
    pattern = '(' + sequence + ')' + fuzzy
    print('Built regex string: {!r}'.format(pattern), file=sys.stderr)
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
    print('Compiled regex: {!r}'.format(pattern), file=sys.stderr)
    compiled_regex = regex.compile(pattern, flags=myflags)
    return compiled_regex

def find_target_matches(compiled_regex, contigs, overlap=False):
    '''Finds all instances of the compiled_regex within the contigs(dict),
    and returns list of matches as 0-indexed
    [contig, start, stop, orientation, oriented-sequence]
    '''
    
    matches = []
    
    for contig in contigs:
        # '+' orientation
        ref = contigs[contig]
        ref_len = len(ref)
        for m in compiled_regex.finditer(ref, overlapped=overlap):
            s = m.start()
            e = m.end()
            matches.append([contig, s, e, '+', ref[s:e]])
        # '-' orientation
        rc_ref = rc(contigs[contig])
        rc_ref_len = len(rc_ref)
        for m in compiled_regex.finditer(rc_ref, overlapped=overlap):
            s = m.start()
            e = m.end()
            #print("match: ", s, e, m.groups(), rc_ref[s:e], ref[ref_len-e:ref_len-s], file=sys.stderr)
            assert m.groups()[0] == rc_ref[s:e]
            assert m.groups()[0] == rc(ref[ref_len-e:ref_len-s]) # m.groups()[0] should be what was matched...
            assert ref_len == rc_ref_len
            matches.append([contig, rc_ref_len-s, rc_ref_len-e, '-', rc_ref[s:e]])
    print("Matching finished", file=sys.stderr)
    
    return matches

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
    m = regex.match('(?:'+seq1+'){e}', seq2, flags=regex.ENHANCEMATCH|regex.IGNORECASE)
    # substitutions, insertions, deletions
    return m.fuzzy_counts

def split_target_sequence(seq, pams):
    """Searches for all PAMs in sequence, and splits accordingly
    Returns: gRNA, PAM
    """
    # Build a regex to only match strings with PAM sites specified in args.pams
    re_pattern = '|'.join(map(lambda x: build_regex_pattern(x)+'$', pams))
    
    # by default, regex finds the longest pattern
    m = regex.search(re_pattern, seq, flags=regex.ENHANCEMATCH|regex.IGNORECASE)
    if m:
        # gRNA, PAM
        return seq[:m.start()], seq[m.start():]
    else:
        # Force finding the PAM...
        #return seq, ''
        l = max(map(len, pams))
        return seq[:-l], seq[-l:]

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