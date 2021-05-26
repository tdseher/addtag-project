#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/nucleotides.py

# List general Python imports
import sys
import random
import difflib
import logging
from itertools import product

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

def kmer_distribution(kmer):
    """Returns a dictionary that is a kmer distribution, with all values set to 0.0"""
    dist = {}
    for i in range(kmer):
        dist = __kmer_distribution_helper(dist)
    return dist

def __kmer_distribution_helper(ddd):
    """Helper function for kmer_distribution"""
    new = {}
    bases = ['A','C','G','T']
    for i in bases:
        if (len(ddd) > 0):
            for d in ddd:
                new[d + i] = 0.0 # 0.0 for floats
        else:
            new[i] = 0.0 # 0.0 for floats
    return new

def make_mirrored_kmer_dist(kmer_length):
    return_dict = {}
    ddd = {}
    kmers = sorted(kmer_distribution(kmer_length))
    for k in kmers:
        rk = rc(k)
        if not ((k in ddd) and (rk in ddd)):
            ddd[k] = 1
            ddd[rk] = 2
            
            temp_list = [0]
            return_dict[k] = temp_list
            return_dict[rk] = temp_list
        #else:
        #   ddd[k] += 4
    return return_dict

def reduce_mirrored_kmer_dist(dist):
    f_added = {}
    for f in sorted(dist):
        if not (rc(f)+' '+f in f_added):
            f_added[f+' '+rc(f)] = dist[f][0]
    return f_added

def normalize_reduced_dist(dist):
    return_dist = {}
    # The following two measures are the same:
    #   print sum(dist.values())
    #   print (((len(seq) - kmer_size) / step_size) + 1)
    if (sum(dist.values()) > 0):
        for k in dist:
            #return_dist[k] = (1.0 * dist[k]) / (((len(seq) - kmer_size) / step_size) + 1)
            return_dist[k] = (1.0 * dist[k]) / sum(dist.values())
    else:
        return_dist = dist
    return return_dist

def normalize_dist(dist):
    return_dist = {}
    # The following two measures are the same:
    #   print sum(dist.values())
    #   print (((len(seq) - kmer_size) / step_size) + 1)
    dist_sum = sum(dist.values())
    if (sum(dist.values()) > 0):
        for k in dist:
            #return_dist[k] = (1.0 * dist[k]) / (((len(seq) - kmer_size) / step_size) + 1)
            return_dist[k] = float(dist[k]) / dist_sum
    else:
        return_dist = dist
    return return_dist

def get_seq_list_dist_stranded(seqs, step_size, kmer_size):
    kmers = kmer_distribution(kmer_size)
    for seq in seqs:
        pos = 0
        while (pos <= len(seq) - kmer_size):
            try:
                kmers[seq[pos:pos+kmer_size]] += 1
            except KeyError: # has iupac ambiguity code
                good_kmers = fix_kmer(seq[pos:pos+kmer_size])
                for i in good_kmers:
                    kmers[i] += 1.0 / len(good_kmers)
            pos += step_size
    return normalize_dist(kmers)

def get_seq_list_dist(seqs, step_size, kmer_size):
    kmers = make_mirrored_kmer_dist(kmer_size)
    for seq in seqs:
        pos = 0
        while (pos <= len(seq) - kmer_size):
            try:
                kmers[seq[pos:pos+kmer_size]][0] += 1
            except KeyError: # has iupac ambiguity code
                good_kmers = fix_kmer(seq[pos:pos+kmer_size])
                for i in good_kmers:
                    kmers[i][0] += 1.0 / len(good_kmers)
            pos += step_size
    return normalize_dist(reduce_mirrored_kmer_dist(kmers))

def get_seq_dist(seq, step_size, kmer_size):
    kmers = make_mirrored_kmer_dist(kmer_size)
    pos = 0
    while (pos <= len(seq) - kmer_size):
        try:
            kmers[seq[pos:pos+kmer_size]][0] += 1
        except KeyError: # has iupac ambiguity code
            good_kmers = fix_kmer(seq[pos:pos+kmer_size])
            for i in good_kmers:
                kmers[i][0] += 1.0 / len(good_kmers)
        pos += step_size
    return normalize_reduced_dist(reduce_mirrored_kmer_dist(kmers))


def fix_kmer(bad_kmer):
    # bad_kmer = ANCN
    iupac = {
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
    fixed_kmers = ['']
    for char in bad_kmer:
        if char in ['a','c','g','t','A','C','G','T']:
            for nc in range(len(fixed_kmers)):
                fixed_kmers[nc] += char
        else:
            # split
            for ool in range(len(fixed_kmers)):
                for num in range(len(iupac[char]) - 1):
                    fixed_kmers.append(fixed_kmers[ool])
            fixed_kmers = sorted(fixed_kmers)
            #print "split", fixed_kmers
            # add
            for pp in range(len(fixed_kmers)):
                fixed_kmers[pp] += iupac[char][pp % len(iupac[char])]
        #print "step", fixed_kmers
    return fixed_kmers

def complement_kmer_distribution(dist):
    """
    Input: dictionary with keys like this "ACAT ATGT" and values with floats.
    sum(dist.values()) should equal 1. This means that the input distribution is normalized
    """
    # If 0 is present, then shift all values by minimum
    cdist = {}
    if 0 in dist.values():
        minv = min([x for x in dist.values() if x > 0])
        for k,v in dist.items():
            cdist[k] = v+minv
    else:
        for k,v in dist.items():
            cdist[k] = v
    
    s = sum([1.0/x for x in cdist.values()])
    for k,v in cdist.items():
        cdist[k] = (1.0/v)/s
    
    return cdist

def sequence_likelihood(seq, dist, step_size=1):
    """
    Computes the likelihood of the sequence given the kmer distribution.
    Output likelihood highly contingent on sequence length
    """
    kmers = {}
    kmer_size = 0
    for k, v in dist.items():
        k1, k2 = k.split(' ')
        kmers[k1] = v
        kmers[k2] = v
        kmer_size = len(k1)
    
    score = 0
    
    # sliding window
    for i in range(0, len(seq)-kmer_size+1, step_size):
        score += kmers[seq[i:i+kmer_size]]
    
    return score
    

#def random_sequence(length=100):
#    return ''.join(random.choice(['A', 'C', 'G', 'T']) for k in range(length))

def random_sequence(length=100, compositions=None):
    '''
    Expects kmer=1
    if None defined, then
    compositions = {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25}
    '''
    if (compositions == None):
        compositions = {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25}
    return ''.join(random_choices(['A', 'C', 'G', 'T'], weights=[compositions['A'], compositions['C'], compositions['G'], compositions['T']], k=length))

def flip(text):
    """Rotates the input text 180 degrees"""
    a = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    f = 'ɐqɔpǝɟɓɥᴉſʞ┐ɯuodbɹsʇnʌʍxʎzⱯʚ'+'\u0186'+'pƎℲϑHIſʞ˥WNOԀÕᴚS┴ՈΛMXʎZ'
    return text.translate(str.maketrans(a, f))[::-1]

def random_choices(population, weights=None, k=1):
    """Return a k sized list of population elements chosen with replacement."""
    #if hasattr(random, "choices"):
    if ('choices' in random.__all__):
        return random.choices(population, weights, k)
    else:
        if (weights != None):
            weight_sum = sum(weights)
            choices = zip(population, weights)
            values = []
            for i in range(k):
                r = random.uniform(0, weight_sum)
                upto = 0
                for c, w in choices:
                    if upto + w >= r:
                        values.append(c)
                        break
                    upto += w
                else:
                    values.append(random.choice(population))
            return values
        else:
            return [random.choice(population) for i in range(k)]

# def rc(seq, kind="dna"):
#     """Returns the reverse-complement of a string containing DNA or RNA characters"""
#     if (kind == "dna"):
#         complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB') # exclude ws, WS
#     elif (kind == "rna"):
#         complements = str.maketrans('acgturymkbdhvACGTURYMKBDHV', 'ugcaayrkmvhdbTGCAAYRKMVHDB') # exclude ws, WS
#     else:
#         raise ValueError("'" + str(kind) + "' is an invalid argument for rc()")
#     return seq.translate(complements)[::-1]

DNA_COMPLEMENTS = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB') # exclude ws/WS, n/N
RNA_COMPLEMENTS = str.maketrans('acgturymkbdhvACGTURYMKBDHV', 'ugcaayrkmvhdbTGCAAYRKMVHDB') # exclude ws/WS, n/N

def rc(seq, kind="dna"):
    """Returns the reverse-complement of a string containing DNA or RNA characters"""
    if (kind == "dna"):
        complements = DNA_COMPLEMENTS
    elif (kind == "rna"):
        complements = RNA_COMPLEMENTS
    else:
        raise ValueError("'" + str(kind) + "' is an invalid argument for rc()")
    return seq.translate(complements)[::-1]

def permutations(length, words='ACGT'):
    return list(''.join(x) for x in product(words, repeat=length))

def random_permutations(length, words='ACGT', n=None):
    if (n >= len(words)**length):
        n = None
    if (n == None):
        return list(''.join(x) for x in product(words, repeat=length))
    #elif ((n < 0) or (n > len(words)**length)):
    #    raise Exception('n is greater than len(words)**length')
    elif (n < 0):
        raise Exception('n is < 0, but it should be >=0')
    else:
        d = {}
        while (len(d) < n):
            d[random_sequence(length)] = 1
        return list(d.keys())

def kmers(k):
    '''
    Make all possible canonical nucleotide strings of length k
    :param k: length of the sequence to return
    :return: list
    '''
    return list(''.join(x) for x in product('ACGT', repeat=k))

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

def lcs(string1, string2, autojunk=False):
    """Find the longest common substring between two strings"""
    matcher = difflib.SequenceMatcher(None, string1, string2, autojunk)
    match = matcher.find_longest_match(0, len(string1), 0, len(string2))
    # Match(a=0, b=15, size=9)
    return match
    # print(string1[match.a: match.a + match.size])
    # print(string2[match.b: match.b + match.size])

def pident(seq1, seq2, min_identity=0.8):
    """
    Returns the percent identity of two strings,
    relative to the longer one (global alignment),
    and permitting for up to 1-min_identity errors.
    
    DOES NOT disambiguate the shorter sequence.
    IUPAC ambiguities in the longer sequence count as matches,
    and ambiguities in the shorter sequence counts as mismatches.
    
    Also DOES NOT do perfect disambiguous calculations:
    If N aligns to R, this SHOULD be a score of 0.5, but
    This function will count it as 0.
    """
    # If lengths of inputs are tied, choose the one with ambiguities
    # as the pattern
    if (len(seq1) == len(seq2)):
        p1 = build_regex_pattern(seq1, capture=False)
        p2 = build_regex_pattern(seq2, capture=False)
        if (len(p1) < len(p2)):
            shorter, longer = seq1, seq2
        else:
            shorter, longer = seq2, seq1
    else:
        # Find which input is shorter, and which is longer
        shorter, longer = sorted([seq1, seq2], key=lambda x:len(x))
    
    # Calculate number of permissable errors
    e = int(len(longer)*max(0, min(1, 1-min_identity)))
    
    # Convert into IUPAC regex pattern
    pattern = build_regex_pattern(longer, max_errors=e, capture=False)
    m = regex.match(pattern, shorter, flags=regex.IGNORECASE|regex.BESTMATCH)
    if m:
        return (len(longer) - sum(m.fuzzy_counts))/len(longer)
    
    return 0

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
    # TODO: Change this to a 'zip(seq1, seq2)' Iterator to improve speed
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
    # TODO: Change this to a 'zip(seq1, seq2)' Iterator to improve speed
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

def get_nt_count(seq, case_sensitive=False):
    '''
    Ambiguity-aware canonical nt residue ('A', 'C', 'G', 'T') counter
    :param seq: string or list containing DNA sequence to count nts
    :return: dict
    '''
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
    
    output = {
        'A': 0.0,
        'C': 0.0,
        'G': 0.0,
        'T': 0.0,
    }
    if case_sensitive: 
        output['a'] = 0.0
        output['c'] = 0.0
        output['g'] = 0.0
        output['t'] = 0.0
    
    for c in seq:
        nts = iupac[c]
        for nt in nts:
            if case_sensitive:
                output[nt] += 1/len(nts)
            else:
                output[nt.upper()] += 1/len(nts)
    
    return output

def get_nt_freq(nts, seq):
    '''
    :param nts: str or list of nucleotides to count. For instance 'GC' or ['G', 'C']
    :param seq: sequence to count
    :return: float
    '''
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

    nt_set = set()
    for nt in nts:
        nt_set.update(iupac[nt])

    nt_count = 0

    for s in seq:
        nt_denom = iupac[s]
        nt_num = nt_set.intersection(nt_denom)
        nt_count += len(nt_num)/len(nt_denom)

    return nt_count/len(seq)

def gc_score(seq):
    """Caclulates the %GC content for seq. IUPAC ambiguities okay."""
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
    
    gc = 0.0
    for i in seq:
        k = iupac[i]
        for j in ['G', 'g', 'C', 'c']:
            if j in k:
                gc += 1.0/len(k)
    
    return 100*gc/len(seq)

def get_gc_freq(seq):
    """
    Ambiguity-aware %GC calculation
    """
    iupac = {
        'a': 0.0,
        'c': 1.0,
        'g': 1.0,
        't': 0.0,
        'u': 0.0,
        'r': 0.5,
        'y': 0.5,
        'm': 0.5,
        'k': 0.5,
        'w': 0.0,
        's': 1.0,
        'b': 2.0/3.0,
        'd': 1.0/3.0,
        'h': 1.0/3.0,
        'v': 2.0/3.0,
        'n': 0.5,
        
        'A': 0.0,
        'C': 1.0,
        'G': 1.0,
        'T': 0.0,
        'U': 0.0,
        'R': 0.5,
        'Y': 0.5,
        'M': 0.5,
        'K': 0.5,
        'W': 0.0,
        'S': 1.0,
        'B': 2.0/3.0,
        'D': 1.0/3.0,
        'H': 1.0/3.0,
        'V': 2.0/3.0,
        'N': 0.5,
    }
    gc_count = 0
    for s in seq:
        gc_count += iupac[s]
    
    return gc_count/len(seq)

def r_score(seq1, seq2, length):
    """
    Input should not include PAM sequence.
    
    nnnnnnnnnnnnnnnnnnnnPAM query
                    iiii--- no mismatch within 4 nt of PAM
                iiiiiiii--- no mismatch within 8 nt of PAM
            iiiiiiiiiiii--- no mismatch within 12 nt of PAM
        iiiiiiiiiiiiiiii--- no mismatch within 16 nt of PAM
    """
    if (ridentities(seq1, seq2) >= length):
        return 1.0
    else:
        return 0.0

#def score_histogram(sequence, algorithms):
#    """Code that scores a gRNA sequence
#    Returns scores"""
#    # two types of off-target scores
#    #  CFD off-target score
#    #  MIT off-target score
#    
#    # Histogram of off-targets:
#    #  For each number of mismatches, the number of off-targets is indicated.
#    #  Example:
#    #   1-3-20-50-60    This means 1 off-target with 0 mismatches, 3 off-targets with 1 mismatch, 20 off-targets with 2 mismatches, etc.
#    #   0-2-5-10-20     These are the off-targets that have no mismatches in the 12 bp adjacent to the PAM. These are the most likely off-targets.
#    #   
#    #   Off-targets are considered if they are flanked by one of the motifs NGG, NAG or NGA.


# So far for DNA only
def build_regex_pattern(iupac_sequence, max_substitutions=0, max_insertions=0, max_deletions=0, max_errors=0, capture=True):
    '''
    Build a regular expression pattern for the nucleotide search, taking IUPAC
    ambiguities into account.
    :param iupac_sequence: Sequence to translate into a (fuzzy) regex pattern
    :param max_substitutions: Max number of permitted substitutions. 'None'/'inf' if any amount permitted.
    :param max_insertions: Max number of permitted insertions. 'None'/'inf' if any amount permitted.
    :param max_deletions: Max number of permitted deletions. 'None'/'inf' if any amount permitted.
    :param max_errors: Max number of permitted errors. 'None'/'inf' if any amount permitted.
    :param capture: Whether or not to use a capturing group or non-capturing group
    :return: String with fuzzy regex pattern
    '''
    
    # Format the fuzzy matching restrictions
    fuzzy = []
    if (max_substitutions == None):
        fuzzy.append('s')
    elif (max_substitutions > 0):
        fuzzy.append('s<={}'.format(max_substitutions))
    
    if (max_insertions == None):
        fuzzy.append('i')
    elif (max_insertions > 0):
        fuzzy.append('i<={}'.format(max_insertions))
    
    if (max_deletions == None):
        fuzzy.append('d')
    elif (max_deletions > 0):
        fuzzy.append('d<={}'.format(max_deletions))
    
    if (max_errors == None):
        fuzzy.append('e')
    elif (max_errors > 0):
        fuzzy.append('e<={}'.format(max_errors))
    
    if (len(fuzzy) > 0):
        fuzzy = '{' + ','.join(fuzzy) + '}'
    else:
        fuzzy = ''
    
    # TODO: Allow ambiguities to match to themselves. For instance:
    #       iupac = {'k': '[kgt]', ...}
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
    #logger.info('Built regex string: {!r}'.format(pattern))
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

def random_disambiguated_sequence(sequence):
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
        
        '.': ['.'],
    }
    out = []
    for s in sequence:
        out.append(random.choice(iupac[s]))
    return ''.join(out)

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
    # Looks like the bugs with this have been fixed!
    m = regex.match(r'(?:'+seq1+r'){e}', seq2, flags=regex.BESTMATCH|regex.IGNORECASE)
    
    # Subs, ins, and dels are all given the same penalty
    
    # In general, (regex.match(seq1, seq2) != regex.match(seq2, seq1))
    # Also, terminal in/dels tend to be ignored (but not always)
    # Also, substitutions are favored over in/dels
    
    # substitutions, insertions, deletions
    return m.fuzzy_counts

def split_target_sequence2(seq, pams):
    """
    Searches for all PAMs in sequence, and splits it. Assumes PAM at 3' end.
    Splits at the shortest PAM...
    Returns gRNA, PAM
    """
    # Build a regex to only match strings with PAM sites specified in args.pams
    re_pattern = r'^(.*)(?:' + '|'.join(map(lambda x: build_regex_pattern(x), pams)) + r')$'
    #m = regex.match(r'^(.*)([ACGT]GG)$', nt)
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

def build_alignment(template_label, template_region_sequencess, template_region_labels, row_labels, row_sequences):
    '''
    a = build_alignment(
      'SEQ',
      ['AAAA', 'CCCCGGGGG', 'TTTTTTTT'],
      ['first', 'second', 'third'],
      ['contig1', 'contig2'],
      ['   ACCCC-GGGGTT', ' AAACCCCGGGGGTTT']
    )
    
    for line in a:
      print(a)
    
             ┌first
            ┌┴─┐┌second─┐┌third─┐
        SEQ AAAACCCCGGGGGTTTTTTTT
    contig1    ACCCC-GGGGTT
    contig2  AAACCCCGGGGGTTT
    
    '''
    pass

def make_labeled_primer_alignments(label_list, sequence_list, contig_name, primer_pair_label_list, primer_pair_list, shifts=None, left_pos=0):
    
    lengths = [len(x) for x in sequence_list]
    output = make_alignment_labels(label_list, lengths)
    output_names = ['' for x in output]
    
    output_names.append(contig_name)
    output.append(''.join(sequence_list))
    
    for i, (pp_label, pp) in enumerate(zip(primer_pair_label_list, primer_pair_list)):
        
        cn = contig_name.split(' ')[1]
        for f_gene, f_locus, f_genome, f_region, f_contig, f_strand, f_start, f_end in pp.forward_primer.locations:
            if (f_contig == cn):
                break
        
        output_names.append(pp_label)
        if (shifts == None):
            #output.append(' '*pp.forward_primer.position + pp.get_formatted())
            #for fseq in pp.get_formatted():
            #    output.append(' '*f_start + fseq[2])
            output.append(' '*(f_start-left_pos) + pp.get_formatted()[0][2])
        else: # Need to improve this so it uses shifts for real
            if pp:
                if (i in [0, 1]):
                    #output.append(' '*pp.forward_primer.position + pp.get_formatted())
                    output.append(' '*(f_start-left_pos) + pp.get_formatted()[0][2])
                elif (i == 2):
                    #output.append(' '*(lengths[0]+lengths[1]) + ' '*pp.forward_primer.position + pp.get_formatted()) # NOT elegant AT ALL
                    output.append(' '*(lengths[0]+lengths[1]) + ' '*(f_start-left_pos) + pp.get_formatted()[0][2])
            else:
                output.append('None')
    
    n = max(len(x) for x in output_names)
    lines = []
    for i in range(len(output)):
        lines.append(output_names[i].ljust(n, ' ') + ' ' + output[i])
    
    return lines

def make_alignment_labels(labels, lengths, edges_open=False):
    """
    Prints multi-line set of exclusive sequences and their labels
    
    Example:
      print_alignment_labels(['us region', 'label', 'us homology', 'insert',
          'ds homology', 'ds region'], [8, 5, 15, 3, 15, 10])
       ┌us region
       │       ┌label              ┌insert           ┌ds region
      ┌┴─────┐┌┴──┐┌─us homology─┐┌┴┐┌─ds homology─┐┌┴───────┐
    """
    lines = ['']
    for i, (label, slen) in enumerate(zip(labels, lengths)):
        
        # Modify the labels
        if (len(label)+2 > slen):
            if (slen == 0):
                if (len(label) > 0):
                    corner = '┌'
                else:
                    corner = ''
            elif (slen == 1):
                if (len(label) > 0):
                    corner = '┌'
                else:
                    corner = ''
            else:
                corner = ' ┌'
            #for j in range(1, len(lines)):
            for j in range(len(lines)-1, 0, -1):
                if (len(lines[j]) <= len(lines[0])):
                    lines[j] += ' '*(len(lines[0])-len(lines[j])) + corner+label
                    break
            else:
                if (len(lines) > 1):
                    new_line = ''
                    for k in range(0, len(lines[0])):
                        if (lines[-1][k] in ['┌', '│', '¦']):
                            if lines[0][k] in ['│', '┤', '┴']:
                                new_line += '│'
                            else:
                                new_line += '¦'
                        else:
                            new_line += ' '
                    if (lengths[i-1] == 0):
                        if (slen == 0):
                            corner = '┌'
                        else:
                            corner = '¦┌'
                    lines.append(new_line + corner+label)
                else:
                    lines.append(' '*len(lines[0]) + corner+label)
        
        # Modify lines[0]
        if (slen == 0):
            pass
        elif (slen == 1):
            if (len(label) > 0):
                lines[0] += '│'
            else:
                lines[0] += '╥'
        elif (slen == 2):
            if (len(label) > 0):
                lines[0] += '┌┤'
            else:
                lines[0] += '┌┐'
        elif (len(label)+2 > slen):
            lines[0] += '┌┴'+'─'*(slen-3)+'┐'
        else:
            line = '┌'+'─'*(slen-2)+'┐'
            pos = slen//2-len(label)//2
            lines[0] += line[:pos] + label + line[pos+len(label):]
    
    #for line in lines[1:]:
    #    print(line)
    #print(lines[0])
    return lines[1:] + [lines[0]]

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
    
    print("=== kmer distribution ===")
    seq = '''\
        GAGTCACGCCAATCACAAATTCCTTTGAAAAACTTGATTCGACCACATTCACAAGTTTGA
        TTGATTTGAAAAACTTGATTCGACACCATCCTGCTGTCCATCCGTGAGCCACACAGATTC
        AGAATTGAGTCGCTGACTAAGCGGTTAGACATACGTGATATTCACCGACTTTGAGAGTCC
        CACTAATCGGCTAGACATACGTAAATTACATAGCTCCCTCCAATACACACCCTACTTACT
        ATTGTCTTTTTTTAACTTTTTCGTAATCTCTACCCATAAAAATACACTTTCCCTCCAAAT
        CTCTAATTTACAACTCAACTGAACTTTAATTAACCTCTACTGCCTTAATTTAAGCTTATT
        TCTTGTCTATCAGCTGTTTCTGTTTCACCATTTTCACAACTTCTCCCCTAGGTGACATTT
        TTTTCTGCTGATTTTTTCTCAAATTCAGCCCAAAAAACTTAAACCAAAACTCAAAATTAC
        AACGCAAACTCTATTTAGAGTGCCCCTACTACCCCTACTGAGTCTTATTTTGAGTTTACC
        ACCGATTTCTGTGCTCCTCCTGTCTCCAGATTTCCGGTCTTCGTTCTTTTTTCGATCGAA
        AACTTTGTAAAACTAAACTAAAAAATTCACTCCATTTGACCAACAAASTGCTCAAAATCA
        GACCAGGCTCACTGCTTCTGCTTTGTCCCTAAAGATTACAAAAGCTACGCTGCAAAAGAA
        CTTAAAATTGCGTTCCATTATAATCTATACACACCCATCTCCTGCTATCACTTCACCTCA
        CGTCCTCCCTGCGCTTGTCCATCCGTGAGTTCAACTACCGCCTCCCTCTTCCCTTGTCCA
        CCCGTGATTCGCCAGTCCCTGGCTCTCCATCTTCCACAGATCCTTCACTTGCTTTCCATT
        GACTATCTTCTTCTCTTGCCCTAGCTTTTGATTTCCATATTCCTTCAACCATTGTACTAA
        CTCTCTCTTTACTCTGTGCTTAACTACTATCTCTCTGATCACCTGGCCTGGCGTTATTCT
        ATTTCCAGTTTTTTTTTTTTTCATTGATCCAACACAACTTCAACTCCCATTCGCTCGGCT
        CTTGACCCCCTTATCCATTCTCTCAGTACTTCCCGATCCCTTTTGTTCTTCATTACCCTT
        TTCTCTGTCTTGCCCTGCTACCCATCCGTGATTYTCCAGCRCTGTTCACTCCCACGTCCC
        CGCTGTTGATTGACATTTCCAATTTCACTGACTTTGTTCCCCTACTTTTGCTCACATTTT
        TCTGTTCTCAAACTCCTCTCTTGAATTCTCAGCTTGCTGTGTCTCCTTCTTGCCATTACA
        ACTGCTTTTCTTCACTTGCTTCCTTCTGCTTTGACAACACTGATCATTGACTTGATTTCA
        TTACTTTTCACAAACCCAGTTTCTAGCTCTATTGACTTCCTCTGCTATCCAGATTTCAAA
        CTTCTTATTGTAACAGTTATAACTGCGTTCTTCATCTCATCTAATTGATTGATTTGTTGT
        CGTTGAAGAAAAGTGATATTTTTTGACCAGCACATTTCTTGTCCAACTTTTTTTCGATGW
        CTTCTCCACACTTTTCTGCCACGTTTTCCCTATTTTTTTTGCCACGTCAGAAAAAAAAAA
        ATTTTTTTCACCACTTTTCTTCCCACCGCCAACAACACCAATGATGTTCTACCTGCCAGA
        GTGCCAGTTCTACATATGTTCCGATTTCCTAGCTCTTCAGATTCAGCAACTCCAACTACC
        AATTTTTGAATTCCCACAATCCAACTAATTCCCCGCCATCTTGCMAACTCAGTCCACAAT
        TTCTGTCCAACYACAAATTTTCAAACTGCAACAACTGTCACTGCCACATGCTATTCAACC
        GGCAAACAWACGAARCTGTAATGATTTCAACAACTGCCATTGATCACTCATTTATCAACC
        ACCAAACACAGCAGCGCAACAGCTTCCACAGTTCTTGTTGCCACGATTTCGGCAACTACG
    '''
    seq = seq.replace('\n', '').replace(' ', '')
    kmer_step_size = 1
    kmer_length = 4
    dist = get_seq_dist(seq, kmer_step_size, kmer_length)
    print('N', 'kmer', 'rc()', 'palindrome', 'frequency', 'count')
    for i, (k, v) in enumerate(dist.items()):
        m = regex.match(r'(.)(.)(.)(.) \1\2\3\4',k)
        print(i, k, 'palindrome' if m else '          ', v, v*(len(seq)-kmer_length+1))
    print('sum =', sum(dist.values()))
    
    print('=== CDIST ===')
    cdist = complement_kmer_distribution(dist)
    print('N', 'kmer', 'rc()', 'palindrome', 'frequency')
    for i, (k, v) in enumerate(cdist.items()):
        m = regex.match(r'(.)(.)(.)(.) \1\2\3\4',k)
        print(i, k, 'palindrome' if m else '          ', v)
    print('sum =', sum(cdist.values()))
    
    seqs = [
        'AGCTACGAGCATCAGGACTACGG',
        'GTCAGACTAGCGACATATATAGG',
        'CCAGTAGAGAATAGAGCCTAGAC',
        'CCATCGAACTACCCTCTCTCAAC'
    ]
    for s in seqs:
        print(s, sequence_likelihood(s, dist), sequence_likelihood(s, cdist))
    
    

if (__name__ == '__main__'):
    test()
