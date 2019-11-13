#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/cigarstrings.py

# List general Python imports

# import non-standard package
import regex

# import AddTag-specific packages

# The CIGAR operations are given in the following table (set '*' if unavailable):
#   Op BAM Description
#   M    0 alignment match (can be a sequence match or mismatch)
#   I    1 insertion to the reference
#   D    2 deletion from the reference
#   N    3 skipped region from the reference
#   S    4 soft clipping (clipped sequences present in SEQ)
#   H    5 hard clipping (clipped sequences NOT present in SEQ)
#   P    6 padding (silent deletion from padded reference)
#   =    7 sequence match
#   X    8 sequence mismatch
#
# Notes:
#  * H can only be present as the first and/or last operation.
#  * S may only have H operations between them and the ends of the CIGAR string.
#  * For mRNA-to-genome alignment, an N operation represents an intron.
#    For other types of alignments, the interpretation of N is not defined.
#  * Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.

def specificate_cigar(cigar, query, subject):
    '''
    'specific' is whether to use '=' and 'X' for match and mismatch, or just 'M'
    
    :param cigar: CIGAR string containing 'M'
    :param query: Aligned query sequence
    :param subject: Aligned subject sequence
    :return: CIGAR string containing '=' or 'X' instead of 'M'
    '''
    out = []
    matches = regex.findall(r'(\d*)([MIDNSHP=X])', cigar)
    qi = 0
    si = 0
    for match_n, match_c in matches:
        if (match_n == ''):
            n = 1
        else:
            n = int(match_n)
        if (match_c == 'M'):
            for i in range(n):
                if (query[qi] == subject[si]):
                    out.append('=')
                else:
                    out.append('X')
                qi += 1
                si += 1
        elif (match_c == 'I'):
            for i in range(n):
                out.append('I')
                qi += 1
        elif (match_c == 'D'):
            for i in range(n):
                out.append('D')
                si += 1
        elif (match_c == 'N'):
            # TODO: Confirm that skipped regions are treated equivalently to deletions
            for i in range(n):
                out.append('N')
                si += 1
        elif (match_c == 'S'):
            # TODO: Confirm that soft-masking affects both query and subject coordinates
            for i in range(n):
                out.append('S')
                qi += 1
                si += 1
        elif (match_c == 'H'):
            # TODO: Confirm that hard-masking affects only subject coordinates
            for i in range(n):
                out.append('H')
                si += 1
        elif (match_c == 'P'):
            # TODO: Confirm that 'silent deletions' are treated equivalently to deletions
            for i in range(n):
                out.append('P')
                si += 1
        elif (match_c in ['=', 'X']):
            for i in range(n):
                out.append(match_c)
                qi += 1
                si += 1
    return ''.join(out)

def collapse_cigar(cigar, abbreviate=False):
    '''
    Preface repeated characters with digit quantifiers
    :param cigar: 
    :param abbreviate Remove preceeding 1 for non-repeated codes
    :return: 
    '''
    # Consolodate consecutive repeating characters
    out = []
    
    prev_n = None
    prev_c = None
    
    matches = regex.findall(r'(\d*)([MIDNSHP=X])', cigar)
    
    for match_n, match_c in matches:
        if (match_n == ''):
            n = 1
        else:
            n = int(match_n)
        if (match_c == prev_c):
            prev_n += n
        elif (prev_c == None):
            prev_n = n
            prev_c = match_c
        else:
            if (abbreviate and (prev_n == 1)):
                pass
            else:
                out.append(str(prev_n))
            out.append(prev_c)
            prev_n = n
            prev_c = match_c
    
    if (prev_c != None):
        if (abbreviate and (prev_n == 1)):
            pass
        else:
            out.append(str(prev_n))
        out.append(prev_c)
    
    return ''.join(out)

def expand_cigar(cigar):
    '''
    Remove all digit quantifiers by fully expanding CIGAR string
    :param cigar: 
    :return: 
    '''
    out = []
    matches = regex.findall(r'(\d*)([MIDNSHP=X])', cigar)
    for match_n, match_c in matches:
        if (match_n == ''):
            n = 1
        else:
            n = int(match_n)
        for i in range(n):
            out.append(match_c)
    return ''.join(out)

def cigar_errors(cigar):
    """
    Returns the number of matches, mismatches, insertions, deletions, and errors
    Only works for specific CIGARS containing '=' or 'X' instead of 'M'
    as defined by the CIGAR string.
    """
    matches = 0
    mismatches = 0
    insertions = 0
    deletions = 0
    
    m = regex.findall(r'(\d*)([=IDX])', cigar)
    for quant, char in m:
        if quant:
            q = int(quant)
        else:
            q = 1
        if (char == 'I'):
            insertions += q
        elif (char == 'D'):
            deletions += q
        elif (char == 'X'):
            mismatches += q
        elif (char == '='):
            matches += q
    errors = mismatches + insertions + deletions
    
    return matches, mismatches, insertions, deletions, errors

def matrix_averages(scores):
    """
    Returns:
     * Average score of entire matrix
     * Average match score
     * Average mismatch score
     * Average insertion score
     * Average deletion score
    """
    total_s = 0
    total_n = 0
    
    match_s = 0
    match_n = 0
    
    mismatch_s = 0
    mismatch_n = 0
    
    ins_s = 0
    ins_n = 0
    
    del_s = 0
    del_n = 0
    
    for k, v in scores.items():
        # These are invalid so we can ignore them
        if ((len(k) == 3) and (k[0] == k[1])):
            continue # Skip this entry
        
        total_s += v
        total_n += 1
        
        # k[0] is the query
        # k[1] is the subject
        if (len(k) == 2):
            if (k[0] == k[1]):
                match_s += v
                match_n += 1
            else:
                mismatch_s += v
                mismatch_n += 1
        else: # (len(k) == 3)
            if (k[0] == '-'): # The query is missing these, so these are deletions
                del_s += v
                del_n += 1
            elif (k[1] == '-'): # The subject is missing these, so these are insertions
                ins_s += v
                ins_n += 1
    
    # entire_table_average, match_average, mismatch_average, insertion_average, deletion_average
    return match_s/match_n, mismatch_s/mismatch_n, ins_s/ins_n, del_s/del_n, total_s/total_n

def cigar2score(cigar, score_matrix):
    """
    Returns the score of the alignment, given the cigar string.
    """
    # matches, mismatches, insertions, deletions, errors = cigar_errors(cigar)
    # avg_match, avg_mismatch, avg_ins, avg_del, avg_total = matrix_averages(score_matrix)
    counts = cigar_errors(cigar)[:4] # remove errors total
    scores = matrix_averages(score_matrix)[:4] # remove total
    
    return sum(a*b for a,b in zip(counts, scores))

def cigar2query_position(cigar):
    """
    Finds the 0-indexed start/end position of the query.
    Uses hard masking positions
    Returns (start, end)
    """
    matches = regex.findall(r'(\d*)([MIDNSHP=X])', cigar)
    start = 0
    end = 0
    started = False
    #ended = False
    for q, char in matches:
        if q:
            quant = int(q)
        else:
            quant = 1
        if char in ['H']:
            if not started:
                start += quant
                end += quant
        elif char in ['M', 'I', 'S', 'P', '=', 'X']:
            started = True
            end += quant
    return start, end

def cigar2query_aligned_length(cigar):
    """
    Returns the length of the query sequence used to generate the
    CIGAR string. (Works for abbreviated CIGAR)
    """
    m = regex.findall(r'(\d*)[MISP=X]', cigar)
    return sum(int(x) if x else 1 for x in m)

def cigar2subject_aligned_length(cigar):
    """
    Returns the length of the subject sequence used to generate the
    CIGAR string. (Works for abbreviated CIGAR)
    """
    m = regex.findall(r'(\d*)[MDSP=X]', cigar)
    return sum(int(x) if x else 1 for x in m)

def reverse_cigar(cigar):
    """
    Reverses the CIGAR string
    """
    matches = [m.group(0) for m in regex.finditer(r'(\d*\D)', cigar)]
    return ''.join(matches[::-1])

# TODO: Fix this function (it only seems to work some of the time),
#       because 'regex' does its calculations in a cryptic way.
def match2cigar(m, specific=False, abbreviated=False):
    '''
    Input is a 'regex.Match' object that uses fuzzy matching.
    Outputs the CIGAR string describing the alignment.
    'specific' is whether to use =/X for match/mismatch, or just M
    'abbreviated' is whether to print the digit for non-repeated codes
    
    Returns incorrect output on certain regex (see regex-align-test.py).
    Otherwise, it is pretty close to correct.
    '''
    subs = m.fuzzy_changes[0]
    ins = m.fuzzy_changes[1]
    dels = m.fuzzy_changes[2]
    
    cigar = ''
    
    i = m.start()
    end = m.end()
    while (i < end):
        processed = False
        
        if i in subs:
            cigar += 'X'
            processed = True
        if i in ins:
            cigar += 'I'
            processed = True
        if i in dels:
            cigar += 'D'
            end += 1
            processed = True
        
        if not processed:
            cigar += '='
        i += 1
    
    #### Idea I haven't tested ####
    # repeats = 0
    # while (i < end):
    #     processed = False
    #     
    #     if i in subs:
    #         cigar += 'X'
    #         processed = True
    #     if i in ins:
    #         cigar += 'I'
    #         processed = True
    #     if i in dels:
    #         cigar += 'D'
    #         processed = True
    #     if not processed:
    #         cigar += '='
    #     
    #     if (repeats > 0):
    #         repeats -= 1
    #     else
    #         i += 1
    #### Idea I haven't tested ####
    
    if (specific == False):
        new_cigar = ''
        for c in cigar:
            if c in ['X', '=']:
                new_cigar += 'M'
            else:
                new_cigar += c
        cigar = new_cigar
    
    # Consolodate consecutive repeating characters
    new_cigar = ''
    matches = regex.finditer(r'(.)\1*', cigar)
    for m in matches:
        if (abbreviated and (len(m.group()) == 1)):
            new_cigar += m.group()[0]
        else:
            new_cigar += str(len(m.group())) + m.group()[0]
    cigar = new_cigar
    return cigar

def alignment2cigar(query, subject, specific=False, abbreviated=False):
    """
    Input is two sequences of equal length.
    Outputs the CIGAR string describing the alignment.
    'specific' is whether to use =/X for match/mismatch, or just M
    'abbreviated' is whether to print the digit for non-repeated codes
    """
    
    # Gracefully fail if the sequences are different lengths
    if (len(query) != len(subject)):
        return None
    
    cigar = ''
    current = ''
    previous = ''
    quantifier = ''
    
    q_iter = iter(query)
    s_iter = iter(subject)
    
    for qs in q_iter:
        ss = next(s_iter)
        
        if ((qs == '-') and (ss == '-')): # These should not happen, but if they dy
            continue # Skip this character
        elif (qs == '-'):
            current = 'D' # Deletion in query
        elif (ss == '-'):
            current = 'I' # Insertion in query
        else: # Aligned
            if (specific):
                if (qs == ss):
                    current = '=' # Match
                else:
                    current = 'X' # Mismatch
            else:
                current = 'M' # Aligned
        
        if (current == previous):
            quantifier += 1
        else:
            if (abbreviated and (quantifier == 1)):
                cigar += previous
            else:
                cigar += str(quantifier) + previous
            previous = current
            quantifier = 1
    
    if (abbreviated and (quantifier == 1)):
        cigar += previous
    else:
        cigar += str(quantifier) + previous
    
    return cigar

def sam_orientation(field):
    """Uses the bit flag column in SAM output to determine the orientation
    Returns
      '+' for forward orientation
      '-' for reverse-complement orientation
    """
    # Typecast if necessary
    if isinstance(field, str):
        field = int(field)
    
    if (field & 16):
        return '-'
    else:
        return '+'

def decode_sam_flags(field, kind='str'):
    """
    Decodes the bit flag column in SAM output into a string or list
    of characters.
    """
    # Flag    Hexadecimal     Binary Decimal  FlagChar  Description                                       New Description
    # 0x0000  0x0                  0       0            read maps?
    # 0x0001  0x1                  1       1  p         the read is paired in sequencing                  template having multiple segments in sequencing
    # 0x0002  0x2                 10       2  P         the read is mapped in a proper pair               each segment properly aligned according to the aligner
    # 0x0004  0x4                100       4  u         the query sequence itself is unmapped             segment unmapped
    # 0x0008  0x8               1000       8  U         the mate is unmapped                              next segment in the template unmapped
    # 0x0010  0x10            1_0000      16  r         strand of this read (query): '+'=0, '-'=1         SEQ being reverse complemented
    # 0x0020  0x20           10_0000      32  R         strand of the mate                                SEQ of the next segment in the template being reversed
    # 0x0040  0x40          100_0000      64  1         the read is the first read in a pair              the first segment in the template
    # 0x0080  0x80         1000_0000     128  2         the read is the second read in a pair             the last segment in the template
    # 0x0100  0x100      1_0000_0000     256  s         the alignment is not primary                      secondary alignment
    # 0x0200  0x200     10_0000_0000     512  f         the read fails platform/vendor quality checks     not passing quality controls
    # 0x0400  0x400    100_0000_0000    1024  d         the read is either a PCR or an optical duplicate  PCR or optical duplicate
    # 0x0800  0x800   1000_0000_0000    2048  a                                                           supplementary alignment
    
    flags = {1:'p', 2:'P', 4:'u', 8:'U', 16:'r', 32:'R', 64:'1', 128:'2',256:'s', 512:'f', 1024:'d', 2048:'a'}
    # Typecast if necessary
    if isinstance(field, str):
        field = int(field)
    
    str_flags = ''
    list_flags = []
    for k in flags:
        if (field & k):
            str_flags += flags[k]
            list_flags.append(flags[k])
    
    if (kind == 'str'):
        return str_flags
    elif (kind == 'list'):
        return list_flags
