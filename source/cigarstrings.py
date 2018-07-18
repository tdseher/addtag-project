#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/aligners/bowtie2.py

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
# Notes
#  * H can only be present as the first and/or last operation.
#  * S may only have H operations between them and the ends of the CIGAR string.
#  * For mRNA-to-genome alignment, an N operation represents an intron. For other types of alignments, the interpretation of N is not defined.
#  * Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.

def get_errors(cigar):
    """
    Returns the number of mismatches, insertions, deletions, and errors
    Only works for specific CIGARS containing =/X instead of M
    as defined by the CIGAR string.
    """
    mismatches = 0
    insertions = 0
    deletions = 0
    
    m = regex.findall(r'(\d*)([IDX])', cigar)
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
    errors = mismatches + insertions + deletions
    
    return mismatches, insertions, deletions, errors

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
    """Decodes the bit flag column in SAM output into a string or list of characters"""
    # Flag    Hexadecimal   Binary Decimal  FlagChar  Description                                       New Description
    # 0x0000  0x0                0       0            read maps?
    # 0x0001  0x1                1       1  p         the read is paired in sequencing                  template having multiple segments in sequencing
    # 0x0002  0x2               10       2  P         the read is mapped in a proper pair               each segment properly aligned according to the aligner
    # 0x0004  0x4              100       4  u         the query sequence itself is unmapped             segment unmapped
    # 0x0008  0x8             1000       8  U         the mate is unmapped                              next segment in the template unmapped
    # 0x0010  0x10           10000      16  r         strand of the query (1 for reverse)               SEQ being reverse complemented
    # 0x0020  0x20          100000      32  R         strand of the mate                                SEQ of the next segment in the template being reversed
    # 0x0040  0x40         1000000      64  1         the read is the first read in a pair              the first segment in the template
    # 0x0080  0x80        10000000     128  2         the read is the second read in a pair             the last segment in the template
    # 0x0100  0x100      100000000     256  s         the alignment is not primary                      secondary alignment
    # 0x0200  0x200     1000000000     512  f         the read fails platform/vendor quality checks     not passing quality controls
    # 0x0400  0x400    10000000000    1024  d         the read is either a PCR or an optical duplicate  PCR or optical duplicate
    # 0x0800  0x800   100000000000    2048  a                                                           supplementary alignment
    
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

def test():
    q = "AAA--AAATCCC"
    s = "AGACGAAA-CCC"
    alignment2cigar(q, s) # '3M2D3M1I3M'
    alignment2cigar(q, s, abbreviated=True) # '3M2D3MI3M'
    alignment2cigar(q, s, specific=True) # '1=1X1=2D3=1I3='
    alignment2cigar(q, s, specific=True, abbreviated=True) # '=X=2D3=I3='

if (__name__ == '__main__'):
    test()
