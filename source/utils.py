#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/utils.py

# Import standard packages
import sys
import os
import gzip
import regex

class sliding_window(object):
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

def flatten(iterable, remove_none=False):
    """Make a flat list out of a list of lists"""
    # Will remove None
    if remove_none:
        return [item for sublist in iterable for item in sublist if item != None]
    else:
        return [item for sublist in iterable for item in sublist]

def lcs(string1, string2):
    """Find the longest common substring between two strings"""
    import difflib
    
    matcher = difflib.SequenceMatcher(None, string1, string2, True)
    match = matcher.find_longest_match(0, len(string1), 0, len(string2))
    # Match(a=0, b=15, size=9)
    return match
    # print(string1[match.a: match.a + match.size])
    # print(string2[match.b: match.b + match.size])

def rc(seq, kind="dna"):
    """Returns the reverse-complement of a string containing DNA or RNA characters"""
    if (kind == "dna"):
        complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB') # exclude ws, WS
    elif (kind == "rna"):
        complements = str.maketrans('acgturymkbdhvACGTURYMKBDHV', 'ugcaayrkmvhdbTGCAAYRKMVHDB') # exclude ws, WS
    else:
        raise ValueError("'" + str(kind) + "' is an invalid argument for rc()")
    return seq.translate(complements)[::-1]

def linear_score(seq1, seq2):
    """Scores low if substitutions near 3' end of the sequence
    To add: insertions and deletions """
    x = list(range(20))
    t = sum(x)
    y = list(map(lambda i: i/t, x))
    score = 1
    for i in range(len(seq1)):
        if (seq1[i] != seq2[i]):
            score -= y[i]
    return score*100

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

def cigar_length(cigar):
    """Returns the length of the sequence specified by the CIGAR string"""
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
    
    m = regex.findall(r'(\d+)[MISP=X]', cigar)
    return sum(map(int, m))

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

def load_sam_file(filename, sep=':'):
    """Read in SAM file.
      sep is the separator for the header. Positions are converted to 0-index
    """
    
    # Sequence Alignment/Map (SAM) format is TAB-delimited. Apart from the
    # header lines, which are started with the '@' symbol, each alignment line
    # consists of:
    #   Col  Field  Description
    #   1    QNAME  Query template/pair NAME
    #   2    FLAG   bitwise FLAG
    #   3    RNAME  Reference sequence NAME
    #   4    POS    1-based leftmost POSition/coordinate of clipped sequence
    #   5    MAPQ   MAPping Quality (Phred-scaled)
    #   6    CIAGR  extended CIGAR string
    #   7    MRNM   Mate Reference sequence NaMe ('=' if same as RNAME)
    #   8    MPOS   1-based Mate POSistion
    #   9    TLEN   inferred Template LENgth (insert size)
    #   10   SEQ    query SEQuence on the same strand as the reference
    #   11   QUAL   query QUALity (ASCII-33 gives the Phred base quality)
    #   12+  OPT    variable OPTional fields in the format TAG:VTYPE:VALUE
    
    # Code to decompress a *.bam file should go here
    
    alignments = {}
    with open(filename, 'r') as flo:
        for line in flo:
            if not line.startswith('@'):
                sline = line.rstrip().split("\t")
                if (len(sline) > 5):
                    feature, source_contig, source_start, source_end = sline[0].split(sep)
                    
                    source = (source_contig, int(source_start), int(source_end))
                    # convert mapping to 0-based index
                    dest = (sline[2], int(sline[3])-1, int(sline[3])-1+cigar_length(sline[5]), sam_orientation(int(sline[1])))
                    try:
                        alignments[feature][source].append(dest)
                    except KeyError:
                        try:
                            alignments[feature][source] = [dest]
                        except KeyError:
                            alignments[feature] = {source: [dest]}
    
    return alignments

def load_fasta_file(filename):
    """Load contig sequences from file into dict()
    Primary sequence headers must be unique
    Can open a *.gz compressed file
    """
    # Future features: parse the non-primary sequence header to find tags
    # specifying whether the contig is linear or circular:
    #  linear=true, linear=yes, linear=1
    #  circular=true, circular=yes, circular=1
    #  etc
    # instead of storing the values as strings, use tuples:
    #  contigs['header'] = ('ACGTGACGA', 'linear')
    
    if filename.endswith('.gz'):
        flo = gzip.open(filename, 'rt')
    else:
        flo = open(filename, 'r')
    
    contigs = {}
    #with open(filename, 'r') as flo:
    with flo:
        name = None
        seq = None
        for line in flo:
            line = line.rstrip()
            if line.startswith('>'):
                if (name != None):
                    contigs[name] = seq
                name = regex.split(r'\s+', line[1:], 1)[0]
                seq = ''
            else:
                seq += line
        if (name != None):
            contigs[name] = seq
    
    print('FASTA file parsed: {!r}'.format(filename), file=sys.stderr)
    return contigs

def load_gff_file(filename, features, tag):
    """Load General Feature Format (GFF) file into dict()
    One line per feature, each containing 9 columns of data, plus optional
    track definition lines.
    
    Converts positions to 0-based index.
    """
    # Fields must be tab-separated. Also, all but the final field in each
    # feature line must contain a value; "empty" columns should be denoted
    # with a '.'
    #  1) seqid - name of the chromosome or scaffold; chromosome names can
    #     be given with or without the 'chr' prefix. Important note: the
    #     seqname must be a standard chromosome name or an identifier such as
    #     a scaffold ID, without any additional content such as species or
    #     assembly.
    #  2) source - name of the program that generated this feature, or the
    #     data source (database or project name)
    #  3) feature - feature type name, e.g. Gene, Variation, Similarity
    #  4) start - Start position of the feature, with sequence numbering
    #     starting at 1.
    #  5) end - End position of the feature, with sequence numbering
    #     starting at 1.
    #  6) score - A floating point value. As in earlier versions of the format,
    #     the semantics of the score are ill-defined. It is strongly
    #     recommended that E-values be used for sequence similarity features,
    #     and that P-values be used for ab initio gene prediction features.
    #     If there is no score, put a '.' (a period) in this field.
    #  7) strand - defined as '+' (forward), '-' (reverse), '.' (unstranded),
    #     '?' (relevant, but unknown).
    #  8) frame - for CDS features, '0', '1' or '2'. '0' indicates that the first base of
    #     the feature is the first base of a codon, '1' that the second base
    #     is the first base of a codon, and so on. Other features can use '.'.
    #  9) attribute - A semicolon-separated list of tag-value pairs, providing
    #     additional information about each feature. A list of feature
    #     attributes in the format tag=value. Multiple tag=value pairs are
    #     separated by semicolons. URL escaping rules are used for tags or
    #     values containing the following characters: ",=;". Spaces are allowed
    #     in this field, but tabs must be replaced with the %09 URL escape.
    #     This field is not required.
    #     Column 9 tags have predefined meanings:
    #       ID - Indicates the unique identifier of the feature. IDs must be
    #            unique within the scope of the GFF file.
    #       Name - Display name for the feature. This is the name to be displayed to the user. Unlike IDs, there is no requirement that the Name be unique within the file.
    #       Alias - A secondary name for the feature. It is suggested that this tag be used whenever a secondary identifier for the feature is needed, such as locus names and accession numbers. Unlike ID, there is no requirement that Alias be unique within the file.
    #       Parent - Indicates the parent of the feature. A parent ID can be used to group exons into transcripts, transcripts into genes, and so forth. A feature may have multiple parents. Parent can *only* be used to indicate a partof relationship.
    #       Target - Indicates the target of a nucleotide-to-nucleotide or protein-to-nucleotide alignment. The format of the value is "target_id start end [strand]", where strand is optional and may be "+" or "-". If the target_id contains spaces, they must be escaped as hex escape %20.
    #       Gap - The alignment of the feature to the target if the two are not collinear (e.g. contain gaps). The alignment format is taken from the CIGAR format described in the Exonerate documentation. http://cvsweb.sanger.ac.uk/cgi-bin/cvsweb.cgi/exonerate?cvsroot=Ensembl). See the GFF3 specification for more information.
    #       Derives_from - Used to disambiguate the relationship between one feature and another when the relationship is a temporal one rather than a purely structural "part of" one. This is needed for polycistronic genes. See the GFF3 specification for more information.
    #       Note - A free text note.
    #       Dbxref - A database cross reference. See the GFF3 specification for more information.
    #       Ontology_term - A cross reference to an ontology term. See the GFF3 specification for more information.
    #     Multiple attributes of the same type are indicated by separating the
    #     values with the comma "," character, as in: 'Parent=AF2312,AB2812,abc-3'
    #     Note that attribute names are case sensitive. "Parent" is not the
    #     same as "parent". All attributes that begin with an uppercase letter
    #     are reserved for later use. Attributes that begin with a lowercase
    #     letter can be used freely by applications. You can stash any
    #     semi-structured data into the database by using one or more
    #     unreserved (lowercase) tags.
    
    feature = features[0]
    
    #contig_lengths = {}
    annotations = {}
    
    with open(filename, 'r') as flo:
        for line in flo:
            line = line.rstrip()
            #if line.startswith("##sequence-region"):
            #    # sample line:
            #    # ##sequence-region contig_1 1 209073
            #    sline = line.split(' ')
            #    contig_lengths[sline[1]] = int(sline[3])
            #else:
            if not line.startswith('#'):
                sline = line.split('\t')
                if ((len(sline) > 6) and (sline[2] == feature)):
                    m = regex.findall(tag + r'=([^;]*)', sline[8])
                    if m:
                        #           gene     contig       start(bp)          end(bp)      strand
                        annotations[m[0]] = (sline[0], int(sline[3])-1, int(sline[4])-1, sline[6])
    
    print('GFF file parsed: {!r}'.format(filename), file=sys.stderr)
    return annotations

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

def load_git_version():
    '''Returns the git version for the current head and master'''
    # Current repository version stored in this file:
    #   .../addtag-project/.git/refs/heads/master
    root = os.path.join(os.path.dirname(__file__), '..')
    filename = os.path.join('.git', 'refs', 'heads', 'master')
    filepath = os.path.join(root, filename)
    version = None
    try:
        with open(filepath, 'r') as flo:
            version = flo.readline().rstrip()
    except FileNotFoundError:
        version = 'missing'
    return version

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

def test():
    #print(load_git_version())
    alignments = load_sam_file(sys.argv[1])
    for a in alignments:
        for b in sorted(alignments[a]):
            print(a, b, alignments[a][b])

if (__name__ == "__main__"):
    test()