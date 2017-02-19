#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/utils.py

# Import standard packages
import sys
import gzip
import regex

class sliding_window(object):
    """Iterator for a sliding window across a nucleotide sequence"""
    def __init__(self, sequence, window, start=0, stop=None, step=1):
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
        if ((self.position >= self.stop) or (self.position + self.window > len(self.sequence))):
            raise StopIteration
        seq = self.sequence[self.position:self.position+self.window]
        self.position += self.step
        return seq

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

def load_fasta_file(filename):
    """Load contig sequences from file into dict()
    Primary sequence headers must be unique
    Can open a *.gz compressed file
    """
    
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

def disambiguate_iupac(iupac_sequence):
    """converts a string containing IUPAC nucleotide sequence to a list
    of non-iupac sequences.
    """
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
    
    # Throw an error if there is a non-cannonical nucleotide
    assert all(map(lambda seq: all(map(lambda x: x in 'AaCcGgTtUu', seq)), sequences)) == True
    
    return sequences

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

