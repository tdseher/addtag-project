#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/utils.py

# Import standard packages
import sys
import os
import gzip
import datetime
import logging

# Import non-standard packages
import regex

logger = logging.getLogger(__name__)

def flatten(iterable, remove_none=False):
    """Make a flat list out of a list of lists"""
    # Will remove None
    if remove_none:
        return [item for sublist in iterable for item in sublist if item != None]
    else:
        return [item for sublist in iterable for item in sublist]

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
    
    #m = regex.findall(r'(\d+)[MISP=X]', cigar)
    m = regex.findall(r'(\d+)[MDSP=X]', cigar)
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

def load_fasta_file(filename):
    """
    Load contig sequences from file into dict()
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
    
    logger.info('FASTA file parsed: {!r}'.format(filename))
    return contigs

def load_gff_file(filename, features, tag):
    """
    Load General Feature Format (GFF) file into dict()
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
                        #           gene     contig       start(bp)          end(bp)    strand
                        annotations[m[0]] = (sline[0], int(sline[3])-1, int(sline[4]), sline[6])
    
    logger.info('GFF file parsed: {!r}'.format(filename))
    return annotations

def load_git_date():
    '''
    Returns the date of the most-recent git commit for the current head
    and master.
    '''
    # Repository versions stored in this file:
    #   .../addtag-project/.git/logs/regs/heads/master
    root = os.path.join(os.path.dirname(__file__), '..')
    filename = os.path.join('.git', 'logs', 'refs', 'heads', 'master')
    filepath = os.path.join(root, filename)
    seconds = []
    try:
        with open(filepath, 'r') as flo:
            for line in flo:
                m = regex.search(' (\d{10,}) ([-+]?\d{4})\t', line) # only works for ints
                if m:
                    secs, off = m.groups()
                    # seconds.append(int(secs) + int(off[:-2])*60*60) # Convert to GMT time
                    seconds.append(int(secs)) # Don't convert time, as it will display in the correct local time
        date = str(datetime.datetime.fromtimestamp(seconds[-1]))
    except FileNotFoundError:
        date = 'missing'
    return date

def load_git_commits():
    '''Returns the git commit number for the current head and master'''
    # Repository versions stored in this file:
    #   .../addtag-project/.git/logs/regs/heads/master
    root = os.path.join(os.path.dirname(__file__), '..')
    filename = os.path.join('.git', 'logs', 'refs', 'heads', 'master')
    filepath = os.path.join(root, filename)
    #filepath = os.path.realpath(os.path.join(root, filename))
    versions = set()
    try:
        with open(filepath, 'r') as flo:
            for line in flo:
                versions.add(line.rstrip().split(' ', 2)[1])
        version = len(versions)
    except FileNotFoundError:
        version = 'missing'
    return version

def load_git_version():
    '''Returns the git hash for the current head and master'''
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

def merge_sets(sets):
    """
    Merge sets into non-overlapping groups.
    """
    groups = []
    pending_sets = sets[:]
    while (len(pending_sets) > 0):
        current, *rest = pending_sets
        current = set(current)
        
        # keep looping through the rest as long as the current set
        # is being updated
        previous_length = -1
        while (previous_length < len(current)):
            previous_length = len(current)
    
            rest2 = []
            # Loop through the remaining sets
            for r in rest:
                # If the remaining set overlaps with the current one, then
                # update the current one
                if (len(current.intersection(set(r))) > 0):
                    current.update(set(r))
                # Otherwise, add the non-overlapping sets to rest
                else:
                    rest2.append(r)     
            rest = rest2
    
        groups.append(current)
        pending_sets = rest
    
    return groups

def load_homologs(filename, sep="\t"):
    """
    Parse homolog file.
    Arguments
     - filename
     - column separator ("\t" by default)
     - row headers (True by default)
    Returns
     - dictionary of homolog groups
    """
    # homologs = {}
    # with open(filename, 'r') as flo:
    #     for line in flo:
    #         sline = line.rstrip().split(sep)
    #         if headers:
    #             row = sline[1:]
    #         else:
    #             row = sline
    #         for h1 in row:
    #             for h2 in row:
    #                 if (h1 != h2):
    #                     homologs.setdefault(h1, set()).add(h2)
    
    feature2gene = {}
    
    sets = []
    with open(filename, 'r') as flo:
        for line in flo:
            sline = line.rstrip().split(sep)
            sets.append(set(sline[1:]))
            for f in sline[1:]:
                feature2gene[f] = sline[0]
    
    groups = merge_sets(sets)
    homologs = {}
    for group in groups:
        for g in group:
            homologs[g] = group
    
    #print("homologs:", file=sys.stderr)
    #for h in homologs:
    #    print(' ', h, homologs[h], file=sys.stderr)
    logger.info('Homologs file parsed: {!r}'.format(filename))
    return homologs, feature2gene

def test():
    #print(load_git_version())
    alignments = load_sam_file(sys.argv[1])
    for a in alignments:
        for b in sorted(alignments[a]):
            print(a, b, alignments[a][b])

if (__name__ == "__main__"):
    test()