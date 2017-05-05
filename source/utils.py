#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/utils.py

# Import standard packages
import sys
import os
import gzip
import datetime

# Import non-standard packages
import regex

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
                    #       contig,   contig_start,    contig_end,                             orientation(+/-)
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
                        #           gene     contig       start(bp)          end(bp)    strand
                        annotations[m[0]] = (sline[0], int(sline[3])-1, int(sline[4]), sline[6])
    
    print('GFF file parsed: {!r}'.format(filename), file=sys.stderr)
    return annotations

def load_git_date():
    '''Returns the date of the most-recent git commit for the current head and master'''
    # Repository versions stored in this file:
    #   .../addtag-project/.git/logs/regs/heads/master
    root = os.path.join(os.path.dirname(__file__), '..')
    filename = os.path.join('.git', 'logs', 'refs', 'heads', 'master')
    filepath = os.path.join(root, filename)
    seconds = []
    try:
        with open(filepath, 'r') as flo:
            for line in flo:
                m = regex.search(' (\d{10,}) ([-+]?\d{4})\tcommit', line) # only works for ints
                if m:
                    secs, off = m.groups()
                    # Convert to GMT time
                    seconds.append(int(secs) + int(off[:-2])*60*60)
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

def load_homologs(filename, sep="\t"):
    """
    Parse homolog file.
    Arguments
     - filename
     - column separator (optional)
    Returns
     - dictionary of homolog groups
    """
    homologs = {}
    with open(filename, 'r') as flo:
        for line in flo:
            sline = line.rstrip().split(sep)
            for h1 in sline[1:]:
                for h2 in sline[1:]:
                    if (h1 != h2):
                        homologs.setdefault(h1, set()).add(h2)
    #print("homologs:", file=sys.stderr)
    #for h in homologs:
    #    print(' ', h, homologs[h], file=sys.stderr)
    print('Homologs file parsed: {!r}'.format(filename), file=sys.stderr)
    return homologs

def generate_excision_spacers(filename, sequences, sep=':'):
    """
    Creates a FASTA file
    >feature:contig:orientation:start..end motif=STR azimuth=N off-target=N alignments=N
    """
    with open(filename, 'w') as flo:
        for s in sequences:
            #feature, contig, orientation, start, end, seq, motif, azimuth, off_target, alignments = line
            print(">" + s.feature + sep + s.contig + sep + s.contig_orientation + \
                sep + str(s.contig_start) + '..' + str(s.contig_end) + \
                ' motif=' + s.contig_motif + \
                ' on-target=' + str(round(s.azimuth, 2)) + \
                ' off-target=' + str(round(s.off_targets['hsuzhang'], 2)) + \
                ' alignments=' + str(len(s.alignments)), file=flo)
            print(s.contig_sequence, file=flo)
    print('Excision spacers FASTA generated: {!r}'.format(filename), file=sys.stderr)
    return filename

def generate_reversion_spacers(filename, sequences, sep=':'):
    """
    Creates a FASTA file
    >id:orientation:start..end:feature:contig motif=STR azimuth=N off-target=N alignments=N
    """
    with open(filename, 'w') as flo:
        for s in sequences:
            #feature, contig, orientation, start, end, seq, motif, azimuth, off_target, alignments = line
            #dDNA_id, feature, contig, orientation, start, end, seq, side, spacer, pam = line
            print(">" + s.dDNA_id + sep + s.contig_orientation + \
                sep + str(s.contig_start) + '..' + str(s.contig_end) + \
                sep + s.feature + sep + s.contig + \
                ' motif=' + s.contig_motif + \
                ' on-target=' + str(round(s.azimuth,2)) + \
                ' off-target=' + str(round(s.off_target_hsuzhang,2)) + \
                ' alignments=' + str(len(s.alignments)), file=flo)
            print(s.contig_sequence, file=flo)
    print('Excision spacers FASTA generated: {!r}'.format(filename), file=sys.stderr)
    return filename

def generate_excision_query(filename, sequences, sep=':'):
    """
    Creates a FASTA file
    >feature:contig:orientation:start..end
    """
    # Create the query, converting list of sequences into a FASTA file
    with open(filename, 'w') as flo:
        for line in sequences:
            feature, contig, orientation, start, end, seq, side, spacer, pam = line
            #print(">" + sep.join(map(str, [feature, contig, orientation, start, end])), file=flo)
            print(">" + feature + sep + contig + sep + orientation + sep + str(start) + '..' + str(end), file=flo)
            print(seq, file=flo)
            #print('+')
            #print('9'*len(seq))
    print('Excision query FASTA generated: {!r}'.format(filename), file=sys.stderr)
    return filename

def generate_reversion_query(filename, sequences, sep=':'):
    """
    Creates a FASTA file
    >id:orientation:start..end:feature:contig
    """
    with open(filename, 'w') as flo:
        for line in sequences:
            #rev_entry = tuple(["dDNA-"+str(i), feature, contig] + (orientation, start, end, filt_seq, side, filt_spacer, filt_pam))
            dDNA_id, feature, contig, orientation, start, end, seq, side, spacer, pam = line
            print(">" + dDNA_id + sep + orientation + sep + str(start) + '..' + str(end) + sep + feature + sep + contig, file=flo)
            print(seq, file=flo)
    print('Reversion query FASTA generated: {!r}'.format(filename), file=sys.stderr)
    return filename

def generate_donor(filename, sequences, sep=':'):
    """
    Creates a FASTA file
    >id feature:contig:orientation:start1..end1:mAT:start2..end2
    """
    with open(filename, 'w') as flo:
        for line in sequences:
            #don_entry = tuple(["dDNA-"+str(i), feature, contig, '+', start1, start2, end1, end2, dDNA])
            dDNA_id, feature, contig, orientation, start1, end1, mAT, start2, end2, seq = line
            print(">" + dDNA_id + ' ' + feature + sep + contig + sep + orientation + sep + str(start1) + '..' + str(end1) + sep + mAT + sep + str(start2) + sep + str(end2), file=flo)
            print(seq, file=flo)
    print('dDNA FASTA generated: {!r}'.format(filename), file=sys.stderr)
    return filename

def test():
    #print(load_git_version())
    alignments = load_sam_file(sys.argv[1])
    for a in alignments:
        for b in sorted(alignments[a]):
            print(a, b, alignments[a][b])

if (__name__ == "__main__"):
    test()