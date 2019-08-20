#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/utils.py

# Import standard packages
import sys
import os
import gzip
import datetime
import logging
import subprocess
import textwrap
import time

# Import non-standard packages
import regex

logger = logging.getLogger(__name__)

class GenBankFile(object):
    def __init__(self, name, sequence, molecule='ds-DNA', shape='linear'):
        self.name = name
        self.sequence = sequence
        self.length = str(len(sequence)) + ' bp'
        self.molecule = molecule
        self.shape = shape
        
        now = time.localtime()
        months = ['UNK', 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
        self.date = '{}-{}-{}'.format(now.tm_mday, months[now.tm_mon], now.tm_year)
        self.annotations = []
        
    def add_annotation(self, feature, start, end, strand, **kwargs):
        """
        features:
          assembly_gap, C_region, CDS, centromere, D-loop, D_segment, exon,
          gap, gene, iDNA, intron, J_segment, mat_peptide, misc_binding,
          misc_difference, misc_feature, misc_recomb, misc_RNA, misc_structure,
          mobile_element, modified_base, mRNA, ncRNA, N_region, old_sequence,
          operon, oriT, polyA_site, precursor_RNA, prim_transcript, primer_bind,
          propeptide, protein_bind, regulatory, repeat_region, rep_origin, rRNA,
          S_region, sig_peptide, source, stem_loop, STS, telomere, tmRNA,
          transit_peptide, tRNA, unsure, V_region, V_segment, variation,
          3'UTR, 5'UTR
        
        start < end
        """
        self.annotations.append([feature, start, end, strand, kwargs])
    def write(self, filename):
        with open(filename, 'w') as flo:
            print('LOCUS       ' + self.name.ljust(24, ' ') + ' ' + self.length + ' ' + self.molecule + '     ' + self.shape + '     ' + self.date, file=flo)
            print('DEFINITION  .', file=flo)
            print('FEATURES             Location/Qualifiers', file=flo)
            for a in self.annotations:
                feature, start, end, strand, details = a
                location = str(start+1) + '..' + str(end) # convert to 1-based index
                if (strand == '-'):
                    location = 'complement(' + location + ')'
                print('     ' + feature.ljust(16, ' ') + location, file=flo)
                for dkey, dvalue in details.items():
                    if isinstance(dvalue, str):
                        print('                     /' + dkey + '="' + dvalue + '"', file=flo)
                    else:
                        print('                     /' + dkey + '=' + str(dvalue), file=flo)
                #print('                     /label="' + label + '"', file=flo)
                #print('                     /ApEinfo_revcolor='+color, file=flo)
                #print('                     /ApEinfo_fwdcolor='+color, file=flo)
            print('ORIGIN', file=flo)
            n = 10
            m = 6
            sep_seqs = [self.sequence[i:i+n] for i in range(0, len(self.sequence), n)]
            for i in range(0, len(sep_seqs), m):
                print(str(i*n+1).rjust(9) + ' ' + ' '.join(sep_seqs[i:i+m]), file=flo)
                    
            print('//', file=flo)
        

def flatten(iterable, remove_none=False):
    """Make a flat list out of a list of lists"""
    # Will remove None
    if remove_none:
        return [item for sublist in iterable for item in sublist if item != None]
    else:
        return [item for sublist in iterable for item in sublist]

def lr_justify(left_text, right_text, width=140):
    if (len(left_text)+len(right_text) < width):
        text = ('{:>'+str(width)+'}').format(right_text)
        return left_text + text[len(left_text):]
    else:
        half = width//2 - 1
        if (len(right_text) > half):
            return left_text[:half] + '.'*(width-half-half) + right_text[-half:]
        else:
            return left_text[:width-len(right_text)-3] + '.'*3 + right_text
    # Some input arguments that cause this to give unexpected output:
    #  Not enough '...' (should be 3)
    #  The right-most text isn't justified
    # >>> lr_justify("", "this is my favorite thing to do in the entire whole wide world", 59)
    # '... the entire whole wide world'
    # >>> lr_justify("this is my favorite thing to do in the entire whole wide world", "something else smells like fish", 60)
    # 'this is my favorite thing to ..mething else smells like fish'

def lr_justify2(left_text, right_text, width=140, colsplit=None):
    if (colsplit == None):
        colsplit = width//2
    else:
        colsplit = max(0, min(colsplit, width))
    
    rts = [' ']*width
    if (len(left_text) + len(right_text) < width):
        for i, l in enumerate(left_text):
            rts[i] = l
        for i, r in enumerate(right_text[::-1]):
            rts[-(i+1)] = r
    else:
        i = 0
        for l in left_text:
            if (i < colsplit):
                rts[i] = l
            else:
                break
            i += 1
        i += 1
        rts[i] = '.'
        i += 1
        rts[i] = '.'
        i += 1
        rts[i] = '.'
        i += 1
        for j, r in enumerate(right_text[::-1]):
            if (width-j > i):
                rts[-(j+1)] = r
            else:
                break
    return ''.join(rts)

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
    
    #m = regex.findall(r'(\d+)[MISP=X]', cigar) # This was giving incorrect lengths...
    m = regex.findall(r'(\d+)[MDSP=X]', cigar) # This was giving correct lengths, even though it violates spec
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

def write_merged_fasta(contig_sequences, filename, columns=80):
    """
    Creates a FASTA file
    """
    wrapper = textwrap.TextWrapper(width=columns, replace_whitespace=False, drop_whitespace=False, expand_tabs=False, break_on_hyphens=False)
    with open(filename, 'w') as flo:
        for name, sequence in contig_sequences.items():
            print('>'+name, file=flo)
            print(wrapper.fill(sequence), file=flo)
    logger.info('Merged FASTA written: {!r}'.format(filename))
    return filename

def load_indexed_fasta_files(filenames):
    """
    Load multiple FASTA files
    Does not require primary sequence headers to be unique
    """
    
    # All input FASTA files are indexed
    # All input FASTA contigs are indexed:
    #
    # For instance, if program is run with the following parameters:
    #  --fasta F1.fasta F2.fasta
    # 
    # fasta_index = {
    # # N: FILENAME
    #   0: 'F1.fasta',
    #   1: 'F2.fasta',
    # }
    # contig_index = {
    # # N: (FASTA_INDEX, CONTIG, HEADER), <-- this tuple structure not implemented yet
    #   0: (0, 'contig1', '>contig1 something'),
    #   1: (0, 'contig2', '>contig2'),
    #   2: (1, 'contig1', '>contig1 nothing'),
    # }
    # 
    # contigs = {
    # # N: SEQUENCE
    #   0: 'ACGACTA...',
    #   1: 'CCAACCA...',
    #   2: 'GGTGTAG...',
    # }
    fasta_index = {} # {N: FILENAME, ...}
    contig_index = {} # {N: CONTIG, ...}
    contig_sequences = {} # {N: SEQUENCE, ...}
    for filename in filenames:
        ci, cs = load_fasta_file(filename)
        shift = len(contig_index)
        for k, v in ci.items():
            contig_index[k+shift] = v
        for k, v in cs.items():
            contig_sequences[k+shift] = v
        fasta_index[len(fasta_index)] = filename
    return fasta_index, contig_index, contig_sequences

def load_fasta_files_into_list(filenames):
    rounds = []
    
    for filename in filenames:
        ci, cs = load_fasta_file(filename)
        
        seqs = {}
        for index, name in ci.items():
            seqs[name] = cs[index]
        rounds.append(seqs)
    
    return rounds

def load_multiple_fasta_files(filenames):
    """
    Load multiple FASTA files
    Requires primary sequence headers to be unique
    """
    contig_index = {} # {N: CONTIG, ...}
    contig_sequences = {} # {CONTIG: SEQUENCE, ...}
    for filename in filenames:
        ci, cs = load_fasta_file(filename)
        shift = len(contig_index)
        
        for index, name in ci.items():
            contig_index[index+shift] = name
            contig_sequences[name] = cs[index]
    
    if (len(contig_sequences) != len(contig_index)):
        raise Exception('Duplicate contig names found in input FASTA files')
    
    return contig_sequences

def check_header_collisions(contig_index):
    """
    Function to warn the user if any contigs have the same name.
    """
    if (len(set(contig_index.values())) != len(contig_index.values())):
        logger.info('Duplicate contig names found in input FASTA files')

def load_fasta_file(filename):
    """
    Load contig sequences from file.
    Can open a *.gz compressed file.
    Returns an index dict() as well as sequence dict()
    """
    contig_index = {}
    contig_sequences = {}
    
    # Automatically process GZIP compressed files
    if filename.endswith('.gz'):
        flo = gzip.open(filename, 'rt')
    else:
        flo = open(filename, 'r')
    
    with flo:
        index = 0
        name = None
        seq = None
        for line in flo:
            line = line.rstrip()
            if line.startswith('>'):
                if (name != None):
                    contig_index[index] = name
                    contig_sequences[index] = seq
                    index += 1
                name = regex.split(r'\s+', line[1:], 1)[0]
                seq = ''
            else:
                seq += line
        if (name != None):
            contig_index[index] = name
            contig_sequences[index] = seq
    
    logger.info('FASTA file parsed: {!r}'.format(filename))
    
    return contig_index, contig_sequences

def old_load_fasta_file(filename):
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

def encode_text(text):
    '''
    GFF3 files are nine-column, tab-delimited, plain text files. Literal use
    of tab, newline, carriage return, the percent (%) sign, and control
    characters must be encoded using RFC 3986 Percent-Encoding; no other
    characters may be encoded. Backslash and other ad-hoc escaping
    conventions that have been added to the GFF format are not allowed.
    The file contents may include any character in the set supported by the
    operating environment, although for portability with other systems, use
    of Latin-1 or Unicode are recommended.
    
    Note that unescaped spaces are allowed within fields, meaning that
    parsers must split on tabs, not spaces. Use of the "+" (plus) character
    to encode spaces is depracated from early versions of the spec and is
    no longer allowed.
    
    Undefined fields are replaced with the "." character, as described in
    the original GFF spec.
    '''
    
    table = {
        '"': '%22',
        '%': '%25',
        '&': '%26',
        "'": '%27',
        '(': '%28',
        ')': '%29',
        ',': '%2C',
        ';': '%3B',
        '=': '%3D',
        '[': '%5B',
        '\\': '%5C',
        ']': '%5D',
        chr(127): '%7F',
    }
    for i in range(33):
        table[chr(i)] = '%' + hex(i)[2:].zfill(2)
    #table[chr(0)] = '%20'
    
    return ''.join(table.get(c,c) for c in text)

def decode_text(text):
    '''
    Replaces %-escaped hexadecimal text with its proper character.
    '''
    change = lambda x: chr(int('0x'+x.group()[1:], 16))
    new_text, number_substitutions = regex.subfn(r'(%[0-9A-F]{2})', change, text)
    return new_text

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

def load_git_revision():
    '''Returns the number of git revisions for the current head and master'''
    command_list = ['git', 'rev-list', 'HEAD']
    working_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
    try:
        cp = subprocess.run(command_list, cwd=working_dir, shell=False, check=True, stdout=subprocess.PIPE)
        return len(cp.stdout.decode().splitlines())
    except (FileNotFoundError, subprocess.CalledProcessError):
        return 'missing'

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

def print_local_file(filename):
    filepath = os.path.join(os.path.dirname(__file__), filename)
    with open(filepath, 'r') as flo:
        for line in flo:
            print(line.rstrip())

def which(cmd, mode=os.F_OK|os.X_OK, path=None, full=False):
    """Given a command, mode, and a PATH string, return the path which
    conforms to the given mode on the PATH, or None if there is no such
    file.
    `mode` defaults to os.F_OK | os.X_OK. `path` defaults to the result
    of os.environ.get("PATH"), or can be overridden with a custom search
    path. When `full` is True, then all instances of cmd on the PATH
    will be returned as a list.
    
    (Modified from Python source)
    """
    # Check that a given file can be accessed with the correct mode.
    # Additionally check that `file` is not a directory, as on Windows
    # directories pass the os.access check.
    def _access_check(fn, mode):
        return (os.path.exists(fn) and os.access(fn, mode)
                and not os.path.isdir(fn))

    # If we're given a path with a directory part, look it up directly rather
    # than referring to PATH directories. This includes checking relative to the
    # current directory, e.g. ./script
    if os.path.dirname(cmd):
        if _access_check(cmd, mode):
            return cmd
        return None

    if path is None:
        path = os.environ.get("PATH", os.defpath)
    if not path:
        return None
    path = path.split(os.pathsep)

    if sys.platform == "win32":
        # The current directory takes precedence on Windows.
        if not os.curdir in path:
            path.insert(0, os.curdir)

        # PATHEXT is necessary to check on Windows.
        pathext = os.environ.get("PATHEXT", "").split(os.pathsep)
        # See if the given file matches any of the expected path extensions.
        # This will allow us to short circuit when given "python.exe".
        # If it does match, only test that one, otherwise we have to try
        # others.
        if any(cmd.lower().endswith(ext.lower()) for ext in pathext):
            files = [cmd]
        else:
            files = [cmd + ext for ext in pathext]
    else:
        # On other platforms you don't have things like PATHEXT to tell you
        # what file suffixes are executable, so just pass on cmd as-is.
        files = [cmd]

    seen = set()
    rets = []
    for dir in path:
        normdir = os.path.normcase(dir)
        if not normdir in seen:
            seen.add(normdir)
            for thefile in files:
                name = os.path.join(dir, thefile)
                if full:
                    if _access_check(name, mode):
                        rets.append(name)
                else:
                    if _access_check(name, mode):
                        return name
    if (len(rets) > 0):
        return rets
    
    return None

def test():
    #print(load_git_version())
    alignments = load_sam_file(sys.argv[1])
    for a in alignments:
        for b in sorted(alignments[a]):
            print(a, b, alignments[a][b])

if (__name__ == "__main__"):
    test()