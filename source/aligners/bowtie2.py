#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/aligners/bowtie2.py

# List general Python imports
import sys
import os
from collections import OrderedDict
import logging
logger = logging.getLogger(__name__)

# import AddTag-specific packages
from .aligner import PairwiseAligner, Record
from ..cigarstrings import sam_orientation, cigar2query_position, cigar2query_aligned_length, cigar2subject_aligned_length, cigar2score
from ..utils import old_load_fasta_file, which
from ..evalues import EstimateVariables, load_scores

class Bowtie2(PairwiseAligner):
    logger = logger.getChild(__qualname__)
    
    def __init__(self):
        super().__init__(
            name="bowtie2",
            authors=['Langmead, Ben', 'Salzberg, Steven L.'],
            title='Fast gapped-read alignment with Bowtie 2',
            journal='Nature Methods',
            issuing='9:357-359',
            year=2012,
            doi='https://doi.org/10.1038/nmeth.1923',
            #citation="Langmead & Salzberg. Fast gapped-read alignment with Bowtie 2. Nature Methods 9, 357-359 (2012).",
            input='fasta',
            output='sam',
            truncated=False,
            #classification='pairwise'
        )
        self.score_matrix = {}
        self.ev = {}
        self.binaries = self.get_binaries(['bowtie2', 'bowtie2-build'], full=False)

    def is_available(self):
        if all([which(x) for x in ['bowtie2', 'bowtie2-build']]):
            return True
        elif (sys.platform.startswith('win') and all([which(x) for x in ['bowtie2.bat', 'bowtie2-build.bat']])):
            return True
        else:
            return False
    
    def index(self, fasta, output_prefix, output_folder, threads, *args, **kwargs):
        """
        Call bowtie2-build on non-compressed FASTA file.
        
        The indexed FASTA will be stored in 'folder'.
        The 'fasta' basename will be used as both 'folder' and index name.
        
        Additionally, compute the Karlin-Altschul statistics parameters
        necessary for E-value calculations
        """
        
        # Karlin-Altschul calculations
        score_matrix = load_scores(os.path.join(os.path.dirname(__file__), 'bowtie2_scores.txt'))
        scontigs = old_load_fasta_file(fasta)
        ev = EstimateVariables(score_matrix, scontigs)
        
        # threads=(os.cpu_count() or 1)
        
        # Find the FASTA basename
        #name = os.path.splitext(os.path.basename(output_filename))[0]
        name = output_prefix
        
        # Make the directory if it does not yet exist
        #try:
        #    os.mkdir(os.path.join(output_folder, name))
        #except FileExistsError:
        #    pass
        os.makedirs(os.path.join(output_folder, name), exist_ok=True)
        
        index_file = os.path.join(output_folder, name, name)
        
        options = OrderedDict([
            ('--threads', threads),
            ('--offrate', 1),
            (fasta, None),
            (index_file, None),
        ])
        outpath = self.process(self.binaries[1], index_file, options) # 'bowtie2-build'
        
        self.ev[outpath] = ev
        self.score_matrix[outpath] = score_matrix
        
        return outpath

    def align(self, query, subject, output_prefix, output_folder, threads, *args, **kwargs):
        """
        Align the query file to the subject file.
        
        query is the FASTA to align
        subject is the bowtie2 index
        
        Returns the path of the SAM file generated
        """
        os.environ['BOWTIE2_INDEXES'] = os.path.dirname(subject)
        output_filename_path = os.path.join(output_folder, output_prefix + '.' + self.output)
        options = OrderedDict([
            ('-p', threads), # Number of processors to use
            ('-S', output_filename_path), # sam file,
            ('-x', os.path.basename(subject)), # Path to the index prefix (excludes file extensions)
            ('-U', query),
            ('-k', 100), # specifies the maximum number of alignments per sequence to return
            ('-N', 1), # Sets the number of mismatches to allowed in a seed alignment
                       # during multiseed alignment. Can be set to 0 or 1. Setting
                       # this higher makes alignment slower (often much slower) but
                       # increases sensitivity.
            ('-L', 10), # Sets the length of the seed substrings to align during
                        # multiseed alignment. Smaller values make alignment slower
                        # but more sensitive.
            ('-D', 30), # Up to <int> consecutive seed extension attempts can "fail"
                        # before Bowtie 2 moves on, using the alignments found so far.
                        # A seed extension "fails" if it does not yield a new best or
                        # a new second-best alignment. This limit is automatically
                        # adjusted up when -k or -a are specified. Default: 15.'
            ('-R', 3), # the maximum number of times Bowtie 2 will "re-seed" reads
                       # with repetitive seeds. When "re-seeding," Bowtie 2 simply
                       # chooses a new set of reads (same length, same number of
                       # mismatches allowed) at different offsets and searches for
                       # more alignments. A read is considered to have repetitive
                       # seeds if the total number of seed hits divided by the
                       # number of seeds that aligned at least once is greater
                       # than 300. Default: 2.
            ('-i', 'S,1,0'), # split query read into seed index every 1 bp.
            ('--n-ceil', 'L,3,0'), # Maximum number of Ns is 3
            ('-f', None), # query reads are FASTA
            ('--end-to-end', None), # disallow soft masking
            ('--mp', '3,2'), # Set mismatch penalty to 3 Default: '6,2'
            ('--rdg', '4,2'), # Sets the read gap open (<int1>) and extend (<int2>)
                              # penalties. A read gap of length N gets a penalty of
                              # <int1> + N * <int2>. Default: '5,3'.
            ('--rfg', '4,2'), # Sets the reference gap open (<int1>) and extend (<int2>)
                              # penalties. A reference gap of length N gets a penalty of
                              # <int1> + N * <int2>. Default: '5,3'.
            ('--score-min', 'L,-0.9,-0.9'), # Sets function governing minimum alignment
                                            # score needed for an alignment to be considered
                                            # good enough to report, where 'read length' = L:
                                            #   score(L, a, b) = a + b * L
                                            # Default: 'L,-0.6,-0.6'
                                            # Switching to 'L,-0.9,-0.9' returns roughly 9 times
                                            # more alignments than the default
            ('--xeq', None), # Use '='/'X', instead of 'M', to specify matches/mismatches in CIGAR string
        ])
        
        outpath = self.process(self.binaries[0], output_filename_path, options) # 'bowtie2'
        
        self.ev[outpath] = self.ev[subject]
        self.score_matrix[outpath] = self.score_matrix[subject]
        
        return outpath
    
    def load(self, filename, *args, **kwargs):
        """
        Yields records from the input SAM file, one at a time.
        BOWTIE2 SAM files have one record per line, and one line per record.
        :param filename: Path of the file to create Records from
        :param args:
        :param kwargs: 
        :return: Yields either the next Record, or a StopIteration Exception if no more alignments
        """
        
        # Sequence Alignment/Map (SAM) format is TAB-delimited. Apart from the
        # header lines, which are started with the '@' symbol, each alignment line
        # consists of:
        #   Col  Field  Description
        #     0  QNAME  Query template/pair NAME
        #     1  FLAG   bitwise FLAG
        #     2  RNAME  Reference sequence NAME
        #     3  POS    1-based leftmost POSition/coordinate of clipped sequence
        #     4  MAPQ   MAPping Quality (Phred-scaled)
        #     5  CIAGR  extended CIGAR string
        #     6  MRNM   Mate Reference sequence NaMe ('=' if same as RNAME)
        #     7  MPOS   1-based Mate POSistion
        #     8  TLEN   inferred Template LENgth (insert size)
        #     9  SEQ    query SEQuence on the same strand as the reference
        #    10  QUAL   query QUALity (ASCII-33 gives the Phred base quality)
        #   11+  OPT    variable OPTional fields in the format TAG:VTYPE:VALUE
        
        # Code to decompress a *.bam file should go here
        
        with open(filename) as flo:
            for line in flo:
                # Ignore header lines in SAM file
                if not line.startswith('@'):
                    # Remove newline character from end of line
                    # Split SAM record str into a list
                    sline = line.rstrip().split("\t")
                    
                    # Require the line to have a record (not a blank line)
                    # Require the record to align to a contig
                    if ((len(sline) > 5) and (sline[2] != '*')):
                        yield self.create_record(sline, filename)
    
    def create_record(self, sline, filename):
        '''
        Converts a split line from a SAM file into a Record object
        :param sline: line.rstrip().split('\t')
        :return: a Record object
        '''
        # Find evalue of the alignment
        ev = self.ev[filename]
        score_matrix = self.score_matrix[filename]
        cigar_score = cigar2score(sline[5], score_matrix)
        evalue = ev.calculate_evalue(cigar_score, sline[9])
        
        # Build Record
        c2sal = cigar2subject_aligned_length(sline[5])
        record = Record(
            sline[0], sline[2], # query_name, subject_name,
            sline[9], None, # query_sequence, subject_sequence,
            cigar2query_position(sline[5]), (int(sline[3])-1, int(sline[3])-1+c2sal), # query_position, subject_position,
            cigar2query_aligned_length(sline[5]), c2sal, # query_length, subject_length,
            int(sline[1]), sline[5], float(sline[4]), evalue, None # flags, cigar, score, evalue, length
        )

        # Return the processed Record
        return record
