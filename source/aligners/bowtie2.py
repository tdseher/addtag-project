#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/aligners/bowtie2.py

# List general Python imports
import sys
import os
import subprocess
from collections import OrderedDict
import logging
logger = logging.getLogger(__name__)

# import non-standard package
import regex

# import AddTag-specific packages
#from .. import utils

if (__name__ == "__main__"):
    from aligner import Aligner, Record
else:
    from .aligner import Aligner, Record

# Treat modules in PACKAGE_PARENT as in working directory
if (__name__ == "__main__"):
    # Relative path for package to import
    PACKAGE_PARENT = '..'
    # Obtain path of currently-running file
    SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
    # Convert to absolute path, and add to the PYTHONPATH
    sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))
    
    from cigarstrings import sam_orientation, cigar2query_position, cigar2query_aligned_length, cigar2subject_aligned_length, cigar2score
    from utils import old_load_fasta_file
    from evalues import EstimateVariables, load_scores
else:
    from ..cigarstrings import sam_orientation, cigar2query_position, cigar2query_aligned_length, cigar2subject_aligned_length, cigar2score
    from ..utils import old_load_fasta_file
    from ..evalues import EstimateVariables, load_scores

class Bowtie2(Aligner):
    logger = logger.getChild(__qualname__)
    
    def __init__(self):
        super().__init__("bowtie2", "Langmead & Salzberg", 2012,
            citation="Langmead & Salzberg. Fast gapped-read alignment with Bowtie 2. Nature Methods 9, 357-359 (2012).",
            input='fasta',
            output='sam',
            truncated=False,
            classification='pairwise'
        )
        self.score_matrix = {}
        self.ev = {}
        self.current_file = None
    
    def index(self, fasta, output_filename, output_folder, threads, *args, **kwargs):
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
        name = os.path.splitext(os.path.basename(output_filename))[0]
        
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
        outpath = self.process('bowtie2-build', index_file, options)
        
        self.ev[outpath] = ev
        self.score_matrix[outpath] = score_matrix
        
        return outpath
    
    def align(self, query, subject, output_filename, output_folder, threads, *args, **kwargs):
        """
        Align the query file to the subject file.
        
        query is the FASTA to align
        subject is the bowtie2 index
        
        Returns the path of the SAM file generated
        """
        output_filename_path = os.path.join(output_folder, output_filename)
        options = OrderedDict([
            ('-p', threads), # Number of processors to use
            ('-S', output_filename_path), # sam file,
            ('-x', subject), # Path to the index prefix (excludes file extensions)
            ('-U', query),
            ('-k', 25), # specifies the maximum number of alignments per sequence to return
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
        
        outpath = self.process('bowtie2', output_filename_path, options)
        
        self.ev[outpath] = self.ev[subject]
        self.score_matrix[outpath] = self.score_matrix[subject]
        
        return outpath
    
    #def load_file(self, filename, *args, **kwargs):
    #    """
    #    Function to prepare an iterator for the SAM file
    #    """
    #    # Each line is a separate record
    #    self.open_file = open(filename, 'r')
    
    def load_record(self, flo, *args, **kwargs):
        """
        Loads the next record of the SAM file.
        
        Since SAM files have one record per line, and one line per record,
        Just iterate until a valid line is found, then parse it.
        
        Returns:
         None if the file is complete
         Or a Record object
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
        
        # Define a non-existing record that will be replaced
        record = None
        
        # Loop through lines in SAM file until a valid, complete record is found
        record_found = False
        while (record_found == False):
            try:
                line = next(flo)
                if not line.startswith('@'):
                    sline = line.rstrip().split("\t")
                    if ((len(sline) > 5) and (sline[2] != '*')):
                        record = sline
                        record_found = True
            except StopIteration:
                break
        
        # Process the record if found
        if record:
            ev = self.ev[self.current_file]
            score_matrix = self.score_matrix[self.current_file]
            cigar_score = cigar2score(record[5], score_matrix)
            evalue = ev.calculate_evalue(cigar_score, record[9])
            record = Record(
                record[0], record[2], # query_name, subject_name,
                record[9], None, # query_sequence, subject_sequence,
                cigar2query_position(record[5]), (int(record[3])-1, int(record[3])-1+cigar2subject_aligned_length(record[5])), # query_position, subject_position,
                cigar2query_aligned_length(record[5]), cigar2subject_aligned_length(record[5]), # query_length, subject_length,
                int(record[1]), record[5], float(record[4]), evalue, None # flags, cigar, score, evalue, length
            )
        
        # Return the processed record, otherwise None if the file is complete
        return record

def test():
    """Code to test the Bowtie2 Aligner subclass"""
    if (len(sys.argv[1:]) != 3):
        print("USAGE python3 bowtie2.py query.fasta subject.fasta folder", file=sys.stderr)
        sys.exit(1)
    query = sys.argv[1]
    subject = sys.argv[2]
    folder = sys.argv[3]
    
    os.makedirs(folder, exist_ok=True)
    logging.basicConfig(filename=os.path.join(folder, 'log.txt'), level=logging.INFO, format='%(message)s') # format='%(levelname)s %(asctime)s: %(message)s'
    print("=== Bowtie 2 ===")
    A = Bowtie2()
    index_path = A.index(subject, 'myindex', folder, 4)
    alignment_path = A.align(query, index_path, 'myoutput.sam', folder, 16)
    print(index_path)
    print(alignment_path)
    
    with open(alignment_path, 'r') as flo:
        record = True
        while (record != None):
            record = A.load_record(flo)
            print(record)

if (__name__ == '__main__'):
    test()
