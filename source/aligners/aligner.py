#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/aligner.py

# Import standard packages
import sys
import os
import shlex
import subprocess
from collections import OrderedDict
import logging
logger = logging.getLogger(__name__)

# import non-standard package
import regex

# import AddTag-specific packages
sys.path.append( os.path.join( os.path.dirname(__file__), os.path.pardir ) )
from utils import flatten

class Aligner(object): # Name of the subclass
    """
    General class that should be subclassed when adding a new alignment program
    """
    def __init__(self, 
        name,
        author,
        year,
        citation=None,
        input='fasta',
        output='sam',
        truncated=False,
        classification='pairwise'
    ):
        """
        Specify general information regarding this new instance of the 
        Aligner class.
        """
        self.name = name                     # Unique name for the algorithm (str). No other Aligner objects should have this name.
        self.author = author                 # Author of the algorithm (str)
        self.year = year                     # Year algorithm published (int)
        self.citation = citation             # Citation (None, str)
        self.input = input                   # Format of expected input (STDIN, fasta, fastq, etc)
        self.output = output                 # Designate the output format the alignment will be in, so AddTag can select the correct parser
        self.truncated = truncated           # Normally, the entire short read is locally aligned to the genome contigs. This tells the parser that what is reported in the SAM file in the 10th column is incomplete.
        self.classification = classification # Type ('high-throughput', 'pairwise', 'multiple-sequence', 'ht', 'pw', 'ms')
        
    def index(self, fasta, output_filename, output_folder, threads, *args, **kwargs):
        """
        Run the program process to index the subject.
        Returns the path of the index generated.
        
        Overload this method.
        """
        return ''
    
    def align(self, query, subject, output_filename, output_folder, threads, *args, **kwargs):
        """
        Align the query file to the subject file.
        Returns the path of the alignment generated (usually, the SAM file).
        
        Overload this method.
        """
        # Example 
        #options = OrderedDict([
        #    ('-q', query),
        #    ('-s', subject),
        #    ('-o', os.path.join(output_folder, output_filename)),
        #    ('-t', threads),
        #])
        #return self.process('align', 'test.sam', options)
        return ''
    
    def process(self, program, output_filename, fixed_options, capture_stdout=False, *args, **kwargs):
        """
        Align the query file to the subject file.
        
        fixed_options (str, list, or dict)
        user_options (str, list, or dict)
        
        Returns the path of the alignment generated (usually, the SAM file)
        """
        #user_options = None
        
        # name='' # the basename to give the generated files
        # threads=(os.cpu_count() or 1)
        # folder=os.getcwd() # generated files will be placed here
        # sep='' # the delimiter to use for generating the sequence headers
        
        #name = os.path.splitext(os.path.basename(input_query_file))[0]
        #sam_file = os.path.join(folder, name + '.sam')
        #error_file = os.path.join(folder, name + '.err')
        error_filename = os.path.splitext(output_filename)[0] + '.err'
        
        # if options is specified, then append to the list of options that will be called
        
        flat_fixed_options = []
        if fixed_options:
            if isinstance(fixed_options, str):
                flat_fixed_options = shlex.split(fixed_options)
            elif isinstance(fixed_options, list):
                flat_fixed_options = fixed_options
            elif isinstance(fixed_options, OrderedDict):
                flat_fixed_options = flatten(fixed_options.items(), remove_none=True)
            elif isinstance(fixed_options, dict):
                flat_fixed_options = flatten(fixed_options.items(), remove_none=True)
        
        #flat_user_options = []
        #if user_options:
        #    if isinstance(user_options, str):
        #        flat_user_options = shlex.split(fixed_options)
        #    elif isinstance(user_options, list):
        #        flat_user_options = user_options
        #    elif isinstance(user_options, OrderedDict):
        #        flat_user_options = flatten(user_options.items(), remove_none=True)
        #    elif isinstance(user_options, dict):
        #        flat_user_options = flatten(user_options.items(), remove_none=True)
        #
        ## Throw an exception if the user tries to set a reserved option
        #for o in flat_fixed_options:
        #    if o in flat_user_options:
        #        raise Exception("Command option {!r} is reserved by AddTag".format(o))
        #
        #command_list = [program] + flat_fixed_options + flat_user_options
        command_list = [program] + list(map(str, flat_fixed_options))
        command_str = ' '.join(command_list)
        
        # Perform alignment
        logger.info('command: {!r}'.format(command_str))
        if capture_stdout:
            with open(error_filename, 'w+') as err_flo:
                with open(output_filename, 'w+') as out_flo:
                    cp = subprocess.run(command_list, shell=False, check=True, stdout=out_flo, stderr=err_flo)
        else:
            with open(error_filename, 'w+') as flo:
                cp = subprocess.run(command_list, shell=False, check=True, stdout=flo, stderr=subprocess.STDOUT)
        
        logger.info('errors: {!r}'.format(error_filename))
        logger.info('output: {!r}'.format(output_filename))
        return output_filename
    
    def load_record(self, flo, *args, **kwargs):
        """
        Takes an open file like object as input.
        Outputs a Record object representing the alignment.
        
        Overload this method.
        """
        return None
    
    def load(self, filename, filetype, *args, **kwargs):
        """
        Parses alignment output and returns list of lists readable by AddTag
        [[contig, start, end, orientation, sequence], ...]
        """
        if (filetype == 'sam'):
            # SAM is TAB-delimited. Apart from the header lines, which are
            # started with the '@' symbol, each alignment line consists of:
            #
            # Col  Field  Description
            #   0  QNAME  Query (pair) NAME
            #   1  FLAG   bitwise FLAG
            #   2  RNAME  Reference sequence NAME
            #   3  POS    1-based leftmost POSition/coordinate of clipped sequence
            #   4  MAPQ   MAPping Quality (Phred-scaled)
            #   5  CIAGR  extended CIGAR string
            #   6  MRNM   Mate Reference sequence NaMe ('=' if same as RNAME)
            #   7  MPOS   1-based Mate POSistion
            #   8  ISIZE  Inferred insert SIZE
            #   9  SEQ    query SEQuence on the same strand as the reference
            #  10  QUAL   query QUALity (ASCII-33 gives the Phred base quality)
            #  11  OPT    variable OPTional fields in the format TAG:VTYPE:VALUE
            pass
        elif (filetype == 'blastn'):
            # Col  Field     Description
            #   0  qaccver   Query accesion.version
            #   1  saccver   Subject accession.version
            #   2  pident    Percentage of identical matches
            #   3  length    Alignment length
            #   4  mismatch  Number of mismatches
            #   5  gapopen   Number of gap openings
            #   6  qstart    Start of alignment in query
            #   7  qend      End of alignment in query
            #   8  sstart    Start of alignment in subject
            #   9  send      End of alignment in subject
            #  10  evalue    Expect value
            #  11  bitscore  Bit score
            pass
        elif (filetype == 'casoffinder'):
            pass
        
        return []
    
    def __repr__(self):
        """
        Return the string representation of the Aligner
        """
        return self.__class__.__name__ + '(' + ', '.join(['name='+repr(self.name), 'author='+repr(self.author), 'year='+repr(self.year)]) + ')'

class Record(object):
    def __init__(self,
            query_name, subject_name,
            query_sequence, subject_sequence,
            query_position, subject_position,
            query_length, subject_length,
            flags, cigar, score, evalue, length):
        
        self.query_name = query_name
        self.subject_name = subject_name
        
        self.query_sequence = query_sequence
        self.subject_sequence = subject_sequence
        
        self.query_position = query_position # (start, end)
        self.subject_position = subject_position # (start, end)
        
        self.query_length = query_length
        self.subject_length = subject_length
        
        self.flags = flags
        self.cigar = cigar
        self.score = score
        self.evalue = evalue
        self.length = length
        
    def __repr__(self):
        """Return string representation of the Record"""
        return self.__class__.__name__ + '(' + ', '.join([
            'query='+self.query_name+':'+'..'.join(map(str, self.query_position)),
            'subject='+self.subject_name+':'+'..'.join(map(str, self.subject_position))
            #'evalue='+str(self.evalue)
            ]) + ')'
