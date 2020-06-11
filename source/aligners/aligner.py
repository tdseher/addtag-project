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

# import AddTag-specific packages
#sys.path.append( os.path.join( os.path.dirname(__file__), os.path.pardir ) )
#from utils import flatten
from ..utils import flatten, which

class Aligner(object): # Name of the subclass
    """
    General class that should be subclassed when adding a new alignment program
    """
    logger = logger.getChild(__qualname__)
    
    def __init__(
        self,
        name,
        authors,
        title,
        journal,
        issuing,
        year,
        doi,
        #citation=None,
        input='fasta',
        output='sam',
        truncated=False,
        #classification='pairwise'
    ):
        """
        Specify general information regarding this new instance of the 
        Aligner class.
        """
        self.name = name                     # Unique name for the algorithm (str). No other Aligner objects should have this name.
        self.authors = authors               # Author of the algorithm (list of str objects)
        self.title = title
        self.journal = journal
        self.issuing = issuing
        self.year = year                     # Year algorithm published (int)
        self.doi = doi
        #self.citation = citation             # Citation (None, str)
        self.input = input                   # Format of expected input (STDIN, fasta, fastq, etc)
        self.output = output                 # Designate the output format the alignment will be in, so AddTag can select the correct parser
        self.truncated = truncated           # Normally, the entire short read is locally aligned to the genome contigs. This tells the parser that what is reported in the SAM file in the 10th column is incomplete.
        #self.classification = classification # Type ('high-throughput', 'pairwise', 'multiple-sequence', 'ht', 'pw', 'ms')
        self.available = self.is_available()

    def is_available(self):
        """
        Determines if the prerequisites for the Aligner have been met.
        For instance, it will query the PATH to determine if the alignment software is available.

        Overload this method
        :return: True or False
        """
        return False
    
    def get_binaries(self, binaries_list, full=False):
        paths = []
        
        if sys.platform.startswith('win'):
            for b in binaries_list:
                for b2 in [b, b+'.bat']:
                    p = which(b2)
                    if full:
                        paths.append(p)
                        break
                    elif p:
                        paths.append(b2)
                        break
                else:
                    paths.append(None)
        
        else:
            for b in binaries_list:
                p = which(b)
                if full:
                    paths.append(p)
                elif p:
                    paths.append(b)
                else:
                    paths.append(None)
        
        return paths
    
    def process(self, program, output_filename, fixed_options, *args, capture_stdout=False, append=False, **kwargs):
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

        if append:
            # Open for reading and writing.  The file is created if it does not
            # exist.  The stream is positioned at the end of the file.  Subse-
            # quent writes to the file will always end up at the then current
            # end of file, irrespective of any intervening fseek(3) or similar.
            mode = 'a+'
        else:
            mode = 'w+'

        # Perform alignment
        self.logger.info('command: {!r}'.format(command_str))
        if capture_stdout:
            with open(error_filename, mode) as err_flo, open(output_filename, 'w+') as out_flo:
                cp = subprocess.run(command_list, shell=False, check=True, stdout=out_flo, stderr=err_flo)
        else:
            with open(error_filename, mode) as flo:
                cp = subprocess.run(command_list, shell=False, check=True, stdout=flo, stderr=subprocess.STDOUT)
        
        self.logger.info('errors: {!r}'.format(error_filename))
        self.logger.info('output: {!r}'.format(output_filename))
        return output_filename
    
    def load(self, filename, *args, **kwargs):
        """
        Takes a file path as input.
        Yields either an alignment Record or a StopIteration Exception
        This should be a Generator function.
        
        Overload this method.
        """
        yield None

    # TODO: create 'Aligner.cleanup()' methods for each 'Aligner'
    def cleanup(self):
        """
        Delete intermediary files

        Overload this method
        """
        pass
    
    def __repr__(self):
        """
        Return the string representation of the Aligner
        """
        return self.__class__.__name__ + '(' + ', '.join(['name='+repr(self.name), 'authors='+repr(self.authors), 'year='+repr(self.year)]) + ')'

# class HighThroughputAligner(Aligner):
#     pass

class PairwiseAligner(Aligner):
    logger = logger.getChild(__qualname__)
    
    def index(self, fasta, output_prefix, output_folder, threads, *args, **kwargs):
        """
        Run the program process to index the subject.
        Returns the path of the index generated.
        
        Overload this method.
        """
        return ''

    def align(self, query, subject, output_prefix, output_folder, threads, *args, **kwargs):
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

class MultipleSequenceAligner(Aligner):
     logger = logger.getChild(__qualname__)
     
     def align(self, sequences_file, output_prefix, output_folder, threads, *args, **kwargs):
         '''
         Align the sequences within 'sequences_file', and write resultant MSA to 
         output_folder/output_prefix/output_extension. 
         :param sequences_file: 
         :param output_prefix: 
         :param output_folder: 
         :param threads: 
         :param args: 
         :param kwargs: 
         :return: path of MSA file generated
         '''
         
         return ''
     

class Record(object):
    __slots__ = [
        'query_name',
        'subject_name',
        
        'query_sequence',
        'subject_sequence',
        
        'query_position',
        'subject_position',
        
        'query_length',
        'subject_length',
        
        'flags',
        'cigar',
        'score',
        'evalue',
        'length',
    ]
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
            'query={}:{}..{}'.format(self.query_name, self.query_position[0], self.query_position[1]),
            'subject={}:{}..{}'.format(self.subject_name, self.subject_position[0], self.subject_position[1]),
            #'query='+self.query_name+':'+'..'.join(map(str, self.query_position)),
            #'subject='+self.subject_name+':'+'..'.join(map(str, self.subject_position)), 
            'cigar={}'.format(self.cigar),
            'score={}'.format(self.score),
            'evalue={}'.format(self.evalue)
            ]) + ')'
