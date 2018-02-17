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
    from aligner import Aligner
else:
    from .aligner import Aligner

# Treat modules in PACKAGE_PARENT as in working directory
if (__name__ == "__main__"):
    # Relative path for package to import
    PACKAGE_PARENT = '..'
    # Obtain path of currently-running file
    SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
    # Convert to absolute path, and add to the PYTHONPATH
    sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))
    
    from utils import flatten
else:
    from ..utils import flatten

class Bowtie2(Aligner):
    def __init__(self):
        super().__init__("bowtie2", "Langmead & Salzberg", 2012,
            citation="Langmead & Salzberg. Fast gapped-read alignment with Bowtie 2. Nature Methods 9, 357-359 (2012).",
            input='fasta',
            output='sam',
            truncated=False,
        )
    
    def index(self, fasta, output_filename, output_folder, threads, *args, **kwargs):
        """
        Call bowtie2-build on non-compressed FASTA file.
        
        The indexed FASTA will be stored in 'folder'.
        The 'fasta' basename will be used as both 'folder' and index name.
        """
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
        
        return self.process('bowtie2-build', index_file, options)        
    
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
        ])
        
        return self.process('bowtie2', output_filename_path, options)

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

if (__name__ == '__main__'):
    test()
