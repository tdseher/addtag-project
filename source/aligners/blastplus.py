#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/aligners/blastplus.py

# List general Python imports
import os
from collections import OrderedDict
import logging
logger = logging.getLogger(__name__)

# import AddTag-specific packages
from .aligner import Aligner, Record
from ..utils import which

# TODO: Consider looking into BLAST+ >= 2.6.0
#       because it has SAM format output
#       $ blastn -query query.fasta -subject subject.fasta -outfmt 17 > output.sam

class BlastPlus(Aligner):
    logger = logger.getChild(__qualname__)
    
    def __init__(self):
        super().__init__(
            name="blastn",
            authors=['Altschul, Stephen F.', 'Gish, Warren', 'Miller, Webb', 'Myers, Eugene W.', 'Lipman, David J.,'],
            title='Basic local alignment search tool',
            journal='Journal of Molecular Biology',
            issuing='215(3):403-410',
            year=1990,
            doi='https://doi.org/10.1016/S0022-2836(05)80360-2',
            #citation="Altschul, et al. Basic local alignment search tool. Journal of Molecular Biology 215, 403-410 (1990).",
            input='fasta',
            output='blastn',
            truncated=False,
            classification='pairwise'
        )

    def is_available(self):
        if all([which(x) for x in ['blastn', 'makeblastdb']]):
            return True
        else:
            return False

    def index(self, fasta, output_prefix, output_folder, threads, *args, **kwargs):
        """
        Call makeblastdb on non-compressed FASTA file.
        
        The indexed FASTA will be stored in 'folder'.
        The 'fasta' basename will be used as both 'folder' and index name.
        """
        # Find the FASTA basename
        #name = os.path.splitext(os.path.basename(output_filename))[0]
        name = output_prefix
        
        # Make the directory if it does not yet exist
        os.makedirs(os.path.join(output_folder, name), exist_ok=True)
        
        # Build the database output filename
        index_file = os.path.join(output_folder, name, name)
        
        options = OrderedDict([
            ('-dbtype', 'nucl'),
            ('-in', fasta),
            ('-hash_index', None),
            ('-out', index_file),
        ])
        
        return self.process('makeblastdb', index_file, options)
        
    def align(self, query, subject, output_prefix, output_folder, threads, *args, **kwargs):
        """
        Aligns sequences using BLAST+
        
        'query' is the FASTA to align
        'subject' is the BLAST+ nucleotide index prefix
        
        Returns the path of the '*.blastn' file generated
        """
        output_filename_path = os.path.join(output_folder, output_prefix + '.' + self.output)
        options = OrderedDict([
            ('-num_threads', threads), # Number of processors to use
            ('-query', query), # Path the the query FASTA file
            ('-db', subject), # Path to the index prefix (excludes file extensions)
            ('-outfmt', 6), # Use BLAST+ tabular output format
            ('-out', output_filename_path), # Path to the output '*.blastn' file to generate
            ('-word_size', 11), # Word size for wordfinder algorithm (length of best perfect match)
            ('-dust', 'no'), # Filter query sequence with DUST (Format: 'yes', 'level window linker', or 'no' to disable)
            #('-soft_masking', 'false') # Apply filtering locations as soft masks (i.e., only for finding initial matches).
        ])
        
        return self.process('blastn', output_filename_path, options)
    
    def load(self, filename, *args, **kwargs):
        """
        Yields the next record of the '*.blastn' file.
        
        Since '-outfmt 6' *.blastn files have one record per line, and one
        line per record, this just iterates until a valid line is found, then
        parses it.
        
        Returns:
         None if the file is complete
         Record object otherwise
        """
        
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
        
        with open(filename) as flo:
            for line in flo:
                # Ignore header/comment lines in BLASTN file
                if not line.startswith('#'):
                    # Remove newline character from end of line
                    # Split BLASTN record str into a list
                    sline = line.rstrip().split("\t")
                    
                    # Require the line to have a record (not a blank line)
                    if (len(sline) > 1):
                        yield self.create_record(sline)
    
    def create_record(self, sline):
        '''
        Converts a split line from a BLASTN file into a Record object
        :param sline: line.rstrip().split('\t')
        :return: a Record object
        '''
        
        # Find 'flags' of the alignment
        flags = 0
        smin, smax = sorted(map(int, sline[8:10]))
        smin -= 1
        if (int(sline[8]) > int(sline[9])):
            flags |= 16
        
        # Build sline
        record = Record(
            sline[0], sline[1], # query_name, subject_name,
            None, None, # query_sequence, subject_sequence,
            (int(sline[6])-1, int(sline[7])), (smin, smax), # query_position, subject_position,
            abs(int(sline[7])-int(sline[6])), abs(int(sline[9])-int(sline[8])), # query_length, subject_length,
            flags, None, float(sline[11]), float(sline[10]), int(sline[3]) # flags, cigar, score, evalue, length
        )
        
        # Return the processed record, otherwise None if the file is complete
        return record
