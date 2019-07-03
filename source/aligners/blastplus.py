#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/aligners/blastplus.py

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

# Requires BLAST+ >= 2.6.0
#  because it has SAM format output
#  $ blastn -query query.fasta -subject subject.fasta -outfmt 17 > output.sam

class BlastPlus(Aligner):
    def __init__(self):
        super().__init__("blastn", "Altschul, et al.", 1990,
            citation="Altschul, et al. Basic local alignment search tool. Journal of Molecular Biology 215, 403-410 (1990).",
            input='fasta',
            output='blastn',
            truncated=False,
            classification='pairwise'
        )
    def index(self, fasta, output_filename, output_folder, threads, *args, **kwargs):
        """
        Call makeblastdb on non-compressed FASTA file.
        
        The indexed FASTA will be stored in 'folder'.
        The 'fasta' basename will be used as both 'folder' and index name.
        """
        # Find the FASTA basename
        name = os.path.splitext(os.path.basename(output_filename))[0]
        
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
        
    def align(self, query, subject, output_filename, output_folder, threads, *args, **kwargs):
        """
        Aligns sequences using BLAST+
        
        'query' is the FASTA to align
        'subject' is the BLAST+ nucleotide index prefix
        
        Returns the path of the '*.blastn' file generated
        """
        output_filename_path = os.path.join(output_folder, output_filename)
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
    
    def load_record(self, flo, *args, **kwargs):
        """
        Loads the next record of the *.blastn file.
        
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
        
        # Define a non-existing record that will be replaced
        record = None
        
        # Loop through lines in SAM file until a valid, complete record is found
        record_found = False
        while (record_found == False):
            try:
                line = next(flo)
                sline = line.rstrip().split("\t")
                #if ((len(sline) > 5) and (sline[2] != '*')):
                record = sline
                record_found = True
            except StopIteration:
                break
        
        # Process the record if found
        if record:
            
            flags = 0
            smin, smax = sorted(map(int, record[8:10]))
            smin -= 1
            if (int(record[8]) > int(record[9])):
                flags |= 16
            
            record = Record(
                record[0], record[1], # query_name, subject_name,
                None, None, # query_sequence, subject_sequence,
                (int(record[6])-1, int(record[7])), (smin, smax), # query_position, subject_position,
                abs(int(record[7])-int(record[6])), abs(int(record[9])-int(record[8])), # query_length, subject_length,
                flags, None, float(record[11]), float(record[10]), int(record[3]) # flags, cigar, score, evalue, length
            )
        
        # Return the processed record, otherwise None if the file is complete
        return record

def test():
    """Code to test the BlastPlus Aligner subclass"""
    if (len(sys.argv[1:]) != 3):
        print("USAGE python3 blastplus.py query.fasta subject.fasta folder", file=sys.stderr)
        sys.exit(1)
    query = sys.argv[1]
    subject = sys.argv[2]
    folder = sys.argv[3]
    
    os.makedirs(folder, exist_ok=True)
    logging.basicConfig(filename=os.path.join(folder, 'log.txt'), level=logging.INFO, format='%(message)s') # format='%(levelname)s %(asctime)s: %(message)s'
    print("=== Blast+ ===")
    A = BlastPlus()
    index_path = A.index(subject, 'myindex', folder, 4)
    alignment_path = A.align(query, index_path, 'myoutput.blastn', folder, 2)
    print(index_path)
    print(alignment_path)
    
    with open(alignment_path, 'r') as flo:
        record = True
        while (record != None):
            record = A.load_record(flo)
            print(record)

if (__name__ == '__main__'):
    test()
