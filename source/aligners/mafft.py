#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/aligners/mafft.py

# List general Python imports
import sys
import os
from collections import OrderedDict
import logging
logger = logging.getLogger(__name__)

# import non-standard package
import regex

# import AddTag-specific packages
from .aligner import MultipleSequenceAligner, Record
from ..utils import which

class Mafft(MultipleSequenceAligner):
    def __init__(self):
        super().__init__(
            name="mafft",
            authors=['Toh, Hiroyuki', 'Katoh, Kazutaka', 'Kuma, Kei-ichi', 'Miyata, Takashi'],
            title='MAFFT version 5: improvement in accuracy of multiple sequence alignment',
            journal='Nucleic Acids Research',
            issuing='33(2):511-518',
            year=2005,
            doi='https://doi.org/10.1093/nar/gki198',
            #citation="Katoh, Kuma, Toh, & Miyata. MAFFT version 5: improvement in accuracy of multiple sequence alignment. Nucleic Acids Research 33: 2, 511-518 (2005).",
            input='fasta',
            output='fasta',
            truncated=False,
            #classification='multiple-sequence'
        )
        self.score_matrix = {}
        self.ev = {}
        self.current_file = None
        self.records = {}
        self.binaries = self.get_binaries(['mafft'], full=False)
    
    def is_available(self):
        """
        Determines if the prerequisites for the Aligner have been met.
        :return: True or False
        """
        if (which('mafft') or which('mafft.bat')):
            return True
        else:
            return False
    
    def align(self, query, output_filename, output_folder, threads, *args, **kwargs):
        """
        Perform multiple sequence alignment of all the sequences within the 'query' file.
        
        'query' is the FASTA file containing multiple sequences that will be aligned
        'subject' is ignored, and should be 'None'.
        
        Returns the path of the CLUSTAL file generated.
        """
        output_filename_path = os.path.join(output_folder, output_filename)
        
        # 'mafft --reorder --anysymbol --maxiterate 1000 --retree 1 --localpair input > output'
        
        # To output to fasta file, then run:
        #   mafft --inputorder --anysymbol --maxiterate 1000 --retree 2 --localpair --thread 4 mafft-input.fasta > mafft-output.fasta
        # To have the final line that summarizes amount of conservation
        #   mafft --clustalout --anysymbol --maxiterate 1000 --retree 2 --localpair --thread 4 mafft-input.fasta > mafft-output.aln
        
        options = OrderedDict([
            ('--inputorder', None),   # Describes output format and sequence ordering ('--reorder', '--clustalout', '--inputorder').
            ('--anysymbol', None),    # Upper/lower case is preserved.  The --anysymbol option is internally equivalent to the --preservecase option. Use unusual characters (e.g., 'U' as selenocysteine in protein sequence; 'i' as inosine in nucleotide sequence). Unusual characters are scored as unknown (not considered in the calculation).
            ('--nuc', None),          # Specify nucleotide (DNA/RNA) sequence type.
            ('--maxiterate', 1000),   # Number cycles of iterative refinement are performed.
            ('--retree', 2),          # Guide tree is built number times in the progressive stage. Valid with 6mer distance.
            ('--localpair', None),    # Use the 'L-INS-i' algorithm. Aligns multiple sequences around a single alignable domain. Flanking sequences are ignored in the pairwise alignment by the Smith-Waterman algorithm.
            ('--thread', threads),    # For all-to-all pairwise comparison, the number of threads specified with the --thread option are used.  The wall-clock time for this stage linearly decreases with the number of threads.
            #('--threadtb', threads), # For progressive alignment stage, the efficiency is not high (requiring large RAM) when using a large number (>>20) of threads. When using many threads for the progressive alignment stage, memory shortage error can occur.  So ts's better to unset this flag in most cases.
            #('--threadit', threads), # For iterative refinement stage, too, the efficiency is not high when using a large number (>>20) of threads.
            (query, None),            # FASTA file containing sequences to be aligned
        ])
        # output is sent to STDOUT
        
        outpath = self.process(self.binaries[0], output_filename_path, options, capture_stdout=True)
        
        # Karlin-Altschul calculations
        #self.ev[outpath] = self.ev[subject]
        #self.score_matrix[outpath] = self.score_matrix[subject]
        
        return outpath
    
    def load(self, filename, *args, **kwargs):
        if (self.output == 'fasta'):
            return self.load_fasta(filename)
        elif (self.output == 'clustal'):
            return self.load_clustal(filename)
    
    def load_fasta(self, filename):
        """
        Yields records from the input FASTA file, one at a time.
        :param filename: Path of the file to create Records from
        :param args:
        :param kwargs:
        :return: Yields either the next Record, or a StopIteration Exception if no more alignments
        """
        
        # Code to decompress a '*.gz' file should go here
        
        # Define empty vars
        header = None
        info = None
        sequence = None
        
        with open(filename) as flo:
            for line in flo:
                # Strip line ending characters
                line = line.rstrip()
                
                # Ignore blank lines
                if len(line) > 0:
                    # If this is the header line
                    m = regex.match(r'^>(\S*)\s*(.*)$', line)
                    if m:
                        if header:
                            yield self.create_record(header, info, sequence)
                        header = m.group(1)
                        info = m.group(2)
                        sequence = ''
                    else:
                        sequence += line
            if header:
                yield self.create_record(header, info, sequence)
    
    def create_record(self, header, info, sequence):
        '''
        Converts a FASTA record to a record object (tuple)
        '''
        # We ignore pairwise sequence attributes
        record = MSARecord(header, info, sequence)
        
        # Return the processed Record
        return record
    
    # TODO: function 'load_clustal' is unfinished and broken
    def load_clustal(self, filename):
        """
        Loads the next record of the CLUSTAL file.
        
        Returns:
         None if the file is complete
         Or a Record object
        """
        
        # Describe ALN/CLUSTAL format here:
        
        # Usually identified with the suffix '.aln')
        # Sequence names can be truncated to 10 or 30 characters
        # Residues may or may not be numbered
        # The order in which the sequences appear in the final alignment:
        #  'aligned'    Determined by the guide tree (default)
        #  'input'      Same order as the input sequences
        
        # Format Specification:
        #   1) The first line in the file must start with the words "CLUSTAL", "CLUSTAL W", or "CLUSTALW".
        #      Other information in the first line is ignored.
        #   2) One or more empty lines.
        #   3) One or more blocks of sequence data. Each block consists of:
        #       * One line for each sequence in the alignment. Each line consists of:
        #          1) the sequence name (up to 15 characters)
        #          2) white space
        #          3) up to 60 sequence symbols
        #          4) optional - white space followed by a cumulative count of residues for the sequences
        #       * A line showing the degree of conservation for the columns of the alignment in this block.
        #       * One or more empty lines.
        #   Some rules about representing sequences:
        #    * Case doesn't matter.
        #    * Sequence symbols should be from a valid alphabet.
        #    * Gaps are represented using hyphens ("-").
        #    * The characters used to represent the degree of conservation are
        #       '*'  -- all residues or nucleotides in that column are identical.
        #       ':'  -- conserved substitutions have been observed.
        #       '.'  -- semi-conserved substitutions have been observed.
        #       ' '  -- no match.
        
        # Code to decompress a '*.gz' file should go here
        
        # Define a non-existing record that will be replaced
        record = None
        
        flo = open(filename, 'r')
        
        # Create empty list to store records
        current_file = id(flo)
        self.records[current_file] = []
        
        # Get all records
        current_record_i = 0
        for i, line in enumerate(flo):
            line = line.rstrip()
            if (i == 0): # Skip this header line
                pass
            elif (len(line) > 0): # Only process the line if it contains data
                m = regex.match(r'(\S*)\s+(\S*)(?:\s+(\S*))?', line)
                if m:
                    try:
                        self.records[current_file][current_record_i] += m.group(1)
                    except IndexError:
                        self.records[current_file].append(m.group(1))
                    current_record_i += 1
                else: # This is the consensus line
                    pass ### Add something here??
            else: # This is an empty line
                current_record_i = 0
        
        
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

class MSARecord(object):
    __slots__ = [
        'header',
        'info',
        'sequence',
    ]
    
    def __init__(self, header, info, sequence):
        self.header = header
        self.info = info
        self.sequence = sequence
    
    def __repr__(self):
        """Return string representation of the Record"""
        return self.__class__.__name__ + '(' + ', '.join([
            'header={}'.format(self.header),
            'sequence={}'.format(self.sequence),
            ]) + ')'

# def test():
#     """Code to test the Bowtie2 Aligner subclass"""
#     if (len(sys.argv[1:]) != 2):
#         print("USAGE python3 mafft.py query.fasta folder", file=sys.stderr)
#         sys.exit(1)
#     query = sys.argv[1]
#     subject = None
#     folder = sys.argv[2]
#     
#     os.makedirs(folder, exist_ok=True)
#     logging.basicConfig(filename=os.path.join(folder, 'log.txt'), level=logging.INFO, format='%(message)s') # format='%(levelname)s %(asctime)s: %(message)s'
#     print("=== MAFFT ===")
#     aligner = Mafft()
#     #index_path = aligner.index(subject, 'myindex', folder, 4)
#     alignment_path = aligner.align(query, None, 'myoutput.fasta', folder, (os.cpu_count() or 1))
#     #print(index_path)
#     print(alignment_path)
#     
#     print('records = [')
#     for record in aligner.load(alignment_path):
#         print('  {}'.format(record))
#     print(']')
# 
# if (__name__ == '__main__'):
#     test()
