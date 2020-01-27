#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/aligners/bwa.py

# List general Python imports
import os
from collections import OrderedDict
import logging
logger = logging.getLogger(__name__)

# Import non-standard packages
import regex

# import AddTag-specific packages
from .aligner import PairwiseAligner, Record
from ..cigarstrings import cigar2query_position, cigar2query_aligned_length, cigar2subject_aligned_length, cigar2score, specificate_cigar, collapse_cigar
from ..utils import old_load_fasta_file, which
from ..evalues import EstimateVariables, load_scores
from ..nucleotides import rc

class Bwa(PairwiseAligner):
    logger = logger.getChild(__qualname__)

    def __init__(self):
        super().__init__(
            name="bwa",
            authors=['Li, Heng', 'Durbin, Richard'],
            title='Fast and accurate short read alignment with Burrows–Wheeler transform',
            journal='Bioinformatics',
            issuing='25(14):1754-1760',
            year=2009,
            doi='https://doi.org/10.1093/bioinformatics/btp324',
            input='fasta', # 'fastq',
            output='sam',
            truncated=False,
            #classification='pairwise'
        )
        self.score_matrix = {}
        self.ev = {}
        self.references = {}

    def is_available(self):
        if which('bwa'):
            return True
        else:
            return False

    def index(self, fasta, output_prefix, output_folder, threads, *args, **kwargs):
        """
        Call 'bwa index' on non-compressed FASTA file.

        The indexed FASTA will be stored in 'folder'.
        The 'fasta' basename will be used as both 'folder' and index name.

        Additionally, compute the Karlin-Altschul statistics parameters
        necessary for E-value calculations
        """
        # Karlin-Altschul calculations
        score_matrix = load_scores(os.path.join(os.path.dirname(__file__), 'bwa_scores.txt'))
        scontigs = old_load_fasta_file(fasta)
        ev = EstimateVariables(score_matrix, scontigs)

        # Find the FASTA basename
        #name = os.path.splitext(os.path.basename(output_filename))[0]
        name = output_prefix

        # Make the directory if it does not yet exist
        os.makedirs(os.path.join(output_folder, name), exist_ok=True)

        index_file = os.path.join(output_folder, name, name)

        options = OrderedDict([
            ('index', None), # Command
            ('-p', index_file), # Prefix of the output database [same as db filename]
            ('-a', 'is'), # Algorithm for constructing BWT index
            (fasta, None),
            #(index_file, None),
        ])
        outpath = self.process('bwa', index_file, options)
        
        self.ev[outpath] = ev
        self.score_matrix[outpath] = score_matrix
        self.references[outpath] = scontigs
        
        return outpath
    
    def align(self, query, subject, output_prefix, output_folder, threads, *args, **kwargs):
        '''
        Align the query file to the subject file.
        :param query: The path of the FASTQ to align
        :param subject: The path of the BWA index prefix
        :param output_filename: The path of the SAM file to generate
        :param output_folder: The directory to store the generated SAM file
        :param threads: Number of processors to use
        :return: The path of the SAM file generated
        '''
        
        #output_basename = os.path.splitext(os.path.basename(output_filename))[0]
        output_basepath = os.path.join(output_folder, output_prefix)
        output_sai_path = output_basepath + '.sai'
        output_sam_path = output_basepath + '.sam'
        
        #output_filename_path = os.path.join(output_folder, output_filename)
        options = OrderedDict([
            ('aln', None),           # Command
            ('-n', 5),               # docs: [-n maxDiff] Maximum edit distance if the value is INT, or the fraction
                                     #       of missing alignments given 2% uniform base error rate if FLOAT.
                                     #       In the latter case, the maximum edit distance is automatically
                                     #       chosen for different read lengths. [0.04]
                                     # bin: max #diff (int) or missing prob under 0.02 err rate (float) [-1.00]
            ('-o', 1),               # Maximum number or fraction of gap opens [1]
            ('-e', -1),              # Maximum number of gap extensions, -1 for k-difference mode
                                     # (disallowing long gaps) [-1]
            ('-i', 0),               # Disallow an indel within INT bp towards the ends [5]
            ('-d', 10),              # docs: Disallow a long deletion within INT bp towards the 3’-end [16]
                                     # bin: maximum occurrences for extending a long deletion [10]
            ('-l', 0),               # docs: Take the first INT subsequence as seed. If INT is larger than the query
                                     #       sequence, seeding will be disabled. For long reads, this option is
                                     #       typically ranged from 25 to 35 for '-k 2'. [inf]
                                     # bin: seed length [0]
            ('-k', 5),               # docs: Maximum edit distance in the seed [2]
                                     # bin: maximum differences in the seed [5]
            ('-t', threads),         # Number of threads (multi-threading mode) [1]
            ('-M', 3),               # [-M misMsc] Mismatch penalty. BWA will not search for suboptimal hits
                                     # with a score lower than (bestScore-misMsc). [3]
            ('-O', 11),              # Gap open penalty [11]
            ('-E', 4),               # Gap extension penalty [4]
            ('-N', None),            # docs: Disable iterative search. All hits with no more than maxDiff
                                     #       differences will be found. This mode is much slower than the default.
                                     # bin: non-iterative mode: search for all n-difference hits (slooow)
            ('-f', output_sai_path), # SAI output file to write output to instead of STDOUT
            (subject, None),         # Prefix path of BWA index
            (query, None),           # FASTQ with unaligned reads
        ])
        
        outpath1 = self.process('bwa', output_sai_path, options)
        
        options = OrderedDict([
            ('samse', None),            # Command
            ('-n', 1000000),            # Maximum number of alignments to output in the XA tag for reads paired
                                        # properly. If a read has more than INT hits, the XA tag will not be
                                        # written. [3]
            ('-f', output_sam_path),    # Write SAM to this file instead of STDOUT
            (subject, None),            # Prefix path of BWA index
            (output_sai_path, None),    # SAI file
            (query, None),              # FASTQ with unaligned reads
        ])
        
        outpath2 = self.process('bwa', output_sam_path, options, append=True)
        
        # Add the output file as an additional key pointing to this database
        self.ev[outpath2] = self.ev[subject]
        self.score_matrix[outpath2] = self.score_matrix[subject]
        self.references[outpath2] = self.references[subject]
        
        return outpath2
    
    def load(self, filename, *args, **kwargs):
        '''
        Yields records from the input SAM file, one at a time.
        BWA SAM files can have multiple Records per line. Thus, each line is parsed to retrieve the list of Records.
        :param filename: Path of the file to create Records from
        :param args:
        :param kwargs: 
        :return: Yields either the next Record, or a StopIteration Exception if no more alignments
        '''
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
                        
                        # parse the line
                        plines = self.parse_line(line)
                        for pline in plines:
                            yield self.create_record(pline, filename)
    
    def create_record(self, sline, filename):
        '''
        Converts a split line from a SAM file into a Record object
        :param sline: line.rstrip().split('\t')
        :return: a Record object
        '''
        # Find evalue of the alignment
        ev = self.ev[filename]
        score_matrix = self.score_matrix[filename]
        
        # TODO: We store an entire copy of the genome in memory, just in order to get better
        #       CIGAR strings so we can get accurate 'evalue's. This a waste of computing resources!
        qseq = sline[9]
        sstart, send = (int(sline[3])-1, int(sline[3])-1+cigar2subject_aligned_length(sline[5]))
        sseq = self.references[filename][sline[2]][sstart:send]
        cigar = collapse_cigar(specificate_cigar(sline[5], qseq, sseq))
        ##### End #####
        
        cigar_score = cigar2score(cigar, score_matrix)
        evalue = ev.calculate_evalue(cigar_score, sline[9])
        
        # Build Record
        record = Record(
            sline[0], sline[2], # query_name, subject_name,
            sline[9], None, # query_sequence, subject_sequence,
            cigar2query_position(sline[5]), (int(sline[3])-1, int(sline[3])-1+cigar2subject_aligned_length(sline[5])), # query_position, subject_position,
            cigar2query_aligned_length(sline[5]), cigar2subject_aligned_length(sline[5]), # query_length, subject_length,
            int(sline[1]), cigar, float(sline[4]), evalue, None # flags, cigar, score, evalue, length
        )
        
        # Return the processed Record
        return record
    
    def parse_line(self, line, outfmt=list):
        '''
        Decode the 'XA' tag of BWA SAM files
        Returns same output as 'xa2multi.pl' from BWA
        :param line: line of text to parse
        :return: list of lines
        '''
        
        # List to store output lines
        output = []
        
        sline = line.rstrip().split('\t')
        
        if issubclass(outfmt, list):
            output.append(sline)
        else:
            output.append(line)
        
        #xa = None
        #xai = None
        #for i, field in enumerate(sline[11:]):
        #    xa = regex.search(r'XA:Z:(\S+)', field)
        #    if xa:
        #        xai = 11+i
        #        break
        
        xa = regex.search(r'XA:Z:(\S+)', line)
        
        if xa:
            #output.append(sline[:xai] + sline[xai+1:])
            
            for m in regex.finditer(r'([^,;]+),([-+]\d+),([^,]+),(\d+);', xa.group(1)):
                contig, pos, cigar, errors = m.groups()
                pos = int(pos)
                
                # TODO: Add Insert size calculation (sline[8] = TLEN/ISIZE)
                #if (sline[6] == contig):
                #    my_chr = '='
                #else:
                #    my_chr = sline[6]
                insert = 0
                
                seq = sline[9]
                qual = sline[10]
                
                # If alternative alignment has other orientation than primary,
                # then print the reverse (complement) of sequence and phred string
                
                flag = 0x100 # Mark as non-primary/supplementary alignment
                if ((pos < 0) ^ ((int(sline[1]) & 0x10)>0)): # If either (pos < 0) OR (read maps in '-' orientation), but not both, and not neither
                    seq = rc(seq)
                    qual = qual[::-1]
                flag |= (int(sline[1]) & 0b111011101001) # Inherit all flags from primary alignment except 0x10 and 0x100
                if (pos < 0):
                    flag |= 0x10 # Mark as '-' strand
                
                if issubclass(outfmt, list):
                    output.append([sline[0], flag, contig, abs(pos), 0, cigar, sline[6], sline[7], insert, seq, qual, 'NM:i:'+errors])
                else:
                    output.append('\t'.join(map(str, [sline[0], flag, contig, abs(pos), 0, cigar, sline[6], sline[7], insert, seq, qual, 'NM:i:'+errors])))
        
        # Return the list of lines/slines
        return output
