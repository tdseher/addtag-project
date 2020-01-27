#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/aligners/casoffinder.py

# List general Python imports
import os
from collections import OrderedDict
import logging
logger = logging.getLogger(__name__)

# Import non-standard package
import regex

# Import AddTag-specific packages
from .aligner import PairwiseAligner, Record
from ..utils import which
from ..motifs import OnTargetMotif, OffTargetMotif
from ..nucleotides import motif_conformation2
from ..cigarstrings import alignment2cigar, cigar2query_position, cigar2query_aligned_length, cigar2subject_aligned_length

class CasOFFinder(PairwiseAligner):
    logger = logger.getChild(__qualname__)
    
    def __init__(self):
        super().__init__(
            name="cas-offinder",
            authors=['Bae, Sangsu', 'Park, Jeongbin', 'Kim, Jin-Soo'],
            title='Cas-OFFinder: a fast and versatile algorithm that searches for potential off-target sites of Cas9 RNA-guided endonucleases',
            journal='Bioinformatics',
            issuing='30(10):1473-1475',
            year=2014,
            doi='https://doi.org/10.1093/bioinformatics/btu048',
            input='cas-offinder-in',
            output='cas-offinder-out',
            truncated=False,
            #classification='pairwise'
        )
        self.score_matrix ={}
        self.ev = {}
    
    def is_available(self):
        if which('cas-offinder'):
            return True
        else:
            return False
    
    def index(self, fasta, output_prefix, output_folder, threads, *args, **kwargs):
        """
        Cas-OFFinder does not require separate indexing of the reference genome.
        """
        return fasta
    
    def convert_query_file(self, input_filename, reference_filename, folder):
        '''
        Convert input FASTA to CAS-OFFINDER file format
        Because CAS-OFFINDER does not keep track of sequence labels,
        for simplicity, each FASTA sequence is given its own CAS-OFFINDER input file
        '''
        # Example input file:
        #   /var/chromosomes/human_hg19
        #   NNNNNNNNNNNNNNNNNNNNNRG
        #   GGCCGACCTGTCGCTGACGCNNN 5
        #   CGCCAGCGTCAGCGACAGGTNNN 5
        #   ACGGCGCCAGCGTCAGCGACNNN 5
        #   GTCGCTGACGCTGGCGCCGTNNN 5
        
        qname = os.path.splitext(os.path.basename(input_filename))[0]
        
        # Dict to hold open filehandles for new output files
#        flo_dict = {}
        
#        out_dict = {}
#        
#        # Consolidate Motif masks
#        for M in OnTargetMotif.motifs + OffTargetMotif.motifs:
#            for spacer_list, pam_list, side in M.parsed_list:
#                for pam in pam_list:
#                    for spacer in spacer_list:
#                        if (side == '>'):
#                            mask = spacer + pam
#                        else:
#                            mask = pam + spacer
#                        
#                        out_dict.setdefault(mask, []).append((M.motif_string, spacer, pam, side))
#        
#        # Write the first two lines for CASOFFINDER input files
#        for i, (mask, m_tuple) in enumerate(out_dict.items()):
#            flo = open(qname+'-'+str(i) + '.casoffinder', 'w+')
#            flo_dict[mask] = flo
#            print(reference_filename, flo)
#            print(mask, flo)
        
#        # Dict to hold the data that will be written to files
#        queries = {}
#        
#        # Determine the specific (spacer, pam, side) tuple that the query will be written to 
#        key = ('N'*len(spacer), pam, side)
#        
#        # Add the query to the list of all queries for the determined file
#        if (side == '>'):
#            queries.setdefault(key, []).append(spacer + 'N'*len(pam))
#        else:
#            queries.setdefault(key, []).append('N'*len(pam) + spacer)
        
        # Filenames produced
        outfiles = []
        
        mismatches = 5
        
        # Parse the FASTA file
        with open(input_filename) as in_flo:
            # Create empty 'header' and 'sequence' vars
            header = None
            sequence = None
            spacer = None
            pam = None
            side = None
            
            # Go through each line
            for line in in_flo:
                line = line.rstrip()
                if (len(line) > 0):
                    # If the current line is a header
                    if line.startswith('>'):
                        if header and sequence:
                            spacer, pam, side = self.motif_parse(motif, sequence)
                            outfiles.append(self.write_file(folder, qname, header, spacer, pam, side, reference_filename, mismatches))
                            
                        # Scrape the 'motif=...' field from the header
                        m = regex.search(r'>(\S*).*\bmotif=(\S*)', line)
                        if m:
                            header = m.group(1)
                            motif = m.group(2)
                        else:
                            header = regex.split(r'\s', line[1:], 1)[0]
                            motif = None
                        sequence = ''
                    
                    # Otherwise, if the current line has sequence data
                    else:
                        # This sequence SHOULD have the PAM
                        sequence += line
            
            # Process the remaining sequence
            if header and sequence:
                spacer, pam, side = self.motif_parse(motif, sequence)
                outfiles.append(self.write_file(folder, qname, header, spacer, pam, side, reference_filename, mismatches))
        
        # Return list of all CAS-OFFINDER inputs generated
        return outfiles
    
#    def process_motif_sequence_pair(self, motif, sequence):
#        # Get the file associated with this 'motif'
#        for m_mask, m_tuple in out_dict.items():
#            for motif_string, spacer_list, pam_list, m_side in m_tuple:
#                if (motif == motif_string):
#                    mask = m_mask
#                    side = m_side
#                    break
#        if (side == '>'):
#            print(sequence + 'N'*len(pam), flo_dict[mask])
#        else:
#            print('N'*len(pam)+sequence, flo_dict[mask])
    
    def motif_parse(self, motif, sequence):
        # Find the Motif object associated with the 'motif=...' text
        for M in OnTargetMotif.motifs + OffTargetMotif.motifs:
            if (motif == M.motif_string):
                spacer_list, pam_list, side = M.parsed_list
                m = motif_conformation2(sequence, side, M.compiled_regex)
                if m:
                    spacer = m[0]
                    for p in pam_list:
                        p_regex = M.build_regex_pattern(p)
                        p_match = regex.match(p_regex, m[1])
                        if p_match:
                            pam = p
                            break
                    #pam = m[1]
                    return spacer, pam, side
        
        # If no Motif object exists, then create one
        M = OnTargetMotif(motif)
        spacer_list, pam_list, side = M.parsed_list
        m = motif_conformation2(sequence, side, M.compiled_regex)
        if m:
            spacer = m[0]
            pam = m[1]
            return spacer, pam, side
    
    def write_file(self, folder, qname, header, spacer, pam, side, subject, mismatches):
        filename = os.path.join(folder, qname+'-'+header+'.cas-offinder-in')
        with open(filename, 'w') as out_flo:
            # Add path to reference FASTA
            print(subject, file=out_flo)
            
            # Print the mask and sequence, mismatches
            if (side == '>'):
                print('{}{}'.format('N'*len(spacer), pam), file=out_flo) # mask
                print('{}{} {}'.format(spacer, 'N'*len(pam), mismatches), file=out_flo) # sequence, mismatches
            else:
                print('{}{}'.format(pam, 'N'*len(spacer)), file=out_flo) # mask
                print('{}{} {}'.format('N'*len(pam), spacer, mismatches), file=out_flo) # sequence, mismatches
        return filename
    
    def align(self, query, subject, output_prefix, output_folder, threads, *args, **kwargs):
        """
        """
        # Folder to store all Cas-OFFinder input files
        #casoffinder_folder = os.path.join(output_folder, 'casoffinder')
        #os.makedirs(casoffinder_folder, exist_ok=True)
        
        # Create a folder to store query files
        #qname = os.path.splitext(os.path.basename(query))[0]
        #folder = os.path.join(os.path.dirname(query), qname + '.cas-offinder')
        folder = os.path.join(output_folder, os.path.splitext(os.path.basename(query))[0] + '.cas-offinder')
        os.makedirs(folder, exist_ok=True)
        
        # qname = os.path.splitext(os.path.basename(query))[0]
        # out_file_paths = {}
        # 
        # # Create a single input for each motif to search
        # for parsed_motif in kwargs['parsed_motifs']:
        #     spacers, pams, side = parsed_motif
        #     for i, s in enumerate(spacers):
        #         for j, p in enumerate(pams):
        #             if (side == '>'):
        #                 motif = s+p
        #             else:
        #                 motif = p+s
        #             
        #             # Build special input file for Cas-OFFinder
        #             # /path/to/reference/fasta/folder
        #             # NNNNNNNNNNNNNNNNNNNNNRG
        #             # GGTGTTGAGAGTGCCATCGANNN 5
        #             # GTCGTTTTGATTTGGATGTANNN 5
        #             # TCAAAATCGTGTTGCGAAATNNN 5
        #             input_filename_path = os.path.join(casoffinder_folder, qname+'.casoffinder-input-'+str(i)+str(j))
        #             output_filename_path = os.path.join(casoffinder_folder, qname+'.casoffinder-output-'+str(i)+str(j))
        #             out_file_paths[(s,p)] = output_filename_path
        #             with open(input_filename_path, 'w') as flo:
        #                 # Add path to reference FASTA
        #                 print(os.path.dirname(subject), file=flo)
        #             
        #                 # Expand abbreviated motif into the full IUPAC string
        #                 print(motif, file=flo)
        #                 
        #                 # Print all the sequences to search the genome for, as well as the maximum difference
        #                 for seq in sequences:
        #                     print(seq + ' ' + str(5), file=flo)
        
        
        inputs = self.convert_query_file(query, subject, folder)
        outputs = []
        
        for infile in inputs:
            outfile = os.path.splitext(infile)[0] + '.cas-offinder-out'
            
            # Options to run Cas-OFFinder with
            options = OrderedDict([
                (infile, None),
                ('G', None), # {C|G|A}[device_id(s)] (C: using CPUs, G: using GPUs, A: using accelerators) omit number to use all like devices
                (outfile, None),
            ])
            
            outputs.append(self.process('cas-offinder', outfile, options))
        
        return folder

    def load(self, filename, *args, **kwargs):
        # Output file format (tab delimited)
        #   Col  Field   Description
        #     0  QSEQ    given query sequence
        #     1  SNAME   FASTA sequence title (if you downloaded it from UCSC or Ensembl, it is usually a chromosome name)
        #     2  SPOS    0-based position of the potential off-target site
        #     3  SSEQ    actual sequence located at the position (mismatched bases noted in lowercase letters)
        #     4  STRAND  indicates forward strand(+) or reverse strand(-) of the found sequence
        #     5  ERRORS  the number of the mismatched bases ('N' in PAM sequence are not counted as mismatched bases)
        
        # An example of output file (tab delimited):
        #   GGCCGACCTGTCGCTGACGCNNN chr8    49679 GGgCatCCTGTCGCaGACaCAGG + 5
        #   GGCCGACCTGTCGCTGACGCNNN chr8   517739 GcCCtgCaTGTgGCTGACGCAGG + 5
        #   GGCCGACCTGTCGCTGACGCNNN chr8   599935 tGCCGtCtTcTCcCTGACGCCAG - 5
        #   GGCCGACCTGTCGCTGACGCNNN chr8  5308348 GGCaGgCCTGgCttTGACGCAGG - 5
        #   GGCCGACCTGTCGCTGACGCNNN chr8  9525579 GGCCcAgCTGTtGCTGAtGaAAG + 5
        #   GGCCGACCTGTCGCTGACGCNNN chr8 12657177 GGCCcACCTGTgGCTGcCcaTAG - 5
        
        # Get all '*.cas-offinder-out' files within 'filename'
        file_list = []
        with os.scandir(filename) as it:
            for entry in it:
                if entry.name.endswith('.cas-offinder-out') and entry.is_file():
                    #file_list.append(entry.name) # Just the basename
                    file_list.append(entry.path) # The full path
        
        for outfile in file_list:
            with open(outfile) as flo:
                for line in flo:
                    # Remove newline character from end of line
                    line = line.rstrip()
                    if (len(line) > 0):
                        # Split CAS-OFFINDER record str into a list
                        sline = line.rstrip().split("\t")
                        
                        # Require the line to have a record (not a blank line)
                        # Require the record to align to a contig
                        yield self.create_record(sline, outfile)
        
    def create_record(self, sline, filename):
        '''
        Converts a split line from a SAM file into a Record object
        :param sline: line.rstrip().split('\t')
        :return: a Record object
        '''
        
        # Look up query name
        qname = os.path.splitext(os.path.basename(filename))[0].split('query-')[1]
        #qname = Targets.index[sline[0]]
        
        # Truncate sname
        sname = regex.split(r'\s', sline[1], 1)[0]
        
        qseq = sline[0]
        sseq = sline[3]
        
        # This CIGAR does not necessarily conform the the MOTIF
        cigar = alignment2cigar(qseq, sseq, specific=True, abbreviated=False)
        
        qlen = cigar2query_aligned_length(cigar)
        slen = cigar2subject_aligned_length(cigar)
        
        sstart = int(sline[2])
        qpos = cigar2query_position(cigar)
        spos = (sstart, sstart+slen)
        
        # Calculate the FLAGS
        # TODO: Currently does not add the 0x100 (non-primary alignment) flag
        flags = 0
        if (sline[4] == '-'):
            flags |= 0x10
        
        # We ignore these for now
        score = None # TODO: Does not calculate a score
        evalue = None # TODO: Does not calculate E-value
        length = None
        
        # Build Record
        record = Record(
            qname, sname, # query_name, subject_name,
            qseq, sseq, # query_sequence, subject_sequence,
            qpos, spos, # query_position, subject_position,
            qlen, slen, # query_length, subject_length,
            flags, cigar, score, evalue, length # flags, cigar, score, evalue, length
        )
        
        # Return the processed Record
        return record
