#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/subroutines/_subroutine_generate_primers.py

# Import standard packages
import sys
import os
import logging
import copy
import random
import math
import time
import multiprocessing
import queue
#from collections import namedtuple

# Import non-standard packages
import regex

# Import included AddTag-specific modules
from . import subroutine
from ..feature import Feature
from .. import utils
from .. import nucleotides
from .. import aligners
from .. import thermodynamics
from ..thermodynamics.oligo import Oligo, Primer, PrimerPair, PrimerDesign
from .. import cache
from ..algorithms import detect_overridden_methods

logger = logging.getLogger(__name__)

# TODO: If the absolute number of potential Primers in a region (UPSTREAM_F, DOWNSTREAM_R, FEATURE_F, FEATURE_R)
#       is small (like <1000), then just do ALL the Primer and PrimerPair calculatons FOR THAT REGION
#       Thus, the other regions don't have to be repeatedly calculated over and over.
#       This probably should come AFTER making it so each 'PrimerPair' has its own 'Cutoff' level.
# TODO: When Primer and PrimerPairs are queued for processing, just set their priority level based on their
#       thermodynamic properties OR by their quasi-likelihood?

class GeneratePrimersParser(subroutine.Subroutine):
    logger = logger.getChild(__qualname__)
    
    def __init__(self, subparsers):
        self.subparsers = subparsers
        
        self.name = 'generate_primers'
        self.description = (
            "description:" "\n"
            "  Design primers for confirming whether each step of genome engineering is" " \n"
            "  successful or not. This does not design multiplexable primers." "\n"
        )
        self.help = "Design primers for confirming whether each step of genome engineering is successful or not."
        self.epilog = (
            "example:" "\n"
            "  You can generate cPCR oligos by running:" "\n"
            "   $ python3 {__program__} {__subroutine__}" "\n"
            "     --fasta genome.fasta" "\n"
            "     --dDNAs ko-dDNA.fasta ki-dDNA.fasta" "\n"
            "     --folder output" "\n"
            "     > {__subroutine__}.out 2> {__subroutine__}.err"
        ).format(**dict(list(subroutine.__dict__.items()) + list({"__subroutine__": self.name}.items()))) # key:value pairs in the latter dict will replace any instances in the prior dict
        
        self.define_parser()
        self.define_arguments()
    
    def define_arguments(self):
        # Add required arguments
        required_group = self.parser.add_argument_group('required arguments')
        required_group.add_argument("--fasta", required=True, nargs="+", metavar="*.fasta", type=str,
            help="FASTA files with contigs to find unique gRNA sites. Input FASTA \
            must not be compressed. User can decide whether ambiguous bases can \
            be chosen for gRNA sites. All FASTA sequences should have unique \
            primary headers (everything between the '>' symbol and the first \
            whitespace should be unique). You should include FASTA of the genome \
            and any plasmids.")
        
        required_group.add_argument("--dDNAs", required=True, nargs="+", metavar="*.fasta", type=str,
            help="A dDNA FASTA file for each subsequent CRISPR/Cas transformation. \
            Typically, the first is KO, and the second is KI. However, any number \
            of serial genome engineering experiments can be specified.")
        
        required_group.add_argument("--folder", required=True, metavar="FOLDER",
            type=str, help="Path of folder to store generated files.")
        
        # Add optional arguments
        self.parser.add_argument("--processors", metavar="N", type=int, default=(os.cpu_count() or 1),
            help="Number of processors to use when performing pairwise sequence alignments.")
        
        aligner_choices = [x.name for x in aligners.pw_aligners]
        self.parser.add_argument("--aligner", type=str, choices=aligner_choices, default='blastn',
            help="Program to calculate pairwise alignments. Please note that the 'addtag' internal aligner is very slow.")
        
        #self.parser.add_argument("--number_pcr_conditions", metavar="N", type=int, default=None,
        #    help="Number of PCR conditions to develop primers for. All amplicons \
        #    within each condition will have similar size (nt) and melting temperatures. \
        #    If unspecified, will default to the number of target features.")
        
        self.parser.add_argument("--case", type=str, default="ignore",
            choices=["ignore", "upper-only", "lower-only", "mixed-lower", "mixed-upper", "mixed-only"],
            help="Restrict generation of primers based on case of nucleotides in input FASTA: \
            ignore - keep all potential primers; \
            upper-only - discard any potential primer with a lower-case character in the genome; \
            lower-only - discard all potential primers with upper-case characters; \
            mixed-lower - discard primers that are all upper-case; \
            mixed-upper - discard primers that are all lower-case; \
            mixed-only - only use primers that have both lower- and upper-case characters.")
        
        self.parser.add_argument("--max_number_designs_reported", metavar="N", type=int, default=5,
            help="The maximum number of final cPCR designs to report.")
        
        self.parser.add_argument("--primer_scan_limit", metavar="N", type=int, default=None, #2*60, # Should be None by default, thus a limit is only enforced if this option is specified.
            help="Number of seconds to limit each primer scan.")
        
        self.parser.add_argument("--primer_pair_limit", metavar="N", type=int, default=None, #5*60, # Should be None by default, thus a limit is only enforced if this option is specified.
            help="Amount of time (in seconds) to limit primer pairings.")
        
        self.parser.add_argument("--mandatory_primers", nargs="+", metavar="*.fasta", type=str, default=[],
            help="A FASTA file for each round containing primer sequences that \
            must be used for that round. Usually, these correspond to flanktags.")
        
        oligo_choices = [x.name for x in thermodynamics.oligos]
        self.parser.add_argument("--oligo", type=str, choices=oligo_choices, default='Primer3',
            help="Program to perform thermodynamic calculations.")
        
        self.parser.add_argument("--max_evalue", metavar="N", type=float, default=0.001,
            help="The maximum accepted alignment E-value to consider a recombination.")
        
        self.parser.add_argument("--min_length", metavar="N", type=int, default=35,
            help="The minimum accepted alignment length to consider a recombination.")
        
        self.parser.add_argument("--skip_round", metavar="N", nargs="+", type=int, default=[],
            help="Skip primer calculations for these rounds.")
        
        self.parser.add_argument("--cycle_start", metavar="N", type=int, default=0,
            help="Stringency to consider a primer.")
        
        a = self.parser.add_argument("--cycle_stop", metavar="N", type=int, default=None,
            action=subroutine.ValidateCycleStop,
            help="Stringency to consider a primer. Defaults to the same value as '--cycle_start'.")
        a.process_action_on_default = True
        
        self.parser.add_argument("--subset_size", metavar="N", type=int, default=1000,
            help="Artificially limit the number of primer pairs that are calculated. \
            This represents the max each primer list could be. The maximum number \
            of paiwise comparisons is the square of this number.")
        
        self.parser.add_argument("--cache", action="store_true", default=False,
            help="If you will be re-running designs on these same input files, \
            Intermediate results can be stored and retrieved to speed up \
            total calculations.")
        
        # Temporary: this expects 2 for each dDNA: a before and an after
#        self.parser.add_argument("--internal_primers_required", metavar="y/n",
#            nargs="+", type=str, default=None, action=subroutine.ValidateInternalPrimersRequired,
#            help="For each genome, starting with input, and each subsequent dDNA, \
#            specify whether internal primers are required (oF/oR). \
#            y - yes, these internal primers are required; \
#            n - no, these internal primers are optional. \
#            (For example, if you have 2 rounds of genome engineering, \
#            and only the wild type (input) genome and the final genome \
#            require internal primers, then you would use these command-line \
#            options: '--fasta genome.fasta --dDNAs ko.fasta ki.fasta \
#            --internal_primers_required y n y')")
        
        self.parser.add_argument("--i_primers_required", metavar="y/n",
            nargs="+", type=str, default=None, action=subroutine.ValidatePrimersRequired,
            help="For each genome, starting with input, and each subsequent dDNA, \
            specify whether internal primers are required (iF/iR). \
            y - yes, these internal primers are required; \
            n - no, these internal primers are optional. \
            (For example, if you have 2 rounds of genome engineering, \
            and only the wild type (input) genome and the final genome \
            require internal primers, then you would use these command-line \
            options: '--fasta genome.fasta --dDNAs ko.fasta ki.fasta \
            --i_primers_required y n y')")
        
        self.parser.add_argument("--o_primers_required", metavar="y/n",
            nargs="+", type=str, default=None, action=subroutine.ValidatePrimersRequired,
            help="For each genome, starting with input, and each subsequent dDNA, \
            specify whether outer primers are required (oF/oR). \
            y - yes, these outer primers are required; \
            n - no, these outer primers are optional. \
            (For example, if you have 2 rounds of genome engineering, \
            and only the wild type (input) genome and the final genome \
            require outer primers, then you would use these command-line \
            options: '--fasta genome.fasta --dDNAs ko.fasta ki.fasta \
            --o_primers_required y n y')")
        
        self.parser.add_argument("--skip_partial_sets", action="store_true", default=False,
            help="Do not perform simulated annealing on primer sets that are missing required PrimerPairs.")
        
        # Default should be multi-allelic (NOT allele-agnostic).
        # If the users want allele-specific, then they should use this
        # command-line option.
        #self.parser.add_argument("--allele-specific", action="store_true", default=False,
        #    help="Report only allele-specific primer designs. \
        #    Either primer amplicons will be diagnostically-different sizes, \
        #    or primer sequences themselves will be different.")
        
        #   'exclusive' (single/specific),
                    #   'all' (multi),
                    #   'any' (agnostic)
        self.parser.add_argument("--specificity", choices=['exclusive', 'all', 'any'], default='all',
            help="Report only *-specific primer designs. Either primer \
            amplicons will be diagnostially-different sizes, or primer \
            sequences themselves will be different. The choices are as \
            follows: 'exclusive' (allele-specific), 'all' (multi-allelic), \
            'any' (allele-agnostic).")
        
        # Get TEMP directory
        if sys.platform.startswith('win'):
            from tempfile import gettempdir
            temp_default = gettempdir()
        else: # 'linux' or 'darwin'
            temp_default = '/dev/shm'
        self.parser.add_argument("--temp_folder", type=str, default=temp_default,
            help="Directory to store temporary files. RAMdisk recommended.")
        
        # Nucleotide matching stuff
        #  - number errors (for fuzzy regex)
        
        # PCR conditions:
        #self.parser.add_argument("--primer_lengths", nargs=2, metavar=("MIN", "MAX"), type=int, default=[19,28],
        #    help="Range of lengths acceptable for primer sequences, inclusive.")
        #self.parser.add_argument("--primer_last5gc_count", nargs=2 metavar=("MIN", "MAX"), type=int, default=[1,3],
        #    help="Range of acceptable number of G and C residues in the last 5 nt of oligo, inclusive.")
        #self.parser.add_argument("--primer_gc_clamp_length", nargs=2 metavar=("MIN", "MAX"), type=int, default=[0,2],
        #    help="Range of acceptable consecutive G and C residues in 3' end of oligo, inclusive.")
        #self.parser.add_argument("--primer_gc", nargs=2, metavar=("MIN", "MAX"), type=float, default=[0.25,0.75],
        #    help="Minimum and maximum acceptable fraction of C and C residues within the primer.")
        #self.parser.add_argument("--primer_max_run_length", metavar="N", type=int, default=4,
        #    help="Maximum number of consecutive identical residues.")
        
        # PCR conditions:
        #  - primer_size (min, max)
        #  - monovalent_cation_concentration
        #  - divalent_cation_concentration
        #  - amplicon_size (min, max)
        #  - template_concentration
        #  - dNTP_concentration
        #  - temperature (25 usually)
        #  - max_tm_difference
        #  - melting_temperature (min, max)
        #  - max_3prime_homology_length
        #  - min_delta_g
        #  - gc (min, max)
        #  - gc_clamp_length (min, max)
        #  - max_run_length
    
    def compute(self, args):
        print("# Confirmation primer design (cPCR).")
        
        # Create the 'cache' subdirectory if it doesn't exist already
        os.makedirs(os.path.join(args.folder, 'cache'), exist_ok=True)
        
        # Set intended amplicon size range
        #amplicon_size = (400, 700) # definition moved to make_primer_set() function
        
        # Hard-coded melting temp range (should be adjustable via command line in the future)
        #tm_range = (53, 57) # definition moved to make_primer_set() function
        
        # Require primer pairs to be within 2 degrees Celcius from each other
        #tm_max_difference = 3.0 # definition moved to make_primer_set() function
            
        # Order in which primer lengths should be interrogated, starting from the left
        # (18 doesn't work reliably with UNAFold)
        #primer_sizes = [20, 21, 19, 22, 23, 18, 24, 25, 26]
        
        #args.number_pcr_conditions
        
        # Define variables to hold slines for output
        output0 = []
        output1 = []
        output2 = []
        output3 = []
        output4 = []
        
        # First, create the engineered genomes:
        #   FILENAME          DESCRIPTION
        #   genome-r0.fasta   wild type genomic DNA
        #   genome-r1.fasta   genomic DNA with dDNA1 incorporated
        #   genome-r2.fasta   genomic DNA with dDNA1 incorporated, then dDNA2 incorporated
        
        # For each round, create certain data structures or files
        genome_contigs_list = []
        dDNA_contigs_list = []
        dDNA_fasta_file_list = []
        dDNA_alignment_file_list = []
        genome_fasta_file_list = []
        
        #Datum = namedtuple('Datum', ['dDNA_r', 'dDNA_contig', 'genome_r', 'genome_contig', 'ush_start', 'ush_end', 'dsh_start', 'dsh_end', 'ins_start', 'ins_end'])
        
        # Create 'r0'
        self.logger.info('Working on round r{}'.format(0))
        # Parse the input FASTA
        genome_contigs_list.append(utils.load_multiple_fasta_files(args.fasta))
        dDNA_contigs_list.append(None) # The 'r0' round of genome engineering has no dDNA
        dDNA_fasta_file_list.append(None)
        dDNA_alignment_file_list.append(None)
        
        # Write the input genome FASTA to 'genome-r0.fasta'
        # Merge input FASTA files into a single one
        genome_fasta_file_list.append(utils.write_merged_fasta(genome_contigs_list[0], os.path.join(args.folder, 'genome-r0.fasta')))
        
        # We need the groups of contig names that correspond to the loci in each round
        # For example:
        #   contig_groups = [ ['chr1', 'chr1-r1[gene1A,gene2]', 'chr1-r1[gene1A,gene2]-r2[gene1B]'] ]
        # Populate contig_groups with the initial, 'r0' contig names
        contig_groups = [[c] for c in genome_contigs_list[0]]
        
        # We link the dDNA contigs together
        # For example:
        #   dDNA_groups = [ ['gene1A', 'gene1B'], ['gene2'] ]
        dDNA_groups = []
        
        # We link the actual genome modifications together
        # For example:
        #   mod_groups = [
        #     [('chr1', 1022, 2000, 2300, 4500), ('chr1-r1[gene1A,gene2]', 1322, 2300, 2600, 4800), ('chr1-r1[gene1A,gene2]-r2[gene1B]', 1055, 2000, 2301, 4000)],
        #     [('chr1', 200, 220, 240, 300), ('chr1-r1[gene1A,gene2]', 200, 240, 240, 280)]
        #  ]
        mod_groups = []
        
        # group_links = {
        #     # (dDNA round, dDNA contig name): gDNA r0                                          gDNA r1
        #     (1, 'gene1A'):                    ((0, 'chr1', 100, 200, 300, 400,),              (1, 'chr1-r1[gene1A,gene2]', 200,300,400,500)),
        #     (2, 'gene1B'):                    ((1, 'chr1-r1[gene1A,gene2]', 200,300,400,500), (2, 'chr1-r1[gene1A,gene2]-r2[gene1B]', 100,200,300,400))
        # }
        # Logic: Because the first tuple of 'gene 1B' overlaps with the second tuple of 'gene1A', these are linked
        group_links = [] # {}
        
        # Record calculated dDNA homology regions
        #dDNA_homology_seqs = {} # key=(dDNA_r, dDNA_contig, genome_contig), value=(s_ush_seq, s_dsh_seq)
        
        # For each dDNA, do some operations
        for j, dDNA_filename in enumerate(args.dDNAs):
            # Variable to count the round of genome engineering
            r = j+1
            
            # Let the user know which round is being processed
            self.logger.info('Working on round r{}'.format(r))
            
            # Load input dDNA FASTA into the list of dicts
            dDNA_contigs_list.append(utils.old_load_fasta_file(dDNA_filename))
            
            # Build an index of the previous round's genome
            #genome_index_file = args.selected_aligner.index(genome_fasta_file_list[r-1], os.path.basename(genome_fasta_file_list[r-1]), args.folder, args.processors)
            genome_index_file = args.selected_aligner.index(genome_fasta_file_list[r-1], self.pathstrip(genome_fasta_file_list[r-1]), args.folder, args.processors)
            
            # Write current dDNA to FASTA
            dDNA_fasta_file_list.append(utils.write_merged_fasta(dDNA_contigs_list[r], os.path.join(args.folder, 'dDNA-r'+str(r)+'.fasta')))
            
            # We align dDNA against the genome+dDNA(prev)
            dDNA_alignment_file_list.append(args.selected_aligner.align(dDNA_fasta_file_list[r], genome_index_file, 'dDNA-r'+str(r)+'-alignment', args.folder, args.processors))
            
            # Parse the alignment to make an intermediary record data structure
            records = self.filter_alignment_records_for_cPCR(args, dDNA_alignment_file_list[r], dDNA_contigs_list[r])
            
            # Iterate through all detected paired cross-over events (i.e. only simulate the "proper" cross-overs for significant alignment pairs)
            # Starting with the right-most position in each sequence, make changes to the genome contig sequence
            # We assume no record pairs overlap
            
            # Keep track of how contig names will change
            name_changes = {}
            
            # Seed the results of this round's engineering with the contigs from the previous round
            # (Will still need to change their names)
            genome_contigs_list.append(copy.deepcopy(genome_contigs_list[r-1]))
            
            # Sort the 'qs_pair' tuples from greatest to least based on position in each chromosome
            record_order = sorted(records, key=lambda x: records[x][0].subject_position, reverse=True)
            
            # Iterate through the sorted records
            for i, qs_pair in enumerate(record_order):
                qname, sname = qs_pair
                record_list = records[qs_pair]
                # Each should have only one upstream and one downstream homology region
                # Should be able to comment this condition out as it is already accounted for in the record population/filtering step
                if (len(record_list) == 2):
                    if (record_list[0].flags & 16 == record_list[1].flags & 16): # Each homology region should have identical strand orientation
                        self.logger.info('Working on {} vs {}'.format(qname, sname))
                        self.logger.info("crossover: " + str((i, qname, sname, record_list)))
                        
                        # Identify which record is upstream, and which is downstream (by query)
                        us_record, ds_record = sorted(record_list, key=lambda x: x.query_position)
                        
                        #new_name = sname+'-r'+str(r) # Add round number to the contig header
                        name_changes.setdefault(sname, []).append(qname)
                        
                        # get the entire contig sequences
                        qcontig = dDNA_contigs_list[r][qname]
                        scontig = genome_contigs_list[r][sname]
                        
                        # Get the QUERY homology coordinates
                        qush_start, qush_end = us_record.query_position[0], us_record.query_position[1]
                        qdsh_start, qdsh_end = ds_record.query_position[0], ds_record.query_position[1]
                        
                        # Get the SUBJECT homology coordinates
                        if (us_record.flags & 16): # Reverse complement means upstream/downstream are swapped
                            self.logger.info('  us_record is reverse complemented')
                            sush_start, sush_end = ds_record.subject_position[0], ds_record.subject_position[1]
                            sdsh_start, sdsh_end = us_record.subject_position[0], us_record.subject_position[1]
                        else: # not reverse-complemented
                            sush_start, sush_end = us_record.subject_position[0], us_record.subject_position[1]
                            sdsh_start, sdsh_end = ds_record.subject_position[0], ds_record.subject_position[1]
                        
                        
                        
                        ######## Get the sequences ########
                        
                        for seqloop in range(2):
                            ######## Manage QUERY sequences ########
                            
                            # Calculate new insert sequence (with homology arms) from dDNA
                            q_hih_seq = qcontig[qush_start:qdsh_end] # us-homology + insert + ds-homology (should be entire contig)
                            
                            # Calculate new insert sequence (excluding homology arms) from dDNA
                            insert_seq = qcontig[qush_end:qdsh_start] # insert
                            
                            # Get the up/down-stream homology regions for the dDNA (not the chromosome)
                            q_ush_seq = qcontig[qush_start:qush_end]
                            q_dsh_seq = qcontig[qdsh_start:qdsh_end]
                            
                            if (us_record.flags & 16): # Reverse complement if necessary
                                # We must reverse-complement these dDNA substrings so their orientation matches the genome
                                q_hih_seq = nucleotides.rc(q_hih_seq)
                                insert_seq = nucleotides.rc(insert_seq)
                                q_ush_seq, q_dsh_seq = nucleotides.rc(q_dsh_seq), nucleotides.rc(q_ush_seq)
                            
                            ######## Manage SUBJECT sequences ########
                            
                            # Extract the feature sequence (excluding homology arms)
                            feature_seq = scontig[sush_end:sdsh_start]
                            
                            # Calculate new upstream & downstream chromosome arms (swap us/ds for subject) (excludes homology arms)
                            us_seq = scontig[:sush_start]
                            ds_seq = scontig[sdsh_end:]
                            
                            # Get the up/down-stream homology regions for the chromosome (not the dDNA)
                            s_ush_seq = scontig[sush_start:sush_end]
                            s_dsh_seq = scontig[sdsh_start:sdsh_end]
                            
                            if (seqloop == 0):
                                # Overview for making homology regions mutually exclusive
                                #   No overlap
                                #     q          ......111111111111......2222222222.....           none
                                #     s      ....11111111111................222222222222.......    none
                                #   query overlap
                                #     q          ........1111111111XXX22222222.......              overlap
                                #     s    .......11111111111111............222222222222.......    none
                                #      q2        ........111111111122222222222.......
                                #      s2  .......11111111111'''............222222222222.......
                                #   subject overlap
                                #     q          ......111111111111......2222222222.....           none
                                #     s    ........11111111111111XXXXX2222222222222222....         overlap
                                #      q2        ......1111111'''''......2222222222.....  
                                #      s2  ........11111111111111222222222222222222222....
                                #   query & subject overlap
                                #     q          ........1111111111XXX22222222.......              overlap
                                #     s    ........11111111111111XXXXX2222222222222222....         overlap
                                #      q2        ........111111111122222222222.......
                                #      s2  ........11111111111111222222222222222222222....
                                
                                # if both query and subject overlap
                                if ((qdsh_start < qush_end) and (sdsh_start < sush_end)):
                                    # No fancy calculations required. Just set the new coordinates
                                    qush_end = qdsh_start
                                    sush_end = sdsh_start
                                # If query overlaps
                                elif (qdsh_start < qush_end):
                                    # Calculate the subject trim
                                    q_overlap = qcontig[qdsh_start:qush_end] # Just the overlapping subsequence of the query alignments
                                    if (us_record.flags & 16):
                                        q_overlap = nucleotides.rc(q_overlap)
                                    
                                    # Code to simplify regex
                                    if (len(q_overlap)//2 > 1):
                                        pattern_string = '(?:'+q_overlap+'){e<'+str(len(q_overlap)//2)+'}$'
                                    else:
                                        pattern_string = '(?:'+q_overlap+')$'
                                    self.logger.info('pattern_string = '+pattern_string)
                                    self.logger.info('s_ush_seq = '+s_ush_seq)
                                    
                                    m = regex.search(pattern_string, s_ush_seq, flags=regex.IGNORECASE|regex.BESTMATCH)
                                    if m:
                                        sush_end = sush_start + m.start()
                                    
                                    # Set the query upstream homology end
                                    qush_end = qdsh_start
                                # If subject overlaps
                                elif (sdsh_start < sush_end):
                                    # Calculate the query trim
                                    s_overlap = scontig[sdsh_start:sush_end]
                                    
                                    # Code to simplify regex
                                    if (len(s_overlap)//2 > 1):
                                        pattern_string = '(?:'+s_overlap+'){e<'+str(len(s_overlap)//2)+'}$'
                                    else:
                                        pattern_string = '(?:'+s_overlap+')$'
                                    self.logger.info('pattern_string = '+pattern_string)
                                    self.logger.info('q_ush_seq = '+q_ush_seq)
                                    
                                    m = regex.search(pattern_string, q_ush_seq, flags=regex.IGNORECASE|regex.BESTMATCH)
                                    if m:
                                        qush_end = qush_start + m.start()
                                    
                                    # Set the subject upstream homology end
                                    sush_end = sdsh_start
                                # If there is no overlap
                                else:
                                    # Exit this for loop because no coordinates were modified
                                    # and there is no need to re-calculate the sequences
                                    break
                                
                                # Otherwise, re-calculate all the sequences
                        
                        # Record the dDNA homology regions for later reference
                        # Data structure limitation that each dDNA can only engineer a single locus per contig (i.e. this code needs improvement)
                        #dDNA_homology_seqs[(r, qname, sname)] = (s_ush_seq, s_dsh_seq)
                        
                        #########
                        # Replace entry in 'genome_contigs' dict with the new cross-over contig
                        # This probably doesn't work with multiple targeted engineering events on the same chromosome
                        #genome_contigs_list[r][sname] = us_seq + q_hih_seq + ds_seq # <=---------------------- this should be done in the 'r' outer loop (earlier draft)
                        
                        # Keep the subject up/down-stream homology regions in case 2 engineered loci have overlapping homology arms
                        # like this:
                        #    site1            site2
                        #   ...-----iii------...........
                        #   .............------iii---...
                        # The result would be:
                        #   ...-----iii--------iii---...
                        #genome_contigs_list[r][sname] = us_seq + s_ush_seq + insert_seq + s_dsh_seq + ds_seq # For some reason, this code doesn't properly incorporate dDNA into the genome
                        genome_contigs_list[r][sname] = us_seq + q_hih_seq + ds_seq # Giving this a try
                        
                        # this messes with calling self.calculate_amplicons() for each query-subject pair.
                        
                        ##### alternate #####
                        # This code only works for a single-locus per chromosome! <=------ need to update this so it works for multiple loci per chromosome
                        #updated_contigs[sname] = us_seq + q_hih_seq + ds_seq
                        ### end alternate ###
                        
                        # Add this parsed record data to the 'group_links'
                        # Key refers to dDNA file and contig name, value refers to gDNA file, contig name, and homology regions
                        #group_links.append(Datum(r, qname, r-1, sname, sush_start, sush_end, sdsh_start, sdsh_end, us_record.query_position[1], ds_record.query_position[0]))
                        group_links.append(Datum(r, qname, r-1, sname, sush_start, sush_end, sdsh_start, sdsh_end, qush_end, qdsh_start))
                        
                else:
                    self.logger.info(qname + ' vs ' + sname + ' has ' + str(len(record_list)) + ' regions of homology (2 needed)')
            
            # Rename the contigs based on the modifications
            # This code only works for a single-locus per chromosome! <=----- need to update this so it works for multiple loci per chromosome
            for sname in name_changes:
                new_name = sname+'-r'+str(r)+'['+','.join(name_changes[sname])+']'
                genome_contigs_list[r][new_name] = genome_contigs_list[r].pop(sname)
                ##### alternate #####
                #genome_contigs.pop(sname)
                #genome_contigs[new_name] = updated_contigs[sname]
                ### end alternate ###
                
                # Add the new contig name to the contig name group
                for cg in contig_groups:
                    if sname in cg:
                        cg.append(new_name)
                        break
                else: # This 'else' statement should never run because we've already seeded 'contig_groups' with 'r0' contig names
                    contig_groups.append([sname, new_name])
                
            
            # Write the new modified genome to FASTA file
            genome_fasta_file_list.append(utils.write_merged_fasta(genome_contigs_list[r], os.path.join(args.folder, 'genome-r'+str(r)+'.fasta')))
        
        # Now that the modified genomes have been calculated, we need to look for shared DNA outside of the homology regions
        
        
        
        
        # This aligns each dDNA-rN against the engineered genome-rN[].
        # That is, it tells us the new location of each dDNA contig in the new genome
        # This is important because multiple edits (insertions/deletions) may have occurred
        # in the same contig, and the locations may have shifted.
        # We loop through the dDNAs as before
        for j, dDNA_filename in enumerate(args.dDNAs):
            # Variable to count the round of genome engineering
            r = j+1
            
            # Let the user know which round is being processed
            self.logger.info('Working on round r{}'.format(r))
            
            # Align each dDNA to the genome it produced
            #genome_index_file = args.selected_aligner.index(genome_fasta_file_list[r], os.path.basename(genome_fasta_file_list[r]), args.folder, args.processors)
            genome_index_file = args.selected_aligner.index(genome_fasta_file_list[r], self.pathstrip(genome_fasta_file_list[r]), args.folder, args.processors)
            temp_alignment = args.selected_aligner.align(dDNA_fasta_file_list[r], genome_index_file, 'dDNA-temp-r'+str(r)+'-alignment', args.folder, args.processors)
            
            # Parse the alignment to make an intermediary record data structure
            records = self.filter_alignment_records2(args, temp_alignment, dDNA_contigs_list[r])
            
            # Find the exact locations of each of the dDNAs within the engineered genomes
            # Iterate through the sorted records
            for i, qs_pair in enumerate(records):
                qname, sname = qs_pair
                record = records[qs_pair]
                
                self.logger.info('Working on {} vs {}'.format(qname, sname))
                
                scontig = genome_contigs_list[r][sname]
                qcontig = dDNA_contigs_list[r][qname]
                
                
                # Calculate new insert sequence (with homology arms) from dDNA
                q_hih_seq = qcontig[record.query_position[0]:record.query_position[1]] # us-homology + insert + ds-homology (should be entire contig)
                
                # Calculate new upstream & downstream chromosome arms (excludes homology arms)
                us_seq = scontig[:record.subject_position[0]]
                ds_seq = scontig[record.subject_position[1]:]
                
                s_hih_seq = scontig[record.subject_position[0]:record.subject_position[1]]
                
                # Do I need code here to ensure the us/ds homology regions are mutually exclusive?
                ##### insert code here? #####
                
                sush_start, sush_end = record.subject_position[0], None
                sdsh_start, sdsh_end = None, record.subject_position[1]
                ins_start, ins_end = None, None
                
                if (record.flags & 16): # Reverse complement if necessary
                    self.logger.info('  us_record is reverse complemented')
                    # We must reverse-complement these dDNA substrings so their orientation matches the genome
                    q_hih_seq = nucleotides.rc(q_hih_seq)
                else: # not reverse-complemented
                    pass
            
                #group_links.setdefault((r, qname), []).append(Datum(r, sname, sush_start, sush_end, sdsh_start, sdsh_end))
                group_links.append(Datum(r, qname, r, sname, sush_start, sush_end, sdsh_start, sdsh_end, ins_start, ins_end))
        
        # Link the loci and store in 'datum_groups'
        datum_groups = self.make_datum_groups(group_links, contig_groups, genome_contigs_list)
        
        ###### New code ######
        self.primer_queue_test(args, datum_groups, genome_contigs_list, dDNA_contigs_list)
        
        
        ######################
        
        # Find the longest common substring in the far-upstream and far-downstream regions
        pcr_regions, pcr_region_positions = self.get_far_lcs_regions(datum_groups, genome_contigs_list)
        
        # Identify the feature/insert sequence where the 'rN-oF', 'rN-oR', 'rN-iF', 'rN-iR' primers should be located
        # and store them in 'insert_seqs'
        #Insert = namedtuple('Insert', ['genome_r', 'genome_contig', 'seq', 'us_seq', 'ds_seq', 'fus_dist', 'fds_dist', 'type'])
        max_primer_length = 35
        
        
        for i, dg in enumerate(datum_groups):
            self.logger.info('Analyzing locus: {}'.format(i))
            
            # Add DatumGroup to 'output0' table
            output0.append([i, dg])
            
            fus_seq, fds_seq = pcr_regions[i]
            
            self.logger.info("Grabbing contig substrings for 'Insert' objects...")
            
            # insert_seqs = [
            #     Insert(qname='ko-dDNA', genome_r=0, genome_contig='chr1',               seq='ACGTAACA') 
            #     Insert(qname='ki-dDNA', genome_r=1, genome_contig='chr1-r1[ko]',        seq='ACGTAACA')
            #     Insert(qname='ki-dDNA', genome_r=2, genome_contig='chr1-r1[ko]-r2[ki]', seq='CGATAAGC')
            # ]
            insert_seqs = []
            
            # Go through every datum, and select all the 'before' ones
            # (because they have the full homology regions specified)
            # And use those to calculate the 'before'=feature, and 'after'=insert
            # sequences (with their up/downstream flanking regions)
            for di, datum in enumerate(dg):
                self.logger.info("Parsing 'Datum' {}".format(di))
                if (datum.ins_start != None):
                    # Add feature
                    insert_seqs.append(Insert(
                        datum.genome_r, # genome_r
                        datum.genome_contig, # genome_contig
                        genome_contigs_list[datum.genome_r][datum.genome_contig][datum.ush_end:datum.dsh_start], # seq
                        genome_contigs_list[datum.genome_r][datum.genome_contig][max(0, datum.ush_end-(max_primer_length-1)):datum.ush_end], # us_seq
                        genome_contigs_list[datum.genome_r][datum.genome_contig][datum.dsh_start:datum.dsh_start+(max_primer_length-1)], # ds_seq
                        datum.ush_end - pcr_region_positions[i][di][1], # fus_dist
                        pcr_region_positions[i][di][2] - datum.dsh_start, # fds_dist
                        'b' # type
                    ))
                    
                    # Add insert
                    insert_seqs.append(Insert(
                        datum.genome_r, # genome_r
                        datum.genome_contig, # genome_contig
                        #datum.dDNA_r, # genome_r <=---------------------- may need to change this to be the genome (not the dDNA)
                        #datum.dDNA_contig, # genome_contig <=------------ may need to change this to be the genome (not the dDNA)
                        dDNA_contigs_list[datum.dDNA_r][datum.dDNA_contig][datum.ins_start:datum.ins_end], # seq
                        dDNA_contigs_list[datum.dDNA_r][datum.dDNA_contig][max(0, datum.ins_start-(max_primer_length-1)):datum.ins_start], # us_seq
                        dDNA_contigs_list[datum.dDNA_r][datum.dDNA_contig][datum.ins_end:datum.ins_end+(max_primer_length-1)], # ds_seq
                        datum.ush_end - pcr_region_positions[i][di][1], # fus_dist
                        pcr_region_positions[i][di][2] - datum.dsh_start, # fds_dist
                        'a' # type
                    ))
                # Essentially, this current loop does the following:
                #   if (di == 0): (a 'before')
                #     Then get the original feature seq (for genome_r=0)
                #     And get the insert seq (for genome_r=1)
                #   if (di == 1), (an 'after')
                #     then skip
                #   if (di == 2): (a 'before')
                #     Then get the feature seq (for genome_r=1)
                #     And get the insert seq (for genome_r=2)
                #   if (di == 3): (an 'after)
                #     then skip
            
            # Now that we've identified the sequences where the shared PCR primers should reside ('sF' & 'sR'),
            # As well as the sequences of the feature/insert and their up/downstream flanking sequences
            # (which are used to find primers that span the junctions),
            # we do the actual cPCR calculations
            
            # Old code for 6-set:
            #initial_pair_list, final_pair_list = self.make_primer_set(args, qname, us_seq, ds_seq, insert_seq, feature_seq, q_hih_seq, s_ush_seq, q_ush_seq, s_dsh_seq, q_dsh_seq)
            self.logger.info('Calculating primers for locus {}'.format(i))
            #                                                                             shared_forward  shared_reverse  features/inserts
#            pair_list, insert_pair_list, round_labels, weighted_d_set_list = self.calculate_them_primers(args, fus_seq,        fds_seq,        insert_seqs)
            pair_list, insert_pair_list, round_labels, weighted_d_set_list = [], [], [], []
            
            # Filter 'pair_list' to get the top 10
            self.logger.info("Sorting Amp 'A/B/C' list...")
            pair_list = sorted(pair_list, reverse=True)[:args.max_number_designs_reported]
            
            # Make Amp-D pp list non-redundant
            self.logger.info("Sorting Amp 'D' list...")
            ppdict = {}
            for wdsl in weighted_d_set_list:
                w, pp_list = wdsl
                k = []
                for pp in pp_list:
                    if pp:
                        k.append((pp.forward_primer.sequence, pp.reverse_primer.sequence))
                    else:
                        k.append((None, None))
                k = tuple(k)
                
                d_wdsl = ppdict.setdefault(k, wdsl)
                if (d_wdsl[0] < w):
                    ppdict[k] = wdsl
                
            weighted_d_set_list = list(ppdict.values())
            
            # Filter to get the top 10
            weighted_d_set_list = sorted(weighted_d_set_list, reverse=True)[:args.max_number_designs_reported]
            
            # Let user know program is still working
            self.logger.info("Appending Locus {} results to output lists...".format(i))
            
            # Add the calculated weights of the sF/sR/oF/oR sets to 'output1' table
            self.logger.info("Appending 'output1'...")
            for ppli, (w, pp_list) in enumerate(pair_list):
                output1.append([ppli, i, w, 'A/B/C'])
            
            for ppli, (w, pp_list) in enumerate(weighted_d_set_list):
                output1.append([ppli, i, w, 'D'])
            
            # Print the primers for 'output2' table
            self.logger.info("Appending 'output2'...")
            for ppli, (w, pp_list) in enumerate(pair_list):
                round_labels_iter = iter(round_labels)
                
                for pp_i, pp in enumerate(pp_list):
                    amp_name, round_n = next(round_labels_iter)
                    #amp_name = '-'
                    #round_n = '-'
                    
                    ra_counter = self.round_counter()
                    for ia in range(len(insert_pair_list)):
                        rac = next(ra_counter)
                        converted_genome_round = self.round_converter(rac)
                        
                        if pp:
                            amp_sizes = pp.in_silico_pcr(genome_contigs_list[converted_genome_round])
                            if (len(amp_sizes) > 0):
                                amp_sizes = ','.join(map(str, sorted(amp_sizes)))
                            else:
                                amp_sizes = '-'
                            tms = ','.join([str(round(ptm, 2)) for ptm in pp.get_tms()])
                            output2.append([ppli, pp_i, i, amp_name, rac, converted_genome_round, pp.forward_primer.get_name(), pp.reverse_primer.get_name(), amp_sizes, tms])
                        
                        #if (pp and (round_n == rac)):
                        #    output2.append([ppli, pp_i, i, amp_name, round_n, '-', pp.forward_primer.name, pp.reverse_primer.name, pp.get_amplicon_size(), pp.get_tms()])
                        else:
                            if (amp_name == 'A'):
                                fname, rname = 'sF', 'sR'
                            elif (amp_name == 'B'):
                                fname, rname = 'sF', 'r'+round_n+'-oR'
                            elif (amp_name == 'C'):
                                fname, rname = 'r'+round_n+'-oF', 'sR'
                            else:
                                fname, rname = '-', '-'
                            output2.append([ppli, pp_i, i, amp_name, rac, converted_genome_round, fname, rname, '-', '-'])
            
            amp_name = 'D'
            for ppli, (w, pp_list) in enumerate(weighted_d_set_list):
                d_counter = self.round_counter()
                for pp_i, pp in enumerate(pp_list):
                    round_n = next(d_counter)
                    
                    ra_counter = self.round_counter()
                    for ip_r in range(len(insert_pair_list)):
                        rac = next(ra_counter)
                        converted_genome_round = self.round_converter(rac)
                        if pp:
                            amp_sizes = pp.in_silico_pcr(genome_contigs_list[converted_genome_round])
                            if (len(amp_sizes) > 0):
                                amp_sizes = ','.join(map(str, sorted(amp_sizes)))
                            else:
                                amp_sizes = '-'
                            tms = ','.join([str(round(ptm, 2)) for ptm in pp.get_tms()])
                            output2.append([ppli, pp_i, i, amp_name, rac, converted_genome_round, pp.forward_primer.get_name(), pp.reverse_primer.get_name(), amp_sizes, tms])
                            
                        #if (pp and (round_n == rac)):
                        #    output2.append([ppli, ip_r, i, amp_name, round_n, '-', pp.forward_primer.name, pp.reverse_primer.name, pp.get_amplicon_size(), pp.get_tms()])
                        else:
                            output2.append([ppli, pp_i, i, amp_name, rac, converted_genome_round, 'r'+round_n+'-iF', 'r'+round_n+'-iR', '-', '-'])
                    
            # Add Primer sequences to 'output3' table
            self.logger.info("Appending 'output3'...")
            for ppli, (w, pp_list) in enumerate(pair_list):
                round_labels_iter = iter(round_labels)
                #non_redundant_primers = {}
                for pp_i, pp in enumerate(pp_list):
                    amp_name, round_n = next(round_labels_iter)
                    if pp:
                        fname, fseq = pp.forward_primer.get_name(), pp.forward_primer.sequence
                        rname, rseq = pp.reverse_primer.get_name(), pp.reverse_primer.sequence
                    else:
                        if (amp_name == 'A'):
                            fname, fseq = 'sF', '-'
                            rname, rseq = 'sR', '-'
                        elif (amp_name == 'B'):
                            fname, fseq = 'sF', '-'
                            rname, rseq = 'r'+round_n+'-oR', '-'
                        elif (amp_name == 'C'):
                            fname, fseq = 'r'+round_n+'-oF', '-'
                            rname, rseq = 'sR', '-'
                        else:
                            fname, fseq, rname, rseq = '-', '-', '-', '-'
                    #non_redundany_primers[fname] = fseq
                    #non_redundany_primers[rname] = rseq
                    output3.append([ppli, pp_i, i, fname, fseq])
                    output3.append([ppli, pp_i, i, rname, rseq])
            
            #amp_name = 'D'
            for ppli, (w, pp_list) in enumerate(weighted_d_set_list):
                d_counter = self.round_counter()
                for pp_i, pp in enumerate(pp_list):
                    round_n = next(d_counter)
                    #converted_genome_round = self.round_converter(round_n)
                    
                    if pp:
                        fname, fseq = pp.forward_primer.get_name(), pp.forward_primer.sequence
                        rname, rseq = pp.reverse_primer.get_name(), pp.reverse_primer.sequence
                    else:
                        fname, fseq = 'r'+round_n+'-iF', '-'
                        rname, rseq = 'r'+round_n+'-iR', '-'
                    
                    output3.append([ppli, pp_i, i, fname, fseq])
                    output3.append([ppli, pp_i, i, rname, rseq])
            
            # Add the primer locations to 'output4' table
            self.logger.info("Appending 'output4'...")
            for ppli, (w, pp_list) in enumerate(pair_list):
                round_labels_iter = iter(round_labels)
                for pp_i, pp in enumerate(pp_list):
                    amp_name, round_n = next(round_labels_iter)
                    if pp:
                        # 'f_locations' and 'r_locations' should have the same length
                        f_locations = self.get_primer_location(genome_fasta_file_list, genome_contigs_list, pp.forward_primer.sequence)
                        r_locations = self.get_primer_location(genome_fasta_file_list, genome_contigs_list, pp.reverse_primer.sequence)
                        for loc_i in range(len(f_locations)):
                            for loc in f_locations[loc_i]:
                                output4.append([ppli, pp_i, i, pp.forward_primer.get_name(), pp.forward_primer.sequence, loc[0], loc[1], loc[2], loc[3], loc[4]])
                            for loc in r_locations[loc_i]:
                                output4.append([ppli, pp_i, i, pp.reverse_primer.get_name(), pp.reverse_primer.sequence, loc[0], loc[1], loc[2], loc[3], loc[4]])
                        #output4.append([ppli, pp_i, i, pp.reverse_primer.name, pp.reverse_primer.sequence, 'File', 'Contig', 'Start', 'End', '-'])
                    else:
                        if (amp_name == 'A'):
                            fname, fseq = 'sF', '-'
                            rname, rseq = 'sR', '-'
                        elif (amp_name == 'B'):
                            fname, fseq = 'sF', '-'
                            rname, rseq = 'r'+round_n+'-oR', '-'
                        elif (amp_name == 'C'):
                            fname, fseq = 'r'+round_n+'-oF', '-'
                            rname, rseq = 'sR', '-'
                        
                        output4.append([ppli, pp_i, i, fname, fseq, '-', '-', '-', '-', '+'])
                        output4.append([ppli, pp_i, i, fname, fseq, '-', '-', '-', '-', '-'])
            
            for ppli, (w, pp_list) in enumerate(weighted_d_set_list):
                d_counter = self.round_counter()
                for pp_i, pp in enumerate(pp_list):
                    round_n = next(d_counter)
                    if pp:
                        # 'f_locations' and 'r_locations' should have the same length
                        f_locations = self.get_primer_location(genome_fasta_file_list, genome_contigs_list, pp.forward_primer.sequence)
                        r_locations = self.get_primer_location(genome_fasta_file_list, genome_contigs_list, pp.reverse_primer.sequence)
                        for loc_i in range(len(f_locations)):
                            for loc in f_locations[loc_i]:
                                output4.append([ppli, pp_i, i, pp.forward_primer.get_name(), pp.forward_primer.sequence, loc[0], loc[1], loc[2], loc[3], loc[4]])
                            for loc in r_locations[loc_i]:
                                output4.append([ppli, pp_i, i, pp.reverse_primer.get_name(), pp.reverse_primer.sequence, loc[0], loc[1], loc[2], loc[3], loc[4]])
                    else:
                        output4.append([ppli, pp_i, i, 'r'+round_n+'-iF', '-', '-', '-', '-', '-', '+'])
                        output4.append([ppli, pp_i, i, 'r'+round_n+'-iR', '-', '-', '-', '-', '-', '-'])
        
        self.logger.info('Printing to STDOUT...')
        
        # Print output header information describing the amplicons
        print('# Amplicon diagram')
        print('#                                          Genome')
        print('# Amplicon  ──upstream─┐┌─homology─┐┌──insert/feature──┐┌─homology─┐┌─downstream──')
        print('#        A   sF ===>····················································<=== sR')
        print('#        B   sF ===>···················<=== rN-oR')
        print('#        C                             rN-oF ===>·······················<=== sR')
        print('#        D                        rN-iF ===>······<=== rN-iR')
        print('#')
        print('# F = forward')
        print('# R = reverse')
        print('# s = shared')
        print('# o = outer')
        print('# i = inner')
        print('# rN = round number')
        
        # Print the various output tables
        # TODO: Put Locus/DatumGroup calculation into its own function:
        #       output0 = calc_loci(args, ...)
        print("# In silico recombination")
        print('# ' + '\t'.join(['Locus', 'DatumGroup']))
        for line in output0:
            print('\t'.join(map(str, line)))
        
        print('# ' + '\t'.join(['Set', 'Locus', 'Weight', 'Amplicon']))
        for line in output1:
            print('\t'.join(map(str, line)))
        
        # TODO: Put in-silico PCR (output2) into its own function
        print('# ' + '\t'.join(['Set', 'Index', 'Locus', 'Amplicon', 'Template', 'Genome', 'F', 'R', 'Sizes', 'Tm']))
        for line in output2:
            print('\t'.join(map(str, line)))
        
        print('# ' + '\t'.join(['Set', 'Index', 'Locus', 'Primer', 'Sequence']))
        for line in output3:
            print('\t'.join(map(str, line)))
        
        print('# ' + '\t'.join(['Set', 'Index', 'Locus', 'Primer', 'Sequence', 'File', 'Contig', 'Start', 'End', 'Strand']))
        for line in output4:
            print('\t'.join(map(str, line)))
        
                        #####################
                        # Cut code was here #
                        #####################
                        
                        # This code output GenBank '*.gb' files
        
        # End 'compute()'
    
    def primer_queue_test_old(self, args, loci, genome_contigs_list, dDNA_contigs_list):
        '''
        New pseudocode for primer set calculations
        needs to handle:
          allele-specific, multi-allele, allele-agnostic
          masked characters
        '''
        
        self.logger.info("Starting 'primer_queue_test()'...")
        # This should be input to this function
        # loci = datum_list/datum_groups
        #loci = ['hapA', 'hapB']
        
        FAR_UPSTREAM = 1
        FAR_DOWNSTREAM = 2
        FEATURE = 4
        INSERT = 8
        
        pairs = [
            ((FAR_UPSTREAM, '+', 'sF'), (FAR_DOWNSTREAM, '-', 'sR')), # sF sR = Amp A
            
            ((FAR_UPSTREAM, '+', 'sF'), (FEATURE, '-', 'oR')),        # sF oR = Amp B
            ((FEATURE, '+', 'oF'), (FAR_DOWNSTREAM, '-', 'sR')),      # oF sR = Amp C
            ((FEATURE, '+', 'iF'), (FEATURE, '-', 'iR')),             # iF iR = Amp D
            
            ((FAR_UPSTREAM, '+', 'sF'), (INSERT, '-', 'oR')),         # sF oR = Amp B
            ((INSERT, '+', 'oF'), (FAR_DOWNSTREAM, '-', 'sR')),       # oF sR = Amp C
            ((INSERT, '+', 'iF'), (INSERT, '-', 'iR')),               # iF iR = Amp D
        ]
        regions = set(x for y in pairs for x in y)
        
        # Populate the 'Primer.sequences' dict to serve as a queue (but don't do any calculations yet)
        # Scan each region, and add the sequence as a 'Primer' object
        for locus, dg in enumerate(loci):
            for datum in dg:
                if (datum.ins_start != None):
                    for region, orientation, name in regions:
                        if (region == FAR_UPSTREAM):
                            start, end = max(0, datum.ush_start-1300), datum.ush_start
                            sequence = genome_contigs_list[datum.genome_r][datum.genome_contig]
                            contig = datum.genome_contig
                            nname = name
                        
                        elif (region == FAR_DOWNSTREAM):
                            start, end = datum.dsh_end, min(datum.dsh_end+1300, len(genome_contigs_list[datum.genome_r][datum.genome_contig]))
                            sequence = genome_contigs_list[datum.genome_r][datum.genome_contig]
                            contig = datum.genome_contig
                            nname = name
                        
                        elif (region == FEATURE):
                            start, end = datum.ush_end, datum.dsh_start
                            sequence = genome_contigs_list[datum.genome_r][datum.genome_contig]
                            contig = datum.genome_contig
                            nname = 'r'+str(datum.genome_r)+'-'+name
                        
                        elif (region == INSERT):
                            start, end = datum.ins_start, datum.ins_end
                            sequence = dDNA_contigs_list[datum.dDNA_r][datum.dDNA_contig]
                            contig = datum.dDNA_contig
                            nname = 'r'+str(datum.genome_r+1)+'-'+name
                                
                        #args.selected_oligo.scan(ins.seq, 'left',  primer_size=primer_length_range, folder=temp_folder, us_seq=ins.us_seq, ds_seq=ins.ds_seq, time_limit=args.primer_scan_limit)
                        #args.selected_oligo.scan(seq, l, r, o)
                        nname = None
                        Primer.scan(sequence, 0, locus, datum.genome_r, region, contig, orientation, start, end, nname, primer_size=(19,36))
        
        self.logger.info("Total 'Primer' objects before filtering: {}".format(len(Primer.sequences)))
        
        # Determine the minimal set of loci/genomes/contigs to constitute correct specificity
        # Filter the primer queue in 'Primer.sequences' so only valid primers remain
        self.logger.info("Removing 'Primer' objects that don't meet desired: specificity='{}'".format(args.specificity))
        
        if (args.specificity == 'exclusive'): # allele-specific
            # Primer must be present in only one locus and not any others
            new_primer_sequences = {}
            for seq, p in Primer.sequences.items():
                if (len(p.get_specificity()) == 1): # Assumes non-multiplex design
                    new_primer_sequences[seq] = p
            Primer.sequences = new_primer_sequences
            
        elif (args.specificity == 'all'): # multi-allelic
            # Primer must be present in all loci
            new_primer_sequences = {}
            for seq, p in Primer.sequences.items():
                if (len(p.get_specificity()) == len(loci)):
                    new_primer_sequences[seq] = p
            Primer.sequences = new_primer_sequences
            
        elif (args.specificity == 'any'): # allele-agnostic
            # Primer can be present in any number of loci
            pass
        
        self.logger.info("Total 'Primer' objects after filtering: {}".format(len(Primer.sequences)))
        
        # Set up 'PrimerPair' objects, and place them in a queue 'PrimerPair.pairs'
        self.logger.info("Setting up 'PrimerPair' objects")
        
        
        clist = [
            #   Name, initial, final, delta # states
            Cutoff("length", (19,28), (19,36), (0,2), separate=True), # 1,5
            Cutoff("last5gc_count", (1,3), (0,4), (-1,1), separate=True), # 2, last5gc
            Cutoff("gc_clamp_length", (1,2), (0,4), (-1,1), separate=True), # 2, gcclamp
            Cutoff("gc", (0.4, 0.6), (0.2, 0.8), (-0.1, 0.1)), # 3, gcfreq
            Cutoff("max_run_length", (4,), (6,), (1,)), # 3, runlen_max
            Cutoff("max_3prime_complementation_length", (3,), (6,), (1,)), # 4, 3primecomplen_max
            Cutoff("min_delta_g", (-4.0,), (-7.0,), (-1.0,)), # 4, deltag_min
            Cutoff("tm", (52,65), (52,65), (0,0)), # 1, tm
            Cutoff("max_tm_difference", (2.5,), (4.0,), (0.5,)), # 4, deltatm_max
        ]
        
        citer = CutoffIterator()
        
        design_found = False
        # Do calculations until either:
        #  * all primers are assayed,
        #  * an adequate design is found,
        #  * or the time limit expires
        cycle_n = 0
        while (design_found == False):
            break # added so these calculations won't be run
            try:
                # Get the initial cutoffs (if this is the first loop)
                # or get the next-most relaxed cutoffs (if this isn't the first loop)
                cutoffs = next(citer)
            except StopIteration:
                break
            cutoffs['o_oligo'] = args.selected_oligo
            cutoffs['folder'] = os.path.join(args.temp_folder, 'addtag', os.path.basename(args.folder))
            
            self.logger.info("cycle {}".format(cycle_n))
            self.logger.info("  cutoffs = {}".format(cutoffs))
            
            # Do 'Primer' calculations until the time limit is reached
            p_tot = 0
            p_checked_now = 0
            p_passed_previously = 0
            p_rejected = 0
            p_skipped = 0
            start_time = time.time()
            time_expired = False
            #while (time.time()-start_time < args.primer_scan_limit):
            #    p = primer_queue.get()
            #    # First, perform cheap calculation on the primers
            #    # If there is time remaining, then perform the next-most expensive calculation
            #    # If there is time remaining, then perform the most expensive calculations
            #    # If there are no good primers, then 'relax' the thresholds
            #    # Then do calculations using the 'relaxed' thresholds
            for pi, (seq, p) in enumerate(Primer.sequences.items()):
                if (pi % 1000 == 0):
                    if ((time.time()-start_time) > args.primer_scan_limit):
                        time_expired = True
                
                if not time_expired:
                    cpass = p.summarize(p.checks)
                    if ((p.checks[0] == None) or (not cpass)):
                        p.progressive_check(cutoffs)
                        p_checked_now += 1
                    elif cpass:
                        p_passed_previously += 1
                    else:
                        p_rejected += 1
                else:
                    p_skipped += 1
                
                p_tot += 1
            
            self.logger.info("  'Primer' objects: checked_now={}, passed_previously={}, not_checked={}, skipped={}, total={}".format(p_checked_now, p_passed_previously, p_rejected, p_skipped, p_tot))
            
            
            ##### Some debug code #####
            sss = 0
            ttt = 0
            nnn = 0
            hhh = 0
            ddd = 0
            ooo = 0
            for seq, p in Primer.sequences.items():
                ttt += 1
                if p.summarize(p.checks):
                    sss += 1
                if ((p.checks[0] != None) and p.summarize(p.checks)):
                    nnn += 1
                if (p.o_hairpin != None):
                    hhh += 1
                if (p.o_self_dimer != None):
                    ddd += 1
                if (p.o_reverse_complement != None):
                    ooo += 1
            self.logger.info("   sss={}, nnn={}, ttt={}, hhh={}, ddd={}, ooo={}".format(sss, nnn, ttt, hhh, ddd, ooo))
            ##### Some debug code #####
            
            
            
            subset_size = 1000
            self.logger.info("  Queueing 'PrimerPair' objects")
            #pair_queue = []
            # pair_queue = [
            #     [                        # Locus 0
            #         [PrimerPair(), ...], # ((FAR_UPSTREAM, '+', 'sF'), (FAR_DOWNSTREAM, '-', 'sR')),
            #         [PrimerPair(), ...], # ((FAR_UPSTREAM, '+', 'sF'), (FEATURE, '-', 'oR'))
            #         ...
            #     ],
            #     ...
            # ]
            
            
            
            # BEGIN Try 3/27/2019
            pair_queue = []
            
            for (f_reg, f_ori, f_name), (r_reg, r_ori, r_name) in pairs:
                self.logger.info("  pair: F={}, R={}".format((f_reg, f_ori, f_name), (r_reg, r_ori, r_name)))
                f_list = []
                r_list = []
                #loci_needed = list(range(len(loci)))
                for seq, p in Primer.sequences.items():
                    if ((p.checks[0] != None) and p.summarize(p.checks)):
                        plocs = [loc[0:2]+loc[3:4] for loc in p.locations] # [locus, region, strand]
                        
                        # Assume multi-allelic
                        f_in = [False] * len(loci)
                        r_in = [False] * len(loci)
                        for locus, dg in enumerate(loci):
                            if ((locus, f_reg, f_ori) in plocs):
                                f_in[locus] = True
                            #else:
                            #    break # Already failed, so might as well stop loop to save computations
                            if ((locus, r_reg, r_ori) in plocs):
                                r_in[locus] = True
                            #else:
                            #    break # Already failed, so might as well stop loop to save computations
                                
                        if all(f_in):
                            f_list.append(seq)
                            p.set_name('r?-'+f_name)
                        if all(r_in):
                            r_list.append(seq)
                            p.set_name('r?-'+r_name)
                self.logger.info("    len(f_list)={}, len(r_list)={}".format(len(f_list), len(r_list)))
                pp_list = []
                for p1 in sorted(f_list, reverse=True)[:subset_size]:
                    for p2 in sorted(r_list, reverse=True)[:subset_size]:
                        #if (p1, p2) not in pp_list: # This is implicit
                        PrimerPair(Primer.sequences[p1], Primer.sequences[p2])
                        
                        pair = (p1, p2)
                        pp = PrimerPair.pairs[pair]
                        pp.progressive_check(cutoffs)
                        
                        if ((pp.checks[0] != None) and Primer.summarize(pp.checks)):
                            pp_list.append(pair)
                pair_queue.append(pp_list)
                self.logger.info("    Added {} 'PrimerPair' objects for pair: F={}, R={}".format(len(pp_list), (f_reg, f_ori, f_name), (r_reg, r_ori, r_name)))
            
            
            
            # Populate sF_sR_paired_primers with primer objects
            sF_sR_paired_primers = []
            
            for pair in pair_queue[0]:
                pp = PrimerPair.pairs[pair]
                #pp.weight = pp.get_weight(locus, FAR_UPSTREAM, FAR_DOWNSTREAM, contig, minimize=True)
                pp.weight = pp.get_weight(minimize=True)
                sF_sR_paired_primers.append(pp)
            
            # Go through the 'pair_queue' one sF sR pair at a time
            design_count = 0
            for loop_i, sF_sR_pair in enumerate(sorted(sF_sR_paired_primers, reverse=True)):
                self.logger.info('  loop {}:'.format(loop_i))
                pp_sources = []
                pp_sources.append([sF_sR_pair])
                
                # Order in 'pairs' variable
                # A,   B,C,D,   B,C,D
                #for ppi, pp_list in pp2d_list:
                for ppi, pp_seq_list in enumerate(pair_queue[1:]):
                    pp_list = []
                    for pair in pp_seq_list:
                        pp_list.append(PrimerPair.pairs[pair])
                    
                    if (ppi in [0, 3]): # B
                        pp_sources.append(self.filter_primer_pairs(pp_list, forward=sF_sR_pair.forward_primer))
                    elif (ppi in [1, 4]): # C
                        pp_sources.append(self.filter_primer_pairs(pp_list, reverse=sF_sR_pair.reverse_primer))
                    elif (ppi in [2, 5]): # D
                        pp_sources.append(pp_list)
                
                current_pp_sources = [len(x) for x in pp_sources]
                self.logger.info('  length of sources: {}'.format(current_pp_sources))
                
                
                # Ideally, similar code should be executed, but before the 'PrimerPair.progressive_check()' is calculated
                #required_pattern = [val for val in args.internal_primers_required for b in range(2)]
                required_pattern = ['y', 'y', 'y', 'n', 'y', 'y', 'n']
                should_mask = False
                should_skip = False
                for req_str, num_pp in zip(required_pattern, current_pp_sources):
                    if ((req_str in ['y', 'Y', '1', 'T', 't', 'TRUE', 'True', 'true']) and (num_pp == 0)):
                        should_skip = True
                        break
                    if ((req_str not in ['y', 'Y', '1', 'T', 't', 'TRUE', 'True', 'true']) and (num_pp > 0)):
                        should_mask = True
                if should_skip:
                    self.logger.info('  skipping...')
                    continue
                if should_mask:
                    self.logger.info('  masking...')
                    continue
                
                
                
                
                # Evaluate if any new primer designs are adequate
                design = PrimerDesign(pp_sources)
                #design.optimize(mode='direct')
                design.optimize()
                optimal = design.optimal
                
                design_count += 1
                
                if (design_count >= subset_size):
                    break
            
            # END Try 3/27/2019
            
            cycle_n += 1
            
        
        # Output the results, starting with the best-found design
        #for locus in datum_list:
        #    sets = calculate_primer_sets()
        #    for s in sets[:max_sets]:
        #        print(s)
        self.logger.info("Function 'primer_queue_test()' completed.")
        
    
    def primer_queue_test(self, args, loci, genome_contigs_list, dDNA_contigs_list):
        '''
        New Code for primer set calculations
        needs to handle:
          allele-specific, multi-allele, allele-agnostic (not yet) (it is hard-coded to multi-allele)
          case-masked characters (not yet)
          ambiguous-masked characters (not yet)
          disambiguation of ambiguous characters (not yet) (This can be handled in the Primer object by taking averages of metrics for each disambiguated sequence)
          masked amplicons (not yet)
        '''
        
        self.logger.info("Starting 'primer_queue_test()'...")
        
        #### BEGIN 3/28 ####
        
        # Goal is to re-jiggle Datum and Insert objects so that the data structure is easy to parse
        # so we want:
        #  the GENE is defined by the dDNA 'contig' names
        #    dDNA1a => dDNA2a => dDNA3a
        #  every GENE can align to multiple loci within each genome/round (r0, r1, r2, etc)
        #  each locus has 3 regions, and 4 region+orientations:
        #    upstream F, downstream R, feature/insert F, feature/insert R
        #  upstream F and downstream R must have shared primers among ALL loci for each gene across ALL genome/rounds
        #  feature/insert F and feature/insert R primers must be shared among all loci on a genome/round-specific manner
        #    i.e. dDNA1+locus0+g0+featureF = dDNA1_locus1+g0+featureF
        #
        # let's try a 2D representation of this
        #  GENE           locus genome region              contig   start  end  orientation  name   group         
        #  dDNA1 => dDNA2 0     r0     upstream F          ...      ...    ...  +            sF       <<<           
        #  dDNA1 => dDNA2 0     r0     downstream R        ...      ...    ...  -            sR            >>>      
        #  dDNA1 => dDNA2 0     r0     feature/insert F    ...      ...    ...  +                               0  
        #  dDNA1 => dDNA2 0     r0     feature/insert R    ...      ...    ...  -                                 0
        #  dDNA1 => dDNA2 0     r1     upstream F          ...      ...    ...  +            sF       <<<           
        #  dDNA1 => dDNA2 0     r1     downstream R        ...      ...    ...  -            sR            >>>      
        #  dDNA1 => dDNA2 0     r1     feature/insert F    ...      ...    ...  +                               1  
        #  dDNA1 => dDNA2 0     r1     feature/insert R    ...      ...    ...  -                                 1
        #  dDNA1 => dDNA2 0     r2     upstream F          ...      ...    ...  +            sF       <<<           
        #  dDNA1 => dDNA2 0     r2     downstream R        ...      ...    ...  -            sR            >>>      
        #  dDNA1 => dDNA2 0     r2     feature/insert F    ...      ...    ...  +                               2  
        #  dDNA1 => dDNA2 0     r2     feature/insert R    ...      ...    ...  -                                 2
        #  dDNA1 => dDNA2 1     r0     upstream F          ...      ...    ...  +            sF       <<<           
        #  dDNA1 => dDNA2 1     r0     downstream R        ...      ...    ...  -            sR            >>>      
        #  dDNA1 => dDNA2 1     r0     feature/insert F    ...      ...    ...  +                               0  
        #  dDNA1 => dDNA2 1     r0     feature/insert R    ...      ...    ...  -                                 0
        #  dDNA1 => dDNA2 1     r1     upstream F          ...      ...    ...  +            sF       <<<           
        #  dDNA1 => dDNA2 1     r1     downstream R        ...      ...    ...  -            sR            >>>      
        #  dDNA1 => dDNA2 1     r1     feature/insert F    ...      ...    ...  +                              1   
        #  dDNA1 => dDNA2 1     r1     feature/insert R    ...      ...    ...  -                                 1
        #  dDNA1 => dDNA2 1     r2     upstream F          ...      ...    ...  +            sF       <<<           
        #  dDNA1 => dDNA2 1     r2     downstream R        ...      ...    ...  -            sR            >>>      
        #  dDNA1 => dDNA2 1     r2     feature/insert F    ...      ...    ...  +                              2   
        #  dDNA1 => dDNA2 1     r2     feature/insert R    ...      ...    ...  -                                 2
        
        
        # datum_groups (finished):
        # [
        #   Datum(dDNA_r=1, dDNA_contig='exDonor-5', genome_r=0, genome_contig='Ca22chr4B_C_albicans_SC5314', ush_start=1044509, ush_end=1044547, dsh_start=1044673, dsh_end=1044712, ins_start=38, ins_end=61), 
        #   Datum(dDNA_r=1, dDNA_contig='exDonor-5', genome_r=1, genome_contig='Ca22chr4B_C_albicans_SC5314-r1[exDonor-5]', ush_start=1044509, ush_end=None, dsh_start=None, dsh_end=1044609, ins_start=None, ins_end=None), 
        #   Datum(dDNA_r=2, dDNA_contig='reDonor-0', genome_r=1, genome_contig='Ca22chr4B_C_albicans_SC5314-r1[exDonor-5]', ush_start=1044372, ush_end=1044547, dsh_start=1044570, dsh_end=1044737, ins_start=172, ins_end=299), 
        #   Datum(dDNA_r=2, dDNA_contig='reDonor-0', genome_r=2, genome_contig='Ca22chr4B_C_albicans_SC5314-r1[exDonor-5]-r2[reDonor-0]', ush_start=1044372, ush_end=None, dsh_start=None, dsh_end=1044838, ins_start=None, ins_end=None)
        # ],
        # [
        #   Datum(dDNA_r=1, dDNA_contig='exDonor-5', genome_r=0, genome_contig='Ca22chr4A_C_albicans_SC5314', ush_start=1044481, ush_end=1044519, dsh_start=1044646, dsh_end=1044685, ins_start=38, ins_end=61), 
        #   Datum(dDNA_r=1, dDNA_contig='exDonor-5', genome_r=1, genome_contig='Ca22chr4A_C_albicans_SC5314-r1[exDonor-5]', ush_start=1044481, ush_end=None, dsh_start=None, dsh_end=1044581, ins_start=None, ins_end=None), 
        #   Datum(dDNA_r=2, dDNA_contig='reDonor-0', genome_r=1, genome_contig='Ca22chr4A_C_albicans_SC5314-r1[exDonor-5]', ush_start=1044347, ush_end=1044519, dsh_start=1044542, dsh_end=1044709, ins_start=172, ins_end=299), 
        #   Datum(dDNA_r=2, dDNA_contig='reDonor-0', genome_r=2, genome_contig='Ca22chr4A_C_albicans_SC5314-r1[exDonor-5]-r2[reDonor-0]', ush_start=1044347, ush_end=None, dsh_start=None, dsh_end=1044813, ins_start=None, ins_end=None)
        # ]
        
        self.logger.info('Starting test of new 2D simplification of the data')
        
        # Define regions as constants for low-memory reference variables
        FAR_UPSTREAM = 1
        FAR_DOWNSTREAM = 2
        FEATURE_F = 4
        FEATURE_R = 8
        
        # Text representation of the regions (for debug output)
        region_decode = {
            FAR_UPSTREAM: 'FAR_UPSTREAM',
            FAR_DOWNSTREAM: 'FAR_DOWNSTREAM',
            FEATURE_F: 'FEATURE_F',
            FEATURE_R: 'FEATURE_R',
        }
        
        search_distance = 800
        
        data = []
        for i, dg in enumerate(loci):
            gene = ','.join(sorted(set(datum.dDNA_contig for datum in dg)))
            locus = i
            temp = None
            temp2 = None
            for datum in dg:
                genome = datum.genome_r
                contig = datum.genome_contig
                
                if (datum.dDNA_r == datum.genome_r):
                    temp = [datum.dDNA_r, datum.genome_r, datum.ush_start, datum.dsh_end]
                    
                    # The last one
                    if (datum.genome_r == len(genome_contigs_list)-1):
                        
                        region = FAR_UPSTREAM
                        start = max(0, datum.ush_start - search_distance)
                        end = datum.ush_start
                        strand = '+'
                        data.append([gene, locus, genome, region, contig, strand, start, end])
                        
                        region = FAR_DOWNSTREAM
                        start = datum.dsh_end
                        end = min(datum.dsh_end+search_distance, len(genome_contigs_list[genome][contig]))
                        strand = '-'
                        data.append([gene, locus, genome, region, contig, strand, start, end])
                        
                        region = FEATURE_F
                        start = datum.ush_start+temp2[2]
                        end = datum.ush_start+temp2[3]
                        strand = '+'
                        data.append([gene, locus, genome, region, contig, strand, start, end])
                        
                        region = FEATURE_R
                        start = datum.ush_start+temp2[2]
                        end = datum.ush_start+temp2[3]
                        strand = '-'
                        data.append([gene, locus, genome, region, contig, strand, start, end])
                
                elif (datum.dDNA_r == datum.genome_r+1):
                    
                    region = FAR_UPSTREAM
                    if (temp and (temp[0] == datum.dDNA_r-1) and (temp[1] == datum.genome_r)):
                        start = min(max(0, temp[2]-search_distance), max(0, datum.ush_start-search_distance))
                        end = min(temp[2], datum.ush_start)
                    else:
                        start = max(0, datum.ush_start - search_distance)
                        end = datum.ush_start
                    strand = '+'
                    data.append([gene, locus, genome, region, contig, strand, start, end])
                    
                    region = FAR_DOWNSTREAM
                    if (temp and (temp[0] == datum.dDNA_r-1) and (temp[1] == datum.genome_r)):
                        start = max(temp[3], datum.dsh_end)
                        end = max(min(temp[3]+search_distance, len(genome_contigs_list[genome][contig])), min(datum.dsh_end+search_distance, len(genome_contigs_list[genome][contig])))
                    else:
                        start = datum.dsh_end
                        end = min(datum.dsh_end+search_distance, len(genome_contigs_list[genome][contig]))
                    strand = '-'
                    data.append([gene, locus, genome, region, contig, strand, start, end])
                    
                    region = FEATURE_F
                    start = datum.ush_end
                    end = datum.dsh_start
                    strand = '+'
                    data.append([gene, locus, genome, region, contig, strand, start, end])
                    
                    region = FEATURE_R
                    start = datum.ush_end
                    end = datum.dsh_start
                    strand = '-'
                    data.append([gene, locus, genome, region, contig, strand, start, end])
                    
                    temp2 = [datum.dDNA_r, datum.genome_r, datum.ins_start, datum.ins_end]
        
        ##### BEGIN OUTPUT 1 #####
        # Prints the REGIONS where the primers will be sought
        print('# Region definitions')
        #                       0        1       2         (2)                 3         (3)                 4         5         6        7
        print('# ' + '\t'.join(['Gene', 'Locus', 'Genome', 'Genome (decoded)', 'Region', 'Region (decoded)', 'Contig', 'Strand', 'Start', 'End']))
        for d in data:
            print('\t'.join(list(map(str, d[:3])) + ['genome-r{}.fasta'.format(d[2]), str(d[3]), region_decode[d[3]]] + list(map(str, d[4:]))), flush=True)
            #print('\t'.join(map(str, d)), flush=True)
        ##### END OUTPUT 1 #####
        
        self.logger.info('Finished testing new 2D simplification of the data')
        
        gene_list = list(sorted(set(d[0] for d in data)))
        locus_list = list(sorted(set(d[1] for d in data)))
        locus_set = set(d[1] for d in data)
        genome_list = list(sorted(set(d[2] for d in data)))
        
        for d in data:
            sequence = genome_contigs_list[d[2]][d[4]]
            if ((d[3] == FAR_UPSTREAM) or (d[3] == FAR_DOWNSTREAM)):
                Primer.scan(sequence, gene=d[0], locus=d[1], genome=d[2], region=d[3], contig=d[4], orientation=d[5], start=d[6], end=d[7], name=None, primer_size=(19,36), case=args.case, min_junction_overlap=None)
            else:
                Primer.scan(sequence, gene=d[0], locus=d[1], genome=d[2], region=d[3], contig=d[4], orientation=d[5], start=d[6], end=d[7], name=None, primer_size=(19,36), case=args.case)
        self.logger.info("Total 'Primer' objects before filtering: {}".format(len(Primer.sequences)))
        
        
        # Load cached Primer objects if instructed to
        if args.cache:
            cache_dir = os.path.join(args.folder, 'cache')
            object_name = 'Primer.sequences.pickle'
            if os.path.exists(os.path.join(cache_dir, object_name)):
                self.logger.info("Loading 'Primer' objects from 'cache' folder")
                cache_data = cache.load_object(object_name, cache_dir)
                
                # This simple nested loop is the fastest way to do this
                for k in Primer.sequences:
                    if k in cache_data:
                        Primer.sequences[k] = cache_data[k]
        
        # Build the required pattern
        #                      A    B    C    D    B    C    D    B    C    D
        # required_pattern = ['y', 'y', 'y', 'y', 'y', 'y', 'n', 'y', 'y', 'y']
        if (args.o_primers_required == None):
            o_primers_required = ['n']*len(genome_list)
        else:
            o_primers_required = args.o_primers_required
        
        if (args.i_primers_required == None):
            i_primers_required = ['n']*len(genome_list)
        else:
            i_primers_required = args.i_primers_required
            
        required_pattern = ['y'] # for A
        for oo, ii in zip(o_primers_required, i_primers_required):
            required_pattern.append(oo) # For B
            required_pattern.append(oo) # For C
            required_pattern.append(ii) # For D
        
        # Make simple dict to lookup the feature size given the gene/locus/contig
        gg2feature_size = {}
        for d in data:
            d_gene = d[0]
            d_locus = d[1]
            d_genome = d[2]
            d_region = d[3]
            d_contig = d[4]
            d_start, d_end = d[6], d[7]
            
            if (d_region == FEATURE_F):
                # This size is too small:
                #   For instance, ADE2 should be 1707 nt, but it is reporting 1706.
                gg2feature_size.setdefault((d_gene, d_locus, d_genome, d_contig), list()).append(d_end-d_start)
        
        self.logger.info("gene-to-feature:")
        for k, v in gg2feature_size.items():
            self.logger.info("  {} {}".format(k, v))
        
        for gene in gene_list:
            self.logger.info("Processing gene: {}".format(gene))
            # Make list of SHARED sF primers, and perform progressive checks on them
            sF_seq_list = []
            sR_seq_list = []
            featureF_seq_lists = [list() for r in genome_list]
            featureR_seq_lists = [list() for r in genome_list]
            
            # Get the reference locations needed for 'sF' and 'sR'
            sF_loc_set = set()
            sR_loc_set = set()
            for d in data:
                if (d[0] == gene):
                    if (d[3] == FAR_UPSTREAM):
                        sF_loc_set.add(tuple(d[:-2]))
                    elif (d[3] == FAR_DOWNSTREAM):
                        sR_loc_set.add(tuple(d[:-2]))
            
            # Get the reference locations needed for 'oF', 'oR', 'iF', and 'iR'
            fF_loc_sets = [set() for r in genome_list]
            fR_loc_sets = [set() for r in genome_list]
            for d in data:
                if (d[0] == gene):
                    if (d[3] == FEATURE_F):
                        fF_loc_sets[d[2]].add(tuple(d[:-2]))
                    if (d[3] == FEATURE_R):
                        fR_loc_sets[d[2]].add(tuple(d[:-2]))
            
            for pi, (seq, p) in enumerate(Primer.sequences.items()):
                # We strip out the start/end
                #   location = (Gene, Locus, Genome, Region, Contig, Strand, Start, End)
                p_loc_set = set(loc[:-2] for loc in p.locations) # Contig doesn't matter for FAR_UPSTREAM or FAR_DOWNSTREAM
                
                
                ###### BEGIN THIS ######
                # These 'intersection' expressions is where the allele-specific calculations should be made
                #if (args.specificity == 'exclusive'): # allele-specific
                #elif (args.specificity == 'all'): # multi-allelic
                #elif (args.specificity == 'any'): # allele-agnostic
                
                # If all required positions are present, then we add it
                if (sF_loc_set.intersection(p_loc_set) == sF_loc_set):
                    sF_seq_list.append(seq)
                
                if (sR_loc_set.intersection(p_loc_set) == sR_loc_set):
                    sR_seq_list.append(seq)
                
                # If all required regions are present, then we add it
                for r in range(len(genome_list)):
                    if (fF_loc_sets[r].intersection(p_loc_set) == fF_loc_sets[r]):
                        featureF_seq_lists[r].append(seq)
                    
                    if (fR_loc_sets[r].intersection(p_loc_set) == fR_loc_sets[r]):
                        featureR_seq_lists[r].append(seq)
                ###### END THIS ######
            
            self.logger.info("Total 'sF' 'Primer' sequences after filtering: {}".format(len(sF_seq_list)))
            self.logger.info("Total 'sR' 'Primer' sequences after filtering: {}".format(len(sR_seq_list)))
            for r in range(len(genome_list)):
                self.logger.info("Total 'iF/oF-r{}' 'Primer' sequences after filtering: {}".format(genome_list[r], len(featureF_seq_lists[r])))
                self.logger.info("Total 'iR/oR-r{}' 'Primer' sequences after filtering: {}".format(genome_list[r], len(featureR_seq_lists[r])))
            
            # Define cutoff start and ends (default)
            clist = [
                #   Name, initial, final, delta # states
                Cutoff("length", (19,28), (19,36), (0,2), separate=True), # 1,5
                Cutoff("last5gc_count", (1,3), (0,4), (-1,1), separate=True), # 2, last5gc
                Cutoff("gc_clamp_length", (1,2), (0,4), (-1,1), separate=True), # 2, gcclamp
                Cutoff("gc", (0.4, 0.6), (0.2, 0.8), (-0.1, 0.1)), # 3, gcfreq
                Cutoff("max_run_length", (4,), (6,), (1,)), # 3, runlen_max
                Cutoff("max_3prime_complementation_length", (3,), (6,), (1,)), # 4, 3primecomplen_max
                Cutoff("min_delta_g", (-4.0,), (-7.0,), (-1.0,)), # 4, deltag_min
                Cutoff("tm", (52,65), (52,65), (0,0)), # 1, tm
                Cutoff("max_tm_difference", (2.5,), (4.0,), (0.5,)), # 4, deltatm_max
                Cutoff("amplicon_size_range", (300,700), (300,700), (0,0)), # TODO: Make AddTag intellegently choose this based on the Feature and Insert sizes 
            ]
            
            # Create Iterator that returns the next cutoff and increments appropriately when 'next()' is called
            citer = CutoffIterator()
            
            design_found = False
            # Do calculations until either:
            #  * all primers are assayed,
            #  * an adequate design is found,
            #  * or the time limit expires
            cycle_n = 0
            
            while (cycle_n < args.cycle_start):
                self.logger.info("gene: {}, cycle: {}".format(gene, cycle_n))
                try:
                    self.logger.info("  Calculating next set of cutoffs")
                    # Get the initial cutoffs (if this is the first loop)
                    # or get the next-most relaxed cutoffs (if this isn't the first loop)
                    cutoffs = next(citer)
                    cycle_n += 1
                except StopIteration:
                    self.logger.info("  No additional cutoffs. Ending loop")
                    break
            
            optimal_designs = []
            while ((design_found == False) and (cycle_n <= args.cycle_stop)):
                
                self.logger.info("gene: {}, cycle: {}".format(gene, cycle_n))
                try:
                    self.logger.info("  Calculating next set of cutoffs")
                    # Get the initial cutoffs (if this is the first loop)
                    # or get the next-most relaxed cutoffs (if this isn't the first loop)
                    cutoffs = next(citer)
                    cycle_n += 1
                except StopIteration:
                    self.logger.info("  No additional cutoffs. Ending loop")
                    break
                
                # Add the invariant cutoff parameters
                cutoffs['o_oligo'] = args.selected_oligo
                cutoffs['folder'] = os.path.join(args.temp_folder, 'addtag', os.path.basename(args.folder))
                
                
                self.logger.info("  cutoffs = {}".format(cutoffs))
                
                # Do 'Primer' calculations until the time limit is reached
                p_queue = [sF_seq_list, sR_seq_list] + featureF_seq_lists + featureR_seq_lists
                #p_queue_labels = ['sF', 'sR'] + ['r'+str(r) for r in genome_list] + ['r'+str(r) for r in genome_list]
                
                
                ###### Start multiprocessing ######
                
                for seq_list in p_queue:
                    primers_to_process_list = [Primer.sequences[s] for s in seq_list]
                    self.mp_setup(args, cutoffs, primers_to_process_list)
                
                # Calculate Tm for Primers if missing
                #batch_tm_fxn = getattr(args.selected_oligo, 'find_tms', None)
                #if callable(batch_tm_fxn):
                
                if 'find_tms' in detect_overridden_methods(Oligo, args.selected_oligo):
                    # Make list of sequences whose Tm should be calculated
                    s_list = []
                    for s, p in Primer.sequences.items():
                        # If the RevComp structure was found
                        if p.o_reverse_complement:
                            # If the Tm has not been calculated
                            if (min(p.o_reverse_complement).melting_temperature == None):
                                s_list.append(s)
                    
                    # Find Tms as a batch operation
                    tm_list = args.selected_oligo.find_tms(s_list)
                    
                    # Apply the Tms to each RevComp Structure
                    for sl, tl in zip(s_list, tm_list):
                        for x in Primer.sequences[sl].o_reverse_complement:
                            x.melting_temperature = tl
                            # TODO: The setting of 'x.melting_temperature' should be done within the 'Primer'/'Structure' object, and not here
                
                ###### End multiprocessing ######
                
#                for seq_list in p_queue:
#                    # For each region, we reset the timer
#                    start_time = time.time()
#                    time_expired = False
#                    
#                    # reset metrics for this region
#                    p_tot = 0
#                    p_checked_now = 0
#                    p_passed_previously = 0
#                    p_rejected = 0
#                    p_skipped = 0
#                    p_newly_passed = 0
#                    
#                    for pi, seq in enumerate(seq_list):
#                        p = Primer.sequences[seq]
#                    #for pi, (seq, p) in enumerate(Primer.sequences.items()):
#                        if (pi % 1000 == 0):
#                            if (args.primer_scan_limit and ((time.time()-start_time) > args.primer_scan_limit)):
#                                time_expired = True
#                        
#                        if not time_expired:
#                            cpass = p.summarize(p.checks)
#                            if ((p.checks[0] == None) or (not cpass)):
#                                p.progressive_check(cutoffs)
#                                p_checked_now += 1
#                                if ((p.checks[0] != None) and p.summarize(p.checks)):
#                                    p_newly_passed += 1
#                            elif cpass:
#                                p_passed_previously += 1
#                            else:
#                                p_rejected += 1
#                        else:
#                            p_skipped += 1
#                        
#                        p_tot += 1
#                
#                    self.logger.info("  'Primer' objects: checked_now={}, passed_previously={}, newly_passed={}, not_checked={}, skipped={}, total={}".format(p_checked_now, p_passed_previously, p_newly_passed, p_rejected, p_skipped, p_tot))
                
                
                ##### Some debug code #####
                sss = 0
                ttt = 0
                nnn = 0
                hhh = 0
                ddd = 0
                ooo = 0
                for seq, p in Primer.sequences.items():
                    ttt += 1
                    if p.summarize(p.checks):
                        sss += 1
                    if ((p.checks[0] != None) and p.summarize(p.checks)):
                        nnn += 1
                    if (p.o_hairpin != None):
                        hhh += 1
                    if (p.o_self_dimer != None):
                        ddd += 1
                    if (p.o_reverse_complement != None):
                        ooo += 1
                self.logger.info("  Summary of all 'Primer' objects:")
                self.logger.info("    Number primers with summarize() passed = {}".format(sss))
                self.logger.info("    Number primers with all checks passed = {}".format(nnn))
                self.logger.info("    Number primers in 'Primer.sequences' = {}".format(ttt))
                self.logger.info("    Number primers with hairpins calculated = {}".format(hhh))
                self.logger.info("    Number primers with self-dimers calculated = {}".format(ddd))
                self.logger.info("    Number primers with Tm calculated = {}".format(ooo))
                ##### Some debug code #####
                
                # Save cached Primer objects if instructed to
                if args.cache:
                    self.logger.info("Saving 'Primer' objects to 'cache' folder")
                    cache_dir = os.path.join(args.folder, 'cache')
                    object_name = 'Primer.sequences.pickle'
                    cache.save_object(Primer.sequences, object_name, cache_dir)
                
                # Set up 'PrimerPair' objects, and place them in a queue 'PrimerPair.pairs'
                self.logger.info("  Queueing 'PrimerPair' objects")
                
                #pair_queue = []
                # pair_queue = [
                #     [                        # Locus 0
                #         [PrimerPair(), ...], # ((FAR_UPSTREAM, '+', 'sF'), (FAR_DOWNSTREAM, '-', 'sR')),
                #         [PrimerPair(), ...], # ((FAR_UPSTREAM, '+', 'sF'), (FEATURE, '-', 'oR'))
                #         ...
                #     ],
                #     ...
                # ]
                
                # Limit number of Primer pairs
                subset_size = args.subset_size # 1000 # 10000x10000 takes ~3.5Gb ram and 1.5 hours to calculate
                
                
                
                # Pair the lists, and include their names and amplicon label
                pairs = [
                    [sF_seq_list, sR_seq_list, 'sF', 'sR', 'A']
                ]
                for r in range(len(genome_list)):
                    pairs.append([sF_seq_list, featureR_seq_lists[r], 'sF', 'r{}-oR'.format(genome_list[r]), 'B'])
                    pairs.append([featureF_seq_lists[r], sR_seq_list, 'r{}-oF'.format(genome_list[r]), 'sR', 'C'])
                    pairs.append([featureF_seq_lists[r], featureR_seq_lists[r], 'r{}-iF'.format(genome_list[r]), 'r{}-iR'.format(genome_list[r]), 'D'])
                
                
                
                
                
                
                # Make list of viable 'PrimerPair' objects
                pp_queue = []
                should_skip = False
                for pi, (seq1_list, seq2_list, lab1, lab2, amp) in enumerate(pairs):
                    p1_list = []
                    p2_list = []
                    
                    # If 'Primer' has passed all checks, then include it
                    for seq in seq1_list:
                        p = Primer.sequences[seq]
                        if ((p.checks[0] != None) and Primer.summarize(p.checks)):
                            p1_list.append(p)
                    for seq in seq2_list:
                        p = Primer.sequences[seq]
                        if ((p.checks[0] != None) and Primer.summarize(p.checks)):
                            p2_list.append(p)
                    
                    # For each pair, we reset the timer
                    start_time = time.time()
                    time_expired = False
                    
                    
                    # The p1 x p2 nested for loop takes a reasonable amount of time
                    # if we limit the length of each list to 10,000 (see below)
                    #   ##### BEGIN SAMPLE #####
                    #   t = 0
                    #   for i in range(10000):
                    #      for j in range(10000):
                    #          t += i*j
                    #   ##### END SAMPLE #####
                    self.logger.info("    Number potential 'PrimerPair' objects for pair: {} = {} x {}".format(len(p1_list) * len(p2_list), len(p1_list), len(p2_list)))
                    pp_seq_list = []
                    primerpairs_to_process_list = []
                    for i1, p1 in enumerate(sorted(p1_list, reverse=True)[:subset_size]):
                        if (args.primer_pair_limit and ((time.time()-start_time) > args.primer_pair_limit)):
                            time_expired = True
                        
                        if (i1 % 100 == 0):
                            self.logger.info("      Processed/Queued {} pairs".format(i1*len(p2_list)))
                        
                        if not time_expired:
                            for i2, p2 in enumerate(sorted(p2_list, reverse=True)[:subset_size]):
                                # If it doesn't exist, add PrimerPair to database.
                                # Otherwise, do nothing.
                                PrimerPair(p1, p2)
                                
                                # Get the 'PrimerPair' object from the dict
                                pair = (p1.sequence, p2.sequence)
                                pp = PrimerPair.pairs.get(pair)
                                
                                if pp:
                                    primerpairs_to_process_list.append(pp)
                                    
#                                    # Run checks
#                                    pp.progressive_check(cutoffs)
#                                    
#                                    # If the 'PrimerPair' passes the checks, then add it
#                                    if ((pp.checks[0] != None) and Primer.summarize(pp.checks)):
#                                        pp_seq_list.append(pair)
                    
                    # Load cached PrimerPair objects if instructed to
                    if args.cache:
                        cache_dir = os.path.join(args.folder, 'cache')
                        object_name = 'PrimerPair.pairs.pickle'
                        if os.path.exists(os.path.join(cache_dir, object_name)):
                            self.logger.info("Loading 'PrimerPair' objects from 'cache' folder")
                            cache_data = cache.load_object(object_name, cache_dir)
                            
                            # This simple nested loop is the fastest way to do this
                            for k in PrimerPair.pairs:
                                if k in cache_data:
                                    PrimerPair.pairs[k] = cache_data[k]
                    
                    ###### Start multiprocessing ######
                    
                    self.mpp_setup(args, cutoffs, primerpairs_to_process_list)
                    
                    ###### End multiprocessing ######
                    
                    for i1, p1 in enumerate(sorted(p1_list, reverse=True)[:subset_size]):
                        for i2, p2 in enumerate(sorted(p2_list, reverse=True)[:subset_size]):
                            pair = (p1.sequence, p2.sequence)
                            pp = PrimerPair.pairs.get(pair)
                            
                            if pp:
                                # If the 'PrimerPair' passes the checks, then add it
                                if ((pp.checks[0] != None) and Primer.summarize(pp.checks)):
                                    pp_seq_list.append(pair)
                        
                    pp_queue.append(pp_seq_list)
                    self.logger.info("      Added {} 'PrimerPair' objects for pair: '{}', '{}', amp={}, required={}".format(len(pp_seq_list), lab1, lab2, amp, required_pattern[pi]))
                    
                    # TODO: Each pp should have its own cutoffs!
                    # Ideally:
                    #   Each pp should have its own cutoffs:
                    #     A, r0:B/C/D, r1:B/C/D, r2:B/C/D
                    #   Thus, 10 cutoffs.
                    #   If this specific set of PrimerPairs has 0 passing,
                    #     then cutoffs should be relaxed.
                    #   Thus, each PrimerPair group would have its own set of cutoffs.
                    #   This would minimize relaxing cutoffs that don't need it.
                    if ((required_pattern[pi] in ['y', 'Y', '1', 'T', 't', 'TRUE', 'True', 'true']) and (len(pp_seq_list) == 0)):
                        should_skip = True
                
                # Continue performing the pp calculations for the current round,
                # Then skip without doing any simulated annealing
                if should_skip: # TODO: should this also be 'args.skip_partial_sets'???
                    self.logger.info("  Skipping simulated annealing")
                    continue
                
                # Now we need to take the 'pp_queue' and slot it into a 'PrimerDesign', and optimize
                # ...do calculations here...
                
                # Populate sF_sR_paired_primers with primer objects
                sF_sR_paired_primers = []
                
                #self.logger.info("    Gene={}, number_features={}, Features={}".format(gene, ))
                
                for pair in pp_queue[0]:
                    pp = PrimerPair.pairs[pair]
                    #pp.weight = pp.get_weight(locus, FAR_UPSTREAM, FAR_DOWNSTREAM, contig, minimize=True)
                    #pp.weight = pp.get_weight(minimize=True)
                    pp.weight = pp.get_weight(minimize=gg2feature_size)
                    sF_sR_paired_primers.append(pp)
                
                # Save cached PrimerPair objects if instructed to
                if args.cache:
                    self.logger.info("Saving 'PrimerPair' objects to 'cache' folder")
                    cache_dir = os.path.join(args.folder, 'cache')
                    object_name = 'PrimerPair.pairs.pickle'
                    cache.save_object(PrimerPair.pairs, object_name, cache_dir)
                
                # Go through the 'pp_queue' one sF sR pair at a time
                design_count = 0
                mpd_list = []
                for loop_i, sF_sR_pair in enumerate(sorted(sF_sR_paired_primers, reverse=True)):
                    self.logger.info('  loop {}:'.format(loop_i))
                    pp_sources = []
                    pp_sources.append([sF_sR_pair])
                    
                    # Order in 'pairs' variable
                    # all  r0       r1      r2
                    # A,   B,C,D,   B,C,D,  B,C,D
                    #      0 1 2    3 4 5   6 7 8
                    
                    #for ppi, pp_list in pp2d_list:
                    for ppi, pp_seq_list in enumerate(pp_queue[1:]):
                        pp_list = []
                        for pair in pp_seq_list:
                            pp_list.append(PrimerPair.pairs[pair])
                        
                        if (ppi % 3 == 0): # B
                            pp_sources.append(self.filter_primer_pairs(pp_list, forward=sF_sR_pair.forward_primer))
                        elif (ppi % 3 == 1): # C
                            pp_sources.append(self.filter_primer_pairs(pp_list, reverse=sF_sR_pair.reverse_primer))
                        elif (ppi % 3 == 2): # D
                            pp_sources.append(pp_list)
                    
                    current_pp_sources = [len(x) for x in pp_sources]
                    self.logger.info('  length of sources: {}'.format(current_pp_sources))
                    
                    
                    # Ideally, similar code should be executed, but before the 'PrimerPair.progressive_check()' is calculated
                    should_mask = False
                    should_skip = False
                    for req_str, num_pp in zip(required_pattern, current_pp_sources):
                        if ((req_str in ['y', 'Y', '1', 'T', 't', 'TRUE', 'True', 'true']) and (num_pp == 0)):
                            should_skip = True
                            break
                        if ((req_str not in ['y', 'Y', '1', 'T', 't', 'TRUE', 'True', 'true']) and (num_pp > 0)):
                            should_mask = True
                    if args.skip_partial_sets and should_skip:
                        self.logger.info('  skipping...')
                        continue
                    if should_mask:
                        self.logger.info('  masking...') # TODO: Re-implement masking (see 'calculate_them_primers()' and 'calculate_them_best_set()'
                        continue
                    
                    # Evaluate if any new primer designs are adequate
                    # This assumes that pp_sources has valid primer pair lists as each element
                    # And thus that an optimal solution is possible
                    design = PrimerDesign(pp_sources)
                    mpd_list.append(design)
#                    design.optimize(mode='direct')
#                    #design.optimize(iterations=3000)
#                    optimal = design.optimal
#                    
#                    if optimal:
#                        optimal_designs.append(optimal)
#                        design_found = True
                    
                    design_count += 1
                    
                    #if (design_count >= subset_size):
                    #    self.logger.info("  Stopped queueing designs for 'optimization' because 'subset_size' has been reached.")
                    #    break
                    if (design_count >= args.max_number_designs_reported):
                        self.logger.info("  Stopped queueing designs for 'optimization' because 'max_number_designs_reported' has been reached.")
                        break
                
                ###### Start multiprocessing ######
                
                mpd_list = self.mpd_setup(args, mpd_list)
                for optimal_set in mpd_list:
                    if optimal_set:
                        optimal_designs.append(optimal_set)
                        design_found = True
                
                ###### End multiprocessing ######
                
                
            
            # Write output for Primer sequences
            if (len(optimal_designs) > 0):
                # Sort first by number of PrimerPair objects, and second by the weight
                ordered_od_list = sorted(optimal_designs, key=lambda x: (x.get_primer_pair_count(), x.weight), reverse=True)
                
                print('# Compatible designs')
                print('# Primer sequences')
                print('# ' + '\t'.join(['Gene', 'Weight', 'Count'] + [x[0] for x in optimal_designs[0].get_primer_list()]))
                for od in ordered_od_list: # od is a PrimerSet object
                    gpl = od.get_primer_list()
                    print('\t'.join([gene, str(od.weight), str(od.get_primer_count())] + [x[1].sequence if x[1] else 'None' for x in gpl]))
                
                print('# PrimerPairs')
                print('# ' + '\t'.join(['', '', 'Amplicon'] + [x[0] for x in optimal_designs[0].get_primer_pair_list()]))
                print('# ' + '\t'.join(['Gene', 'Weight', 'Count'] + [x[1] for x in optimal_designs[0].get_primer_pair_list()]))
                for od in ordered_od_list: # od is a PrimerSet object
                    gppl = od.get_primer_pair_list()
                    print('\t'.join([gene, str(od.weight), str(od.get_primer_pair_count())] + [str(x[2]) for x in gppl]))
            else:
                print('# No compatible designs')
            
        self.logger.info("Function 'primer_queue_test()' completed.")
    
    def mp_setup(self, args, cutoffs, data):
        '''
        Use multiprocessing to run 'progressive_checks()' on each 'Primer' object.
        '''
        
        # Create queue that will hold 'Primer' objects
        mp_queue = multiprocessing.JoinableQueue()
        results_queue = multiprocessing.Queue()
        
        # Create locks
        log_lock = multiprocessing.Lock()
        data_lock = multiprocessing.Lock()
        
        # Create a number of processes equal to input cores
        self.logger.info('Creating mp jobs...')
        jobs = []
        for i in range(args.processors):
            w = PrimerWorker(args, cutoffs, mp_queue, results_queue, log_lock, data_lock)
            jobs.append(w)
            w.start()
        
        # After the jobs are created, we populate the queue
        self.logger.info('Adding data to qeueue...')
        for d in data:
            mp_queue.put(d)
        
        # Add a poison pill for each processor
        self.logger.info('Adding kill signals to qeueue...')
        for i in range(args.processors):
            mp_queue.put(None)
        
        # exitcode - The child’s exit code will be None if the process has not yet terminated.
        # A negative value -N indicates that the child was terminated by signal N.
        self.logger.info('Waiting for threads to finish...')
        while any([j.exitcode == None for j in jobs]):
        #while any([j.is_alive() for j in jobs]):
            try:
                seq, results = results_queue.get(timeout=1) # May need to use try/except block, or include a timeout (in seconds)
                p = Primer.sequences[seq]
                p.mp_update(results)
            except queue.Empty:
                self.logger.info('queue empty')
        
        # Wait for all processes to end
        for j in jobs:
            j.join()
        
        self.logger.info('mp finished')
    
    def mpp_setup(self, args, cutoffs, data):
        '''
        Use multiprocessing to run 'progressive_checks()' on each 'PrimerPair' object.
        '''
        # Each process can do computations on a single entry to 'pairs'
        
        # Create queue that will hold 'Primer' objects
        mp_queue = multiprocessing.JoinableQueue()
        results_queue = multiprocessing.Queue()
        
        # Create locks
        log_lock = multiprocessing.Lock()
        data_lock = multiprocessing.Lock()
        
        # Create a number of processes equal to input cores
        self.logger.info('Creating mpp jobs...')
        jobs = []
        for i in range(args.processors):
            w = PrimerPairWorker(args, cutoffs, mp_queue, results_queue, log_lock, data_lock)
            jobs.append(w)
            w.start()
        
        # After the jobs are created, we populate the queue
        self.logger.info('Adding data to qeueue...')
        for d in data:
            mp_queue.put(d)
        
        # Add a poison pill for each processor
        self.logger.info('Adding kill signals to qeueue...')
        for i in range(args.processors):
            mp_queue.put(None)
        
        # exitcode - The child’s exit code will be None if the process has not yet terminated.
        # A negative value -N indicates that the child was terminated by signal N.
        self.logger.info('Waiting for threads to finish...')
        while any([j.exitcode == None for j in jobs]):
        #while any([j.is_alive() for j in jobs]):
            try:
                pair, results = results_queue.get(timeout=1) # May need to use try/except block, or include a timeout (in seconds)
                pp = PrimerPair.pairs.get(pair)
                if pp:
                    pp.mp_update(results)
            except queue.Empty:
                self.logger.info('queue empty')
        
        # Wait for all processes to end
        for j in jobs:
            j.join()
        
        self.logger.info('mpp finished')
    
    def mpd_setup(self, args, data):
        '''
        Use multiprocessing to run 'oligo.PrimerDesign.optimize()'.
        '''
        
        # Create queue that will hold 'Primer' objects
        mp_queue = multiprocessing.JoinableQueue()
        results_queue = multiprocessing.Queue()
        
        # Create locks
        log_lock = multiprocessing.Lock()
        data_lock = multiprocessing.Lock()
        
        # Create a number of processes equal to input cores
        self.logger.info('Creating mpd jobs...')
        jobs = []
        for i in range(args.processors):
            w = PrimerDesignWorker(args, mp_queue, results_queue, log_lock, data_lock)
            jobs.append(w)
            w.start()
        
        # After the jobs are created, we populate the queue
        self.logger.info('Adding data to qeueue...')
        for d in data:
            mp_queue.put(d)
        
        # Add a poison pill for each processor
        self.logger.info('Adding kill signals to qeueue...')
        for i in range(args.processors):
            mp_queue.put(None)
        
        # exitcode - The child’s exit code will be None if the process has not yet terminated.
        # A negative value -N indicates that the child was terminated by signal N.
        self.logger.info('Waiting for threads to finish...')
        data = []
        while any([j.exitcode == None for j in jobs]):
        #while any([j.is_alive() for j in jobs]):
            try:
                d = results_queue.get(timeout=1) # May need to use try/except block, or include a timeout (in seconds)
                data.append(d)
            except queue.Empty:
                self.logger.info('queue empty')
        
        # Wait for all processes to end
        for j in jobs:
            j.join()
        
        self.logger.info('mpd finished')
        
        return data
    
    def make_datum_groups(self, group_links, contig_groups, genome_contigs_list):
        ######## Begin for linking the loci ########
        
        self.logger.info('contig_groups:')
        for cg in contig_groups:
            self.logger.info('  ' + str(cg))
        
        def in_same_contig_group(contig1, contig2, cgroups=contig_groups):
            for g in cgroups:
                if ((contig1 in g) and (contig2 in g)):
                    return True
            else:
                return False
        
        
        
        # Make a copy of 'group_links' that we can pop stuff from
        datum_list = copy.deepcopy(group_links)
        
        self.logger.info('datum_list:')
        for datum in datum_list:
            self.logger.info('  ' + str(datum))
        
        
        # For each (round, locus), we need the (contig name, upstream_homology_length, downstream_homology_length)
        # loop through these:
        
        # For each cross-over event, there will be a datum
        # Typically, there should be 1 Datum object for haploid genomes,
        # and 2 Datum objects for diploid genomes
        
        
        # Populate initial groups based on if there is crossing over in the genome-r0
        datum_groups = []
        to_pop = []
        for di, datum in enumerate(datum_list):
            if datum.genome_contig in genome_contigs_list[0].keys(): # The 0th index refers to genome-r0
                datum_groups.append([datum])
                to_pop.append(di)
        for di in sorted(to_pop, reverse=True):
            datum_list.pop(di)
        to_pop = []
        
        self.logger.info('datum_groups (r0):')
        for dg in datum_groups:
            self.logger.info('  ' + str(dg))
        
        self.logger.info('datum_list:')
        for datum in datum_list:
            self.logger.info('  ' + str(datum))
        
        # Ideally, linking should be done by genomic position (after the blast/alignment results)
        
        # Populate the intermediary groups
        # New try:
        # If there is ANY OVERLAP AT ALL between
        #      round n,    AFTER: ush_start-dsh_end
        # and  round n+1, BEFORE: ush_start-dsh_end
        # Then these loci should be linked!
        no_more_after = False
        no_more_before = False
        groups = []
        wcount = 0
        while (not no_more_after):
            # TODO: Is using 'wcount' actually necessary? can't we just use 'len(genome_contigs_list)'?
            # Specify an error if for some reason this doesn't work
            wcount += 1
            if (wcount > 10000):
                #raise Exception('Too many iterations taken to link loci between engineered genomes.')
                break # Leave the while loop (this is triggered when there is only the 2 genomes in 'genome_contigs_list'?
            
            for target_genome_round in range(1, len(genome_contigs_list)-1):
                to_pop = []
                # Select the first datum ("after")
                datum_after = None
                for di, datum in enumerate(datum_list):
                    if (datum.genome_r == target_genome_round):
                        # If the rounds are the same then the alignment is "after engineering"
                        if (datum.dDNA_r == datum.genome_r):
                            datum_after = datum
                            to_pop.append(di)
                            break
                else:
                    no_more_after = True
                
                # Select the second datum ("before")
                datum_before = None
                if (datum_after != None):
                    for di, datum in enumerate(datum_list):
                        if (datum.genome_r == target_genome_round):
                            # If the rounds are different then the alignment is "before engineering"
                            if (datum.dDNA_r != datum.genome_r):
                                if ((datum.genome_contig == datum_after.genome_contig) and
                                    (datum.dDNA_contig != datum_after.dDNA_contig)):
                                    datum_before = datum
                                    to_pop.append(di)
                                    break
                    else:
                        no_more_before = True
                
                else: # This means (datum_after == None)
                    # If there is no "before" to pair with selected "after"
                    # then this is a terminal genome engineering
                    # Not sure if this code would work if a locus was "skipped" between rounds as follows:
                    #    genome-r0: locus A wt
                    #    genome-r1: locus A ko
                    #    genome-r2: locus A ko (no change)
                    #    genome-r3: locus A ki
                    
                    # Since this is a terminal, then we need to select the correct "before" using homology
                    self.logger.info("DEBUG: This means (datum_after == None)")
                
                # Write some debug log messages
                self.logger.info("DEBUG:")
                self.logger.info("     datum_after={}".format(datum_after))
                self.logger.info("    datum_before={}".format(datum_before))
                self.logger.info("   no_more_after={}".format(no_more_after))
                self.logger.info("  no_more_before={}".format(no_more_before))
                
                if ((datum_after != None) and (datum_before != None)):
                    # Now that the two datum are selected, determine whether or not they overlap
                    # Feature.overlap_coverage(start1, end1, start2, end2, index_base=0)
                    overlap = Feature.overlap_coverage(datum_after.ush_start, datum_after.dsh_end, datum_before.ush_start, datum_before.dsh_end)
                    
                    # If they overlap
                    if (overlap > 0):
                        # Then add them to the same group
                        group = [datum_after, datum_before]
                        groups.append(group)
                        
                        for di in sorted(to_pop, reverse=True):
                            datum_list.pop(di)
                    # Otherwise, they are non-overlapping engineering events on the same genome contig
                    else:
                        # Do nothing
                        self.logger.info('DEBUG: This means that (overlap <= 0)')
        
        self.logger.info('groups:')
        for gi, (datum_after, datum_before) in enumerate(groups):
            self.logger.info('  i={}'.format(gi))
            self.logger.info('      datum_after={}'.format(datum_after))
            self.logger.info('     datum_before={}'.format(datum_before))
        
        self.logger.info('datum_list:')
        for datum in datum_list:
            self.logger.info('  {}'.format(datum))
        
        # Assuming we have 'group' with the [datum_after, datum_before]
        # Now we match homology regions with r0
        for dgi, dg in enumerate(datum_groups):
            # group = [
            #     Datum(dDNA_r=1, dDNA_contig='exDonor-55', genome_r=1, genome_contig='Ca22chr3A_C_albicans_SC5314-r1[exDonor-55]', ush_start=943366, ush_end=None, dsh_start=None, dsh_end=943468)
            #     Datum(dDNA_r=2, dDNA_contig='reDonor-0', genome_r=1, genome_contig='Ca22chr3A_C_albicans_SC5314-r1[exDonor-55]', ush_start=942963, ush_end=943420, dsh_start=943414, dsh_end=943996)
            # ]
            
            last_datum = dg[-1]
            ush_len = abs(last_datum.ush_end - last_datum.ush_start)
            dsh_len = abs(last_datum.dsh_end - last_datum.dsh_start)
            
            dg_ush_seq = genome_contigs_list[last_datum.genome_r][last_datum.genome_contig][last_datum.ush_start:last_datum.ush_start+ush_len]
            dg_dsh_seq = genome_contigs_list[last_datum.genome_r][last_datum.genome_contig][last_datum.dsh_end-dsh_len:last_datum.dsh_end]
            
            for datum_after, datum_before in groups:
                # For 'datum_after', the fields 'datum_after.ush_end' and 'datum_after.dsh_start' both equal 'None'
                g_ush_seq = genome_contigs_list[datum_after.genome_r][datum_after.genome_contig][datum_after.ush_start:datum_after.ush_start+ush_len]
                g_dsh_seq = genome_contigs_list[datum_after.genome_r][datum_after.genome_contig][datum_after.dsh_end-dsh_len:datum_after.dsh_end]
                
                ##### DEBUG CODE #####
                self.logger.info("  dgi = {}".format(dgi))
                self.logger.info("    ush_len = {}, 0.9*ush_len = {}".format(ush_len, 0.9*ush_len))
                self.logger.info("    dsh_len = {}, 0.9*dsh_len = {}".format(dsh_len, 0.9*dsh_len))
                us_lcs = nucleotides.lcs(dg_ush_seq, g_ush_seq)
                ds_lcs = nucleotides.lcs(dg_dsh_seq, g_dsh_seq)
                self.logger.info("    nucleotides.lcs(dg_ush_seq, g_ush_seq) = {}".format(us_lcs))
                self.logger.info("    nucleotides.lcs(dg_dsh_seq, g_dsh_seq) = {}".format(ds_lcs))
                self.logger.info("    us_lcs > 0.9*ush_len = {}".format(us_lcs.size > 0.9*ush_len))
                self.logger.info("    ds_lcs > 0.9*dsh_len = {}".format(ds_lcs.size > 0.9*dsh_len))
                us_pident = nucleotides.pident(dg_ush_seq, g_ush_seq)
                ds_pident = nucleotides.pident(dg_dsh_seq, g_dsh_seq)
                self.logger.info("    us_pident = {}, us_pident > 0.9 = {}".format(us_pident, us_pident > 0.9))
                self.logger.info("    ds_pident = {}, ds_pident > 0.9 = {}".format(ds_pident, ds_pident > 0.9))
                ######################
                
                # If the rounds match, the dDNA_contigs match, and the genome_contigs match
                # and the homology regions match
                if ((last_datum.dDNA_r == datum_after.dDNA_r) and
                    (last_datum.dDNA_contig == datum_after.dDNA_contig) and
                    in_same_contig_group(last_datum.genome_contig, datum_after.genome_contig) and
                    (nucleotides.pident(dg_ush_seq, g_ush_seq) > 0.9) and
                    (nucleotides.pident(dg_dsh_seq, g_dsh_seq) > 0.9)):
                    #(nucleotides.lcs(dg_ush_seq, g_ush_seq).size > 0.9*ush_len) and # Old code
                    #(nucleotides.lcs(dg_dsh_seq, g_dsh_seq).size > 0.9*dsh_len)): # Old code
                    # Then populate their missing fields
                    #if (datum_after.ush_end == None):
                    #    datum_after.ush_end = group[-1].ush_start+ush_len
                    #if (datum_after.dsh_start == None):
                    #    datum_after.dsh_start = group[-1].dsh_end-dsh_len
                    
                    # Then add them to the same group
                    dg.append(datum_after)
                    dg.append(datum_before)
                    self.logger.info('group [datum_after, datum_before] appended to current datum_group.')
        
        self.logger.info('datum_groups:')
        for dg in datum_groups:
            self.logger.info('  {}'.format(dg))
        
        # Populate the final/terminal group 'rN'
        # Using homology! (like we did with r0 above)
        for dgi, dg in enumerate(datum_groups):
            # group = [
            #     Datum(dDNA_r=1, dDNA_contig='exDonor-55', genome_r=1, genome_contig='Ca22chr3A_C_albicans_SC5314-r1[exDonor-55]', ush_start=943366, ush_end=None, dsh_start=None, dsh_end=943468)
            #     Datum(dDNA_r=2, dDNA_contig='reDonor-0', genome_r=1, genome_contig='Ca22chr3A_C_albicans_SC5314-r1[exDonor-55]', ush_start=942963, ush_end=943420, dsh_start=943414, dsh_end=943996)
            # ]
            last_datum = dg[-1]
            ush_len = abs(last_datum.ush_end - last_datum.ush_start)
            dsh_len = abs(last_datum.dsh_end - last_datum.dsh_start)
            
            dg_ush_seq = genome_contigs_list[last_datum.genome_r][last_datum.genome_contig][last_datum.ush_start:last_datum.ush_start+ush_len]
            dg_dsh_seq = genome_contigs_list[last_datum.genome_r][last_datum.genome_contig][last_datum.dsh_end-dsh_len:last_datum.dsh_end]
            
            to_pop = []
            for di, datum in enumerate(datum_list):
                d_ush_seq = genome_contigs_list[datum.genome_r][datum.genome_contig][datum.ush_start:datum.ush_start+ush_len]
                d_dsh_seq = genome_contigs_list[datum.genome_r][datum.genome_contig][datum.dsh_end-dsh_len:datum.dsh_end]
                
                ##### DEBUG CODE #####
                self.logger.info("  dgi = {}".format(dgi))
                self.logger.info("    ush_len = {}, 0.9*ush_len = {}".format(ush_len, 0.9*ush_len))
                self.logger.info("    dsh_len = {}, 0.9*dsh_len = {}".format(dsh_len, 0.9*dsh_len))
                us_pident = nucleotides.pident(dg_ush_seq, d_ush_seq)
                ds_pident = nucleotides.pident(dg_dsh_seq, d_dsh_seq)
                self.logger.info("    us_pident = {}, us_pident > 0.9 = {}".format(us_pident, us_pident > 0.9))
                self.logger.info("    ds_pident = {}, ds_pident > 0.9 = {}".format(ds_pident, ds_pident > 0.9))
                ######################
                
                # If the rounds match, the dDNA_contigs match, and the genome_contigs match
                # and the homology regions match
                if ((last_datum.dDNA_r == datum.dDNA_r) and
                    (last_datum.dDNA_contig == datum.dDNA_contig) and
                    in_same_contig_group(last_datum.genome_contig, datum.genome_contig) and
                    (nucleotides.pident(dg_ush_seq, d_ush_seq) > 0.9) and
                    (nucleotides.pident(dg_dsh_seq, d_dsh_seq) > 0.9)):
                    #(nucleotides.lcs(dg_ush_seq, d_ush_seq).size > 0.9*ush_len) and # Old code
                    #(nucleotides.lcs(dg_dsh_seq, d_dsh_seq).size > 0.9*dsh_len)): # Old code
                    
                    # Then add them to the same group
                    dg.append(datum)
                    to_pop.append(di)
            
            for di in sorted(to_pop, reverse=True):
                datum_list.pop(di)
        
        # Please note:
        #   This "finished" 'datum_groups' could use some work.
        #   The data structure of these linked 'Datum' objects should allow for branching
        #   For instance:
        #     wt1(1,0)──ko1(1,1)──ko1(2,1)──ki1(2,2)
        #      └────────ko1(1,1)──ko1(2,1)──ki1(2,2)
        #   I forgot why, but it is important
        #    ------A-----A-----  contig has 2 locations that should be targeted
        #                        Should PCR be allele/paralog-specific?
        #                        Or should the PCR try to encompass both alleles/paralogs?
        self.logger.info('datum_groups (finished):')
        for dg in datum_groups:
            self.logger.info('  ' + str(dg))
        
        self.logger.info('datum_list (finished):')
        if (len(datum_list) == 0):
            self.logger.info('  EMPTY')
        else:
            for datum in datum_list:
                self.logger.info('  ' + str(datum))
        
        ######## End code for linking the loci ########
        
        return datum_groups
    
    def new_sf_sr_primers(self, args):
        """
        This function should replace 'get_far_lcs_regions()'
        """
        # If allele-specific, then do this for each locus.
        # If NOT allele-specific, then include all alleles in the same 'Datum' (will this even work?)
        
        # The way it should work:
        #   for flank in ['us', 'ds']:
        #     primers_found = defaultdict() # starts every key at value=0
        #     for insert
        #       # Scan region using a sliding window to find putative primers
        #       for primer_size in range(19, 32):
        #         primer = ...
        #         primers_found[primer] += 1
        #     omnipresent_primers = []
        #     for p, c in primers_found.items()
        #       if (c == len(inserts))
        #         omnipresent_primers.append(p)
        
        # A drawback--need a way to keep track of the primer's position...
        # Advantages over old LCS version: does not require long stretches of homology
        pass
    
    def get_far_lcs_regions(self, datum_groups, genome_contigs_list):
        
        # Let's find the longest common substring in the far_upstream and far_downstream regions of each
        # We also need the distance (in nt) between this LCS and the feature/insert
        # (for both the upstream and downstream)
        # (so we can do proper amplicon size calculations later)
        pcr_regions = []
        
        # pcr_region_positions[datum_group index][datum index] = [fus_start, fus_end, fds_start, fds_end]
        pcr_region_positions = []
        
        for dg in datum_groups:
            # Find the upstream LCS
            us_dist = 600
            us_done = False
            fus0 = None
            while (us_done == False):
                fus0 = genome_contigs_list[dg[0].genome_r][dg[0].genome_contig][max(0, dg[0].ush_start-us_dist):dg[0].ush_start]
                for di, datum in enumerate(dg[1:]):
                    if ((us_dist > dg[0].ush_start) or (us_dist > datum.ush_start)):
                        self.logger.info('too far')
                        us_done = True
                    fus = genome_contigs_list[datum.genome_r][datum.genome_contig][max(0, datum.ush_start-us_dist):datum.ush_start]
                    m = nucleotides.lcs(fus0, fus)
                    self.logger.info('us_dist = ' + str(us_dist))
                    self.logger.info('   fus0 = ' + fus0)
                    self.logger.info(('fus'+str(di+1)).rjust(7) +' = ' + fus)
                    self.logger.info('      m = ' + str(m))
                    if (m.size < 200):
                        us_dist += 100
                        break
                    else:
                        fus0 = fus0[m.a:m.a+m.size]
                else:
                    us_done = True
            
            self.logger.info('far_upstream_dist: ' + str(us_dist))
            self.logger.info('far_upstream_seq >= 200: ' + fus0)
            
            # Find the downstream LCS
            ds_dist = 600
            ds_done = False
            fds0 = None
            while (ds_done == False):
                fds0 = genome_contigs_list[dg[0].genome_r][dg[0].genome_contig][dg[0].dsh_end:dg[0].dsh_end+ds_dist]
                for di, datum in enumerate(dg[1:]):
                    if ((dg[0].dsh_end+ds_dist > len(genome_contigs_list[dg[0].genome_r][dg[0].genome_contig])) or (datum.dsh_end+ds_dist > len(genome_contigs_list[datum.genome_r][datum.genome_contig]))):
                        self.logger.info('too far')
                        ds_done = True
                    fds = genome_contigs_list[datum.genome_r][datum.genome_contig][datum.dsh_end:datum.dsh_end+ds_dist]
                    m = nucleotides.lcs(fds0, fds)
                    self.logger.info('us_dist = ' + str(ds_dist))
                    self.logger.info('   fus0 = ' + fds0)
                    self.logger.info(('fus'+str(di+1)).rjust(7) +' = ' + fds)
                    self.logger.info('      m = ' + str(m))
                    if (m.size < 200):
                        ds_dist += 100
                        break
                    else:
                        fds0 = fds0[m.a:m.a+m.size]
                else:
                    ds_done = True
            
            self.logger.info('far_downstream_dist: ' + str(ds_dist))
            self.logger.info('far_downstream_seq >= 200: ' + fds0)
            
            # Add PCR regions to data structure
            pcr_regions.append((fus0, fds0))
            
            # Need to calculate the position of 'sF' and 'sR' within each genome associated with the Datum
            pcr_region_positions.append([])
            for datum in dg:
                pcr_region_positions[-1].append([])
                
                # Add sF_start, sF_end
                temp_var = max(0, datum.ush_start-us_dist)
                fus_start = genome_contigs_list[datum.genome_r][datum.genome_contig][temp_var:datum.ush_start].rindex(fus0) + temp_var
                fus_end = fus_start + len(fus0)
                pcr_region_positions[-1][-1].append(fus_start)
                pcr_region_positions[-1][-1].append(fus_end)
                
                # Add sR_start, sR_end
                fds_start = genome_contigs_list[datum.genome_r][datum.genome_contig][datum.dsh_end:datum.dsh_end+ds_dist].index(fds0) + datum.dsh_end
                fds_end = fds_start + len(fds0)
                pcr_region_positions[-1][-1].append(fds_start)
                pcr_region_positions[-1][-1].append(fds_end)
            
            self.logger.info("checking 'sF' region:")
            self.logger.info("              fus0: " + fus0)
            for prpi2, prp2 in enumerate(pcr_region_positions[-1]):
                self.logger.info(str(prpi2) + ' ' + str(prp2[:2]).rjust(16) + ": " + genome_contigs_list[dg[prpi2].genome_r][dg[prpi2].genome_contig][prp2[0]:prp2[1]] + ' ' + str(dg[prpi2]))
            self.logger.info("checking 'sR' region:")
            self.logger.info("              fds0: " + fds0)
            for prpi2, prp2 in enumerate(pcr_region_positions[-1]):
                self.logger.info(str(prpi2) + ' ' + str(prp2[2:]).rjust(16) + ": " + genome_contigs_list[dg[prpi2].genome_r][dg[prpi2].genome_contig][prp2[2]:prp2[3]] + ' ' + str(dg[prpi2]))
        
        return pcr_regions, pcr_region_positions
    
    def calculate_them_primers(self, args, far_us_seq, far_ds_seq, insert_seqs):
        """
        """
        primer_length_range = (19,32)
        
        # Limit for the number of primers that go into the pair() function.
        # At most, there will be 1000*1000 = 1 millon possible pairs
        subset_size = 1000
        
        self.logger.info('Scanning the regions (shared upstream (sF), shared downstream (sR), feature/insert (rN-oF,rN-oR,rN-iF,rN-iR) for all decent primers')
        
        # Make the 'temp_folder' path
        temp_folder = os.path.join('/dev/shm/addtag', os.path.basename(args.folder))
        
        self.logger.info("Scanning far upstream for 'sF' primers:")
        sF_list = sorted(
            args.selected_oligo.scan(far_us_seq, 'left', primer_size=primer_length_range, folder=temp_folder, time_limit=args.primer_scan_limit),
            key=lambda x: x.weight,
            reverse=True
        )
        self.logger.info('  len(sF_list) = {}'.format(len(sF_list)))
        self.logger.info('  sF: skipping {}/{} calculated primers'.format(max(0, len(sF_list)-subset_size), len(sF_list)))
        sF_list = sF_list[:subset_size]
        # Need to add a condition that if (len(sF_list) < N), then it should re-evaluate the far-upstream region to try to get
        # more sequence to search.
        
        self.logger.info("Scanning far downstream for 'sR' primers:")
        sR_list = sorted(
            args.selected_oligo.scan(far_ds_seq, 'right', primer_size=primer_length_range, folder=temp_folder, time_limit=args.primer_scan_limit),
            key=lambda x: x.weight,
            reverse=True
        )
        self.logger.info('  len(sR_list) = {}'.format(len(sR_list)))
        self.logger.info('  sR: skipping {}/{} calculated primers'.format(max(0, len(sR_list)-subset_size), len(sR_list)))
        sR_list = sR_list[:subset_size]
        # Need to add condition that if there aren't enough putative 'sR' primers, then the sR region is expanded
        # Otherwise, the script can just end prematurely
        
        insert_list = []
        for ins in insert_seqs:
            self.logger.info("Scanning feature/insert for 'rN-oF', 'rN-oR', 'rN-iF', 'rN-iR' primers:")
            iF_list = sorted(
                args.selected_oligo.scan(ins.seq, 'left',  primer_size=primer_length_range, folder=temp_folder, us_seq=ins.us_seq, ds_seq=ins.ds_seq, time_limit=args.primer_scan_limit),
                key=lambda x: x.weight,
                reverse=True
            )
            self.logger.info('  len(iF_list) = {}'.format(len(iF_list)))
            self.logger.info('  insert_F: skipping {}/{} calculated primers'.format(max(0, len(iF_list)-subset_size), len(iF_list)))
            iF_list = iF_list[:subset_size]
            
            iR_list = sorted(
                args.selected_oligo.scan(ins.seq, 'right', primer_size=primer_length_range, folder=temp_folder, us_seq=ins.us_seq, ds_seq=ins.ds_seq, time_limit=args.primer_scan_limit),
                key=lambda x: x.weight,
                reverse=True
            )
            self.logger.info('  len(iR_list) = {}'.format(len(iR_list)))
            self.logger.info('  insert_R: skipping {}/{} calculated primers'.format(max(0, len(iR_list)-subset_size), len(iR_list)))
            iR_list = iR_list[:subset_size]
            
            insert_list.append([iF_list, iR_list])
        
        # Select the best subset of primers so the pairing doesn't take too long
        # This will cap the maximum number of primer pairs for each region span
        # to be 100 x 100 = 10,000
        
        # Do the pair calculations that involve 'sF' and 'sR'
        self.logger.info("Calculating: 'sF' 'sR' paired primers...")
        sF_sR_paired_primers = args.selected_oligo.pair(sF_list, sR_list, intervening=0, folder=temp_folder, time_limit=args.primer_pair_limit)
        
        for pp in sF_sR_paired_primers:
            # Re-calculate the PrimerPair weights to prefer the smallest amplicon sizes
            pp.weight = pp.get_weight(minimize=True)
            
            # Re-name the primers so when they are printed, they are easy to distinguish
            pp.forward_primer.set_name('sF')
            pp.reverse_primer.set_name('sR')
        
        # Sort by weight
        sF_sR_paired_primers = sorted(
            sF_sR_paired_primers,
            key=lambda x: x.get_joint_weight(),
            reverse=True
        )
        self.logger.info('  len(sF_sR_paired_primers) = {}'.format(len(sF_sR_paired_primers)))
        
        # pair_list[ins index] = [sF_oR_paired_primers, oF_sR_paired_primers, iF_iR_paired_primers]
        pair_list = []
        insert_pair_list = []
        
        #for i, (iF_list, iR_list) in enumerate(insert_list):
        for i in range(len(insert_seqs)):
            iF_list, iR_list = insert_list[i]
            ins = insert_seqs[i]
            
            self.logger.info('Pairing primers: ' + str(i))
        
            # If there is a hard constraint for primers that should be used
            # That is, if flanktag primers should be used
            if (len(args.mandatory_primers) > 0):
                pass
            # Otherwise, no flanktags are specified
            else:
                # sF rN-oR
                self.logger.info("  Calculating: 'sF' 'r"+str(ins.genome_r)+ins.type+"-oR' paired_primers...")
                sF_oR_paired_primers = sorted(
                    args.selected_oligo.pair(sF_list, iR_list, intervening=ins.fus_dist, folder=temp_folder, time_limit=args.primer_pair_limit),
                    key=lambda x: x.get_joint_weight(),
                    reverse=True
                )
                for pp in sF_oR_paired_primers:
                    pp.forward_primer.set_name('sF')
                    pp.reverse_primer.set_name('r'+str(ins.genome_r)+ins.type+'-oR')
                self.logger.info('  len(sF_oR_paired_primers) = {}'.format(len(sF_oR_paired_primers)))
            
                # rN-0F sR
                self.logger.info("  Calculating: 'r"+str(ins.genome_r)+ins.type+"-oF' 'sR' paired_primers...")
                oF_sR_paired_primers = sorted(
                    args.selected_oligo.pair(iF_list, sR_list, intervening=ins.fds_dist, folder=temp_folder, time_limit=args.primer_pair_limit),
                    key=lambda x: x.get_joint_weight(),
                    reverse=True
                )
                for pp in oF_sR_paired_primers:
                    pp.forward_primer.set_name('r'+str(ins.genome_r)+ins.type+'-oF')
                    pp.reverse_primer.set_name('sR')
                self.logger.info('  len(oF_sR_paired_primers) = {}'.format(len(oF_sR_paired_primers)))
                
                # rN-iF rN-iR
                self.logger.info("  Calculating: 'r"+str(ins.genome_r)+ins.type+"-iF' 'r"+str(ins.genome_r)+ins.type+"-iR' paired_primers...")
                iF_iR_paired_primers = sorted(
                    args.selected_oligo.pair(iF_list, iR_list, intervening=0, same_template=True, folder=temp_folder, time_limit=args.primer_pair_limit),
                    key=lambda x: x.get_joint_weight(),
                    reverse=True
                )
                for pp in iF_iR_paired_primers:
                    pp.forward_primer.set_name('r'+str(ins.genome_r)+ins.type+'-iF')
                    pp.reverse_primer.set_name('r'+str(ins.genome_r)+ins.type+'-iR')
                self.logger.info('  len(iF_iR_paired_primers) = {}'.format(len(iF_iR_paired_primers)))
                
                # Add calculated primer pairs to the list
                pair_list.append([sF_oR_paired_primers, oF_sR_paired_primers, ins])
                insert_pair_list.append(iF_iR_paired_primers)
        
        # If a primer in feature/insert of one round ALSO is identical
        # with a primer in another round,
        # Then it should only be weighted for a SINGLE entry
        starting_set, finished_set, round_labels = self.calculate_them_best_set(args, sF_sR_paired_primers, pair_list)
        
        d_starting_set, d_finished_set = self.calculate_amp_d_set(args, insert_pair_list)
        
        return finished_set, insert_pair_list, round_labels, d_finished_set
    
    def get_primer_location(self, filenames, genomes, primer_sequence):
        #[filename, contig, start, end, strand]
        locations = []
        for r in range(len(filenames)):
            locations.append([])
            for contig, seq in genomes[r].items():
                for m in regex.finditer(primer_sequence, seq, flags=regex.IGNORECASE, overlapped=True):
                    locations[-1].append([filenames[r], contig, m.start(), m.end(), '+'])
                for m in regex.finditer(nucleotides.rc(primer_sequence), seq, flags=regex.IGNORECASE, overlapped=True):
                    locations[-1].append([filenames[r], contig, m.start(), m.end(), '-'])
        return locations
    
    def round_counter(self):
        """
        Returns '0b', '0a', '1b', '1a', '2b', '2a', ...
        """
        r = 0
        m = 'b'
        while True:
            yield str(r) + m
            if (m == 'b'):
                m = 'a'
            else:
                m = 'b'
                r += 1
    
    def round_converter(self, count):
        """
        Converts '0b' to 0 (the index of the genome)
                 '0a'    1
                 '1b'    1
                 '1a'    2
        Returns int values
        """
        
        if (count[-1] == 'b'):
            return int(count[:-1])
        else:
            return int(count[:-1])+1
    
    def flatten_primer_pair_list(self, pair_list):
        p_set = []
        for pp in pair_list:
            if pp:
                p_set.append(pp.forward_primer)
                p_set.append(pp.reverse_primer)
            else:
                p_set.append(None)
                p_set.append(None)
        return p_set
    
    def optimize_pp_list_by_weight(self, args, pair_list, semirandom=True):
        """
        This will effectively minimize the number of primers sequences needed.
        """
        
        if semirandom:
            # Get ranked-random set
            pp_set = [self.random_primer_by_weight(pp_list) for pp_list in pair_list]
        else:
            # Get top-ranked set
            pp_set = [pp_list[0] if (len(pp_list) > 0) else None for pp_list in pair_list]
        
        p_group_weight = args.selected_oligo.p_group_weight(self.flatten_primer_pair_list(pp_set))
        pp_group_weight = args.selected_oligo.pp_group_weight(pp_set)
        joint_weight = p_group_weight * pp_group_weight
        
        starting_set = (joint_weight, pp_set) # initial set
        
        # iteratively improve until a local maxima is found
        # By continually swapping out the least-weighted component primer
        weight_delta = 1
        while (weight_delta > 0):
            # Get list of indices from smallest weight to largest weight (excluding 'sF' and 'sR')
            # If a primer is 'None', then it is right-most in the list.
            # Left-most element is the smallest weight. Weights increase as the index of the list increases.
            
            wi_order = self.rank_order([pp.get_joint_weight() if pp else math.inf for pp in pp_set])
            
            # If primer is 'None', then its weight will be 'math.inf', and the index
            # corresponding to it will be the last element in wi_order:
            #   rank_order([1, 2, 2.5, 0.001, math.inf]) # [3, 0, 1, 2, 4]
            
            # Start at the worst, and go to the next-worst, then next, then next
            for wi in wi_order:
                # Copy the list of primer pairs
                pp_set2 = pp_set[:]
                
                # Swap the worst-performing primer with an alternative.
                # If the alternative gives a better weight, then keep it
                # and break out of the loop
                for source_pp in pair_list[wi]:
                    pp_set2[wi] = source_pp
                    
                    new_p_group_weight = args.selected_oligo.p_group_weight(self.flatten_primer_pair_list(pp_set2))
                    new_pp_group_weight = args.selected_oligo.pp_group_weight(pp_set2)
                    new_joint_weight = new_p_group_weight * new_pp_group_weight
                    
                    weight_delta = new_joint_weight - joint_weight
                    
                    if (weight_delta > 0):
                        # replace values for next iteration
                        pp_set = pp_set2
                        joint_weight = new_joint_weight
                        self.logger.info('          t-set: ' + str((joint_weight, pp_set)))
                        break
                else:
                    # If no pp swap in 'source_pp' is higher-weighted, then return 0 rather than the last-calculated 'weight_delta'
                    weight_delta = 0
                
                # If the pp swap already produces a higher-weighted primer set, then skip the rest of the 'wi'
                if (weight_delta > 0):
                    break
        
        finished_set = (joint_weight, pp_set)
        
        return starting_set, finished_set
    
    def random_choices(self, population, weights, k=1):
        """
        Return a k sized list of population elements chosen with replacement.
        """
        weight_sum = sum(weights)
        choices = zip(population, weights)
        values = []
        for i in range(k):
            r = random.uniform(0, weight_sum)
            upto = 0
            for c, w in choices:
                if upto + w >= r:
                    values.append(c)
                    break
                upto += w
            else:
                values.append(random.choice(population))
        return values
    
    def random_primer_by_weight(self, pairs):
        if (len(pairs) == 0):
            return None
        elif ('choices' in random.__all__):
            return random.choices(pairs, [x.get_joint_weight() for x in pairs])[0]
        else:
            return self.random_choices(pairs, [x.get_joint_weight() for x in pairs])[0]
    
    def nsum(self, n):
        return n*(n+1)/2
    
    def rank_order(self, x, reverse=False, shift=0):
        return [y+shift for y in sorted(range(len(x)), key=x.__getitem__, reverse=reverse)]
    
    def rank(self, x, reverse=False, shift=0):
        return [z[0]+shift for z in sorted(enumerate(sorted(enumerate(x), key=lambda w: w[1], reverse=reverse)), key=lambda y: y[1][0])]
    
    def rank_dist(self, x):
        s = nsum(len(x))
        return [y/s for y in rank(x, reverse=False, shift=1)]
    
    #def product(x):
    #    z = 1
    #    for y in x:
    #        z *= y
    #    return z
    
    def filter_primer_pairs(self, pairs, forward=None, reverse=None):
        ooo = []
        
        if (forward and reverse):
            for pp in pairs:
                if ((pp.forward_primer.sequence == forward.sequence) and (pp.reverse_primer.sequence == reverse.sequence)):
                    ooo.append(pp)
        elif forward:
            for pp in pairs:
                if (pp.forward_primer.sequence == forward.sequence):
                    ooo.append(pp)
        elif reverse:
            for pp in pairs:
                if (pp.reverse_primer.sequence == reverse.sequence):
                    ooo.append(pp)
        else:
            ooo = pairs
        
        return ooo
    
    def pp_sets_equal(self, set1, set2):
        seqs1 = [(pp.forward_primer.sequence, pp.reverse_primer.sequence) if pp else (None, None) for pp in set1[1]]
        seqs2 = [(pp.forward_primer.sequence, pp.reverse_primer.sequence) if pp else (None, None) for pp in set2[1]]
        return seqs1 == seqs2
    
    def which_pp_equal(self, pp_set, shift=1):
        t_num = (len(pp_set)-1)//2
        for t1 in range(t_num):
            pp1L = pp_set[1+(2*t1)]
            if pp1L:
                pp1Ls = (pp1L.forward_primer.sequence, pp1L.reverse_primer.sequence)
            else:
                pp1Ls = None
            
            pp1R = pp_set[1+(2*t1)+1]
            if pp1R:
                pp1Rs = (pp1R.forward_primer.sequence, pp1R.reverse_primer.sequence)
            else:
                pp1Rs = None
            
            for t2 in range(t_num):
                if (t1 < t2):
                    pp2L = pp_set[1+(2*t2)]
                    if pp2L:
                        pp2Ls = (pp2L.forward_primer.sequence, pp2L.reverse_primer.sequence)
                    else:
                        pp2Ls = None
                        
                    pp2R = pp_set[1+(2*t2)+1]
                    if pp2R:
                        pp2Rs = (pp2R.forward_primer.sequence, pp2R.reverse_primer.sequence)
                    else:
                        pp2Rs = None
                    
                    Ltest = None if (pp1Ls == pp2Ls == None) else pp1Ls == pp2Ls
                    Rtest = None if (pp1Rs == pp2Rs == None) else pp1Rs == pp2Rs
                    self.logger.info('PrimerPair equals ({} vs {}): L={}, R={}'.format(t1, t2, Ltest, Rtest))
    
    def calculate_them_best_set(self, args, sF_sR_paired_primers, pair_list, max_iterations=10000):
        """
        Takes input primer sets and calculates their weights
        
        This starts with the best from a sorted list of 'sF' 'sR' primer pairs
        Then it optimizes by swapping the worst component a finite number of times
        (or until the delta drops below a certain amount)
        """
        
        # This should be in a loop...
        #sF_oR_paired_primers, oF_sR_paired_primers, iF_iR_paired_primers, ins = pair_list[0]
        
        # Set initial iteration count to zero
        iteration_count = 0
        
        # Create lists to hold n-set of primer pairs, and their joint weight
        starting_set = [] # initial set
        finished_set = [] # final set
        
        # Make a list that describes the order of the output
        round_labels = [] # ['0b', '0a', '1b', '1a', '0b',      '0b',      '0a',      '0a',      '1b',      '1b',      '1a',      '1a']
        #                    ^^^^^^^^ sF sR ^^^^^^^, sF r0b-oR, r0b-oF sR, sF r0b-oR, r0b-oF sR, sF r0b-oR, r0b-oF sR, sF r0b-oR, r0b-oF sR
        #                    'A',  'A',  'A',  'A',  'B',       'C',       'B',       'C',       'B',       'C',       'B',       'C'
        
        #for sF_oR_paired_primers, oF_sR_paired_primers, ins in pair_list:
        #    round_labels.append(('A', str(ins.genome_r)+ins.type))
        round_labels.append(('A', None))
        
        for sF_oR_paired_primers, oF_sR_paired_primers, ins in pair_list:
            for amplicon_side in ['B', 'C']:
                round_labels.append((amplicon_side, str(ins.genome_r)+ins.type))
        
        
        sF_sR_iterator = iter(sF_sR_paired_primers)
        
        if ((max_iterations != None) and (len(sF_sR_paired_primers) > max_iterations)):
            self.logger.info('sF_sR_paired_primers: skipping {}/{} calculated primer pairs'.format(max(0, len(sF_sR_paired_primers)-max_iterations), len(sF_sR_paired_primers)))
        
        while (iteration_count < max_iterations):
            # Increment the count
            iteration_count += 1
            
            # Create random 6-sets, with each pair's probability determined by its joint weight
            #uf_dr_pair = random_primer_by_weight(sF_sR_paired_primers) # comment this out to prevent random sampling of this one
            try:
                sF_sR_pair = next(sF_sR_iterator)
            except StopIteration:
                break
            
            self.logger.info('loop {}:'.format(iteration_count))
            
            # Populate set with all primer pairs that have 'sF' and 'sR'
            pp_sources = []
            
            for sF_oR_paired_primers, oF_sR_paired_primers, ins in pair_list:
                # Add the 'oR' from 'sF-oR' pair. These should already be sorted by weight from highest to lowest.
                pp_sources.append(self.filter_primer_pairs(sF_oR_paired_primers, forward=sF_sR_pair.forward_primer))
                
                # Add the 'oF' from 'oF-sR' pair. These should already be sorted by weight from highest to lowest.
                pp_sources.append(self.filter_primer_pairs(oF_sR_paired_primers, reverse=sF_sR_pair.reverse_primer))
            
            self.logger.info('  length of sources: ' + str([len(x) for x in pp_sources]))
            
            # Only proceed with calculations (finally adding this primer to
            # 'starting_set' and 'finished_set') if it has potential internal
            # primers available, according to 'args.internal_primers_required'
            
            # Turns list ['y', 'n', 'y'] into ['y', 'y', 'n', 'n', 'y', 'y']
            if args.internal_primers_required:
                required_pattern = [val for val in args.internal_primers_required for b in range(2)]
                current_pp_sources = [len(x) for x in pp_sources]
                should_mask = False
                should_skip = False
                for req_str, num_pp in zip(required_pattern, current_pp_sources):
                    if ((req_str in ['y', 'Y', '1', 'T', 't', 'TRUE', 'True', 'true']) and (num_pp == 0)):
                        should_skip = True
                        break
                    if ((req_str not in ['y', 'Y', '1', 'T', 't', 'TRUE', 'True', 'true']) and (num_pp > 0)):
                        should_mask = True
                if should_skip:
                    self.logger.info('  skipping...')
                    continue
            
            # If required_list=['y', 'n', 'n', 'y']
            # And current_pp_sources=[4, 3, 1, 0, 1, 0, 10, 9]
            # Then by default the program will always try to include the 1s, which aren't required
            # To make it okay for the problem to ignore these, then we can either
            #  1) Duplicate, then mask current_pp_sources so it artificially looks like this: [4, 3, 0, 0, 0, 0, 10, 9]
            #  2) or we can modify the algorithm to allow discarding elements if there is a 'n' in 'required_pattern'
            # Let's go with (1)
            
            # Oldish code to test
            test_starting, test_final = self.optimize_pp_list_by_weight(args, [[sF_sR_pair]]+pp_sources, semirandom=False)
            self.logger.info('  test-starting: ' + str(test_starting))
            self.logger.info('     test-final: ' + str(test_final))
            # Report whether the "optimal" PrimerPairs chosen are identical
            self.which_pp_equal(test_final[1])
            starting_set.append(test_starting)
            finished_set.append(test_final)
            
            # Brand new code to test
            pp_set = PrimerDesign([[sF_sR_pair]]+pp_sources)
            for oo in range(2):
                pp_set.optimize()
            self.logger.info('NEW >= OLD: ' + str(pp_set.optimal.weight >= test_final[0]))
            
            
            # If necessary, we "mask" the non-required PrimerPairs, then calculate.
            if args.internal_primers_required:
                if should_mask:
                    masked_pp_sources = []
                    for req_str, pp_s in zip(required_pattern, pp_sources):
                        if (req_str in ['y', 'Y', '1', 'T', 't', 'TRUE', 'True', 'true']):
                            masked_pp_sources.append(pp_s)
                        else:
                            masked_pp_sources.append([])
                    masked_starting, masked_final = self.optimize_pp_list_by_weight(args, [[sF_sR_pair]]+masked_pp_sources, semirandom=False)
                    self.logger.info('masked-starting: ' + str(masked_starting))
                    self.logger.info('   masked-final: ' + str(masked_final))
                    # Report whether the "optimal" PrimerPairs chosen are identical
                    self.which_pp_equal(masked_final[1])
                    starting_set.append(masked_starting)
                    finished_set.append(masked_final)
                    
                    # Brand new code to test
                    pp_set = PrimerDesign([[sF_sR_pair]]+masked_pp_sources)
                    for oo in range(2):
                        pp_set.optimize()
                    self.logger.info('NEW >= OLD: ' + str(pp_set.optimal.weight >= masked_final[0]))
            
            # Report the number of PrimerPairs that are shared among the inserts
            t_num = len(pp_sources)//2
            for t1 in range(t_num):
                t1L_set = set([(x.forward_primer.sequence, x.reverse_primer.sequence) for x in pp_sources[2*t1] if x])
                t1R_set = set([(x.forward_primer.sequence, x.reverse_primer.sequence) for x in pp_sources[(2*t1)+1] if x])
                for t2 in range(t_num):
                    if (t1 < t2):
                        t2L_set = set([(x.forward_primer.sequence, x.reverse_primer.sequence) for x in pp_sources[2*t2] if x])
                        t2R_set = set([(x.forward_primer.sequence, x.reverse_primer.sequence) for x in pp_sources[(2*t2)+1] if x])
                        tL_num = len(t1L_set.intersection(t2L_set))
                        tL_den = len(t1L_set.union(t2L_set))
                        tR_num = len(t1R_set.intersection(t2R_set))
                        tR_den = len(t1R_set.union(t2R_set))
                        self.logger.info('number shared insert PrimerPairs ({} vs {}): L={}/{}, R={}/{}'.format(t1, t2, tL_num, tL_den, tR_num, tR_den))
        
        return starting_set, finished_set, round_labels
    
    def which_d_equal(self, pp_set):
        """
        Does all pairwise comparisons of PrimerPairs within list 'pp_set'
        to identify which pairs have identical sequences.
        """
        for i, ppi in enumerate(pp_set):
            if ppi:
                i_seq = (ppi.forward_primer.sequence, ppi.reverse_primer.sequence)
            else:
                i_seq = (None, None)
            
            for j, ppj in enumerate(pp_set):
                if (i < j):
                    if ppj:
                        j_seq = (ppj.forward_primer.sequence, ppj.reverse_primer.sequence)
                    else:
                        j_seq = (None, None)
                    
                    comparison = None if (i_seq == j_seq == (None, None)) else (i_seq == j_seq)
                    self.logger.info("Amp 'D' PrimerPair {} equals PrimerPair {}: {}".format(i, j, comparison))
    
    def calculate_amp_d_set(self, args, insert_pair_list, max_iterations=1000):
        iteration_count = 0
        
        starting_set = []
        finished_set = [] # final set
        
        while (iteration_count < max_iterations):
            # Increment the count
            iteration_count += 1
            self.logger.info('loop {}:'.format(iteration_count))
            
            self.logger.info('  length of sources: ' + str([len(x) for x in insert_pair_list]))
            if args.internal_primers_required:
                current_pp_sources = [len(x) for x in insert_pair_list]
                should_mask = False
                should_skip = False
                for req_str, num_pp in zip(args.internal_primers_required, current_pp_sources):
                    if ((req_str in ['y', 'Y', '1', 'T', 't', 'TRUE', 'True', 'true']) and (num_pp == 0)):
                        should_skip = True
                        break
                    if ((req_str not in ['y', 'Y', '1', 'T', 't', 'TRUE', 'True', 'true']) and (num_pp > 0)):
                        should_mask = True
                if should_skip:
                    self.logger.info('  skipping...')
                    continue
                
            d_starting, d_final = self.optimize_pp_list_by_weight(args, insert_pair_list, semirandom=True)
            self.logger.info('  test-starting: ' + str(d_starting))
            self.logger.info('    Amp-D-final: ' + str(d_final))
            # Report whether the "optimal" PrimerPairs chosen are identical
            self.which_d_equal(d_final[1])
            starting_set.append(d_starting)
            finished_set.append(d_final)
            
            # Brand new code to test
            pp_set = PrimerDesign(insert_pair_list)
            for oo in range(2):
                pp_set.optimize()
            self.logger.info('NEW >= OLD: ' + str(pp_set.optimal.weight >= d_final[0]))
            
            # If necessary, we "mask" the non-required PrimerPairs, then calculate.
            if args.internal_primers_required:
                if should_mask:
                    masked_insert_pair_list = []
                    for req_str, pp_list in zip(args.internal_primers_required, insert_pair_list):
                        if (req_str in ['y', 'Y', '1', 'T', 't', 'TRUE', 'True', 'true']):
                            masked_insert_pair_list.append(pp_list)
                        else:
                            masked_insert_pair_list.append([])
                    masked_starting, masked_final = self.optimize_pp_list_by_weight(args, masked_insert_pair_list, semirandom=True)
                    self.logger.info('masked-starting: ' + str(masked_starting))
                    self.logger.info('   masked-final: ' + str(masked_final))
                    # Report whether the "optimal" PrimerPairs chosen are identical
                    self.which_d_equal(masked_final[1])
                    starting_set.append(masked_starting)
                    finished_set.append(masked_final)
                    
                    # Brand new code to test
                    pp_set = PrimerDesign(masked_insert_pair_list)
                    for oo in range(2):
                        pp_set.optimize()
                    self.logger.info('NEW >= OLD: ' + str(pp_set.optimal.weight >= masked_final[0]))
        
        return starting_set, finished_set
    
    def filter_alignment_records_for_cPCR(self, args, alignment_filename, dDNA_contigs):
        """
        Parse the alignment to identify the exogenous DNA, as well as pairs of left/right flanking homology regions
        Note: This will only allow for each dDNA to align to a single place on each contig.
        If a single dDNA could target multiple loci per contig, then this function will not suffice
        Thus, this code needs modifications for broader applications.
        """
        # Allow for alignments to omit up to 3 nt at edges (of homology regions)
        permitted_edge = 3
        
        # Make a set of invalid records that will be removed because they failed some quality threshold
        invalid_records = set()
        
        # Make dict of all decent records
        records = {}
        
        # Read the alignment file
        for record in args.selected_aligner.load(alignment_filename):
            # Process the record
            # Require certain e-value and length for a "significant" alignment
            if ((record.evalue <= args.max_evalue) and (record.length >= args.min_length)):
                # Require alignment to occur at the termini of the query (either at the extreme beginning, or extreme end)
                #if ((record.query_position[0] < permitted_edge) or (record.query_position[1] > len(dDNA_contigs[record.query_name])-permitted_edge)):
                #    records.setdefault((record.query_name, record.subject_name), []).append(record)
                qs_pair = (record.query_name, record.subject_name)
                
                # If the alignment is on the edges of the query
                if ((record.query_position[0] < permitted_edge) or (record.query_position[1] > len(dDNA_contigs[record.query_name])-permitted_edge)):
                    if (qs_pair not in records):
                        records[qs_pair] = [None, None]
                
                # Arrange the records according to their flanking: records[(qname, sname)] = [left_record, right_record]
                if (record.query_position[0] < permitted_edge): # left
                    if (records[qs_pair][0] == None):
                        records[qs_pair][0] = record
                    else:
                        self.logger.info(str(qs_pair) + ' has too many valid alignments')
                        invalid_records.add(qs_pair)
                
                elif (record.query_position[1] > len(dDNA_contigs[record.query_name])-permitted_edge): # right
                    if (records[qs_pair][1] == None):
                        records[qs_pair][1] = record
                    else:
                        self.logger.info(str(qs_pair) + ' has too many valid alignments')
                        invalid_records.add(qs_pair)
        
        self.logger.info("before filtering:")
        for kkk, vvv in records.items():
            self.logger.info("  "+str(kkk) + " " + str(vvv))
        
        # Queue for removal the (query, subject) pairs that only have 1 of the 2 required alignments
        for qs_pair, record_pair in records.items():
            if (None in record_pair):
                invalid_records.add(qs_pair)
        
        # Remove the (query, subject) pairs that had too many valid alignments
        for bad_key in invalid_records:
            records.pop(bad_key)
            self.logger.info('Removing invalid record: ' + str(bad_key))
        
        self.logger.info("after filtering:")
        for kkk, vvv in records.items():
            self.logger.info("  "+str(kkk) + " " + str(vvv))
        
        return records
    
    def filter_alignment_records2(self, args, alignment_filename, dDNA_contigs):
        records = {}
        
        # Allow for alignments to omit up to 3 nt at edges (of homology regions)
        permitted_edge = 3
        
        # Make a set of invalid records that will be removed because they failed some quality threshold
        invalid_records = set()
        
        # Read the alignment file
        for record in args.selected_aligner.load(alignment_filename):
            # Process the record
            # Require certain e-value and length for a "significant" alignment
            if ((record.evalue <= args.max_evalue) and (record.length >= args.min_length)):
                # Require alignment to occur at the termini of the query (either at the extreme beginning, or extreme end)
                #if ((record.query_position[0] < permitted_edge) or (record.query_position[1] > len(dDNA_contigs[record.query_name])-permitted_edge)):
                #    records.setdefault((record.query_name, record.subject_name), []).append(record)
                qs_pair = (record.query_name, record.subject_name)
                
                # The entire query should align to the subject
                if (record.query_position[0] < permitted_edge < len(dDNA_contigs[record.query_name])-permitted_edge < record.query_position[1]):
                    if (qs_pair not in records):
                        records[qs_pair] = record
                    else:
                        self.logger.info(str(qs_pair) + ' has too many valid alignments')
                        invalid_records.add(qs_pair)
        
        for bad_qs_pair in invalid_records:
            records.pop(bad_qs_pair)
            self.logger.info('Removing invalid record: ' + str(bad_qs_pair))
        
        return records
    
    def pathstrip(self, path):
        return os.path.splitext(os.path.basename(path))[0]

class PrimerWorker(multiprocessing.Process):
    def __init__(self, args, cutoffs, queue, results_queue, log_lock, data_lock):
        super().__init__()
        self.args = args
        self.cutoffs = copy.deepcopy(cutoffs)
        
        self.queue = queue
        self.results_queue = results_queue
        self.log_lock = log_lock
        self.data_lock = data_lock
    
    def run(self):
        # Will need to change temp directory to be process-specific
        # This must be done after the process has been invoked (i.e. in the run method)
        self.cutoffs['folder'] = os.path.join(self.cutoffs['folder'], str(self.pid))
        
        while True:
            # Get the next 'Primer' object from the queue
            p = self.queue.get()
            
            # 'None' is the poison pill that terminates this process
            if (p == None):
                self.queue.task_done()
                break
            
            cpass = p.summarize(p.checks)
            if ((p.checks[0] == None) or (not cpass)):
                p.progressive_check(self.cutoffs)
                self.results_queue.put((p.sequence, p.mp_dump()))
            
            self.queue.task_done()
        
        return None

class PrimerPairWorker(multiprocessing.Process):
    def __init__(self, args, cutoffs, queue, results_queue, log_lock, data_lock):
        super().__init__()
        self.args = args
        self.cutoffs = copy.deepcopy(cutoffs)
        
        self.queue = queue
        self.results_queue = results_queue
        self.log_lock = log_lock
        self.data_lock = data_lock
    
    def run(self):
        # Will need to change temp directory to be process-specific
        # This must be done after the process has been invoked (i.e. in the run method)
        self.cutoffs['folder'] = os.path.join(self.cutoffs['folder'], str(self.pid))
        
        while True:
            # Get the next 'PrimerPair' object from the queue
            pp = self.queue.get()
            
            # 'None' is the poison pill that terminates this process
            if (pp == None):
                self.queue.task_done()
                break
            
            # Run checks
            pp.progressive_check(self.cutoffs)
            
            pair = (pp.forward_primer.sequence, pp.reverse_primer.sequence)
            self.results_queue.put((pair, pp.mp_dump()))
            
            self.queue.task_done()
        
        return None

class PrimerDesignWorker(multiprocessing.Process):
    def __init__(self, args, queue, results_queue, log_lock, data_lock):
        super().__init__()
        self.args = args
        self.queue = queue
        self.results_queue = results_queue
        self.log_lock = log_lock
        self.data_lock = data_lock
    
    def run(self):
        while True:
            # Get the next 'Primer' object from the queue
            d = self.queue.get()
            
            # 'None' is the poison pill that terminates this process
            if (d == None):
                self.queue.task_done()
                break
            
            d.optimize(mode='direct', lock=self.log_lock)
            #design.optimize(iterations=3000, lock=self.log_lock)
            
            self.results_queue.put(d.optimal)
            
            self.queue.task_done()
        
        return None

class Datum(object):
    """
    Class representing the location of homology regions and inserts within
    contigs.
    The general pattern for each locus is:
      Datum(dDNA_r=1, ..., genome_r=0, ...)
      Datum(dDNA_r=1, ..., genome_r=1, ...)
      Datum(dDNA_r=2, ..., genome_r=1, ...)
      Datum(dDNA_r=2, ..., genome_r=2, ...)
    """
    __slots__ = [
        'dDNA_r',
        'dDNA_contig',
        'genome_r',
        'genome_contig',
        'ush_start',
        'ush_end',
        'dsh_start',
        'dsh_end',
        'ins_start',
        'ins_end'
    ]
    
    def __init__(self, dDNA_r, dDNA_contig, genome_r, genome_contig, ush_start, ush_end, dsh_start, dsh_end, ins_start, ins_end):
        self.dDNA_r = dDNA_r
        self.dDNA_contig = dDNA_contig
        self.genome_r = genome_r
        self.genome_contig = genome_contig
        self.ush_start = ush_start
        self.ush_end = ush_end
        self.dsh_start = dsh_start
        self.dsh_end = dsh_end
        self.ins_start = ins_start
        self.ins_end = ins_end
        
    def __repr__(self):
        return self.__class__.__name__ + '(' + ', '.join([s + '=' + repr(getattr(self, s)) for s in self.__slots__]) + ')'

class Insert(object):
    """
    Class representing the insert region of a genome editing event, with
    far-upstream and far-downstream distances.
    """
    __slots__ = [
        'genome_r',
        'genome_contig',
        'seq',
        'us_seq',
        'ds_seq',
        'fus_dist',
        'fds_dist',
        'type'
    ]
    
    def __init__(self, genome_r, genome_contig, seq, us_seq, ds_seq, fus_dist, fds_dist, type):
        self.genome_r = genome_r
        self.genome_contig = genome_contig
        self.seq = seq
        self.us_seq = us_seq
        self.ds_seq = ds_seq
        self.fus_dist = fus_dist
        self.fds_dist = fds_dist
        self.type = type # 'b' for before, 'a' for after
    
    def __repr__(self):
        return self.__class__.__name__ + '(' + ', '.join([s + '=' + repr(getattr(self, s)) for s in self.__slots__]) + ')'

class Cutoff(object):
    # List holding all currently-defined cutoffs
    cutoffs = []
    
    # Index of the cutoff to modify next
    ci = 0
    
    def __init__(self, name, initial, final, delta, separate=False):
        self.name = name
        self.initial = initial
        self.final = final
        self.delta = delta
        self.current = initial
        
        self.separate=separate
        self.separate_current=0
        
        self.cutoffs.append(self)
    
    @classmethod
    def _increment(cls):
        # Increment
        cls.ci += 1
        if (cls.ci >= len(cls.cutoffs)):
            cls.ci = 0
    
    @classmethod
    def loosen(cls):
        cutoff = cls.cutoffs[cls.ci]
        
        new = []
        for i, (c, d, f) in enumerate(zip(cutoff.current, cutoff.delta, cutoff.final)):
            if cutoff.separate:
                if (i == cutoff.separate_current):
                    n = c+d
                    if math.isclose(n, f):
                        n = f
                    elif ((d < 0) and (n < f)):
                        n = f
                    elif ((d > 0) and (n > f)):
                        n = f
                else:
                    n = c
                
            else:
                n = c+d
                if math.isclose(n, f):
                    n = f
                if ((d < 0) and (n < f)):
                    n = f
                elif ((d > 0) and (n > f)):
                    n = f
                
            new.append(n)
        new = tuple(new)
        
        if cutoff.separate:
            cutoff.separate_current += 1
            if (cutoff.separate_current >= len(cutoff.current)):
                cutoff.separate_current = 0
        
        
        # If there was no change, try again with the same CI
        if ((new == cutoff.current) and (new != cutoff.final)):
            finished = cls.loosen()
        
        
        elif (new == cutoff.current):
            if all(x.current == x.final for x in Cutoff.cutoffs):
            #if all(math.isclose(x.current[i], x.final[i]) for x in Cutoff.cutoffs for i in range(len(x.current))): # 2
                # Then 'final' was reached for all cutoffs
                finished = True
            else:
                cls._increment()
                finished = cls.loosen()
        else:
            # Replace current
            cutoff.current = new
            
            # Increment
            cls._increment()
            
            # Set return value
            finished = False
        
        return finished
    
    @classmethod
    def test(cls):
        clist = [
            #   Name, initial, final, delta # states
            #cls("length_min", (19,), (19,), (0,)), # 1
            #cls("length_max", (28,), (36,), (2,)), # 5
            cls("length", (19,28), (19,36), (0,2), separate=True), # 1,5
            #cls("last5gc_min", (1,), (0,), (-1,)), # 2
            #cls("last5gc_max", (3,), (4,), (1,)), # 2
            cls("last5gc", (1,3), (0,4), (-1,1), separate=True), # 2
            #cls("gcclamp_min", (1,), (0,), (-1,)), # 2
            #cls("gcclamp_max", (2,), (4,), (1,)), # 3
            cls("gcclamp", (1,2), (0,4), (-1,1), separate=True), # 2
            cls("gcfreq", (0.4, 0.6), (0.2, 0.8), (-0.1, 0.1)), # 3
            cls("runlen_max", (4,), (6,), (1,)), # 3
            cls("3primecomplen_max", (3,), (6,), (1,)), # 4
            cls("deltag_min", (-4.0,), (-7.0,), (-1.0,)), # 4
            cls("tm", (52, 65), (52, 65), (0, 0)), # 1
            # cls("amplicon", (300,700), (300,700), (0,0)), # 1
            cls("deltatm_max", (2.5,), (4.0,), (0.5,)), # 4
        ]
        finished = False
        i = 0
        while (finished == False):
            
            print(i, [x.current for x in clist], end=' ')
            finished = cls.loosen()
            print(finished)
            
            # Convert to dict
            d = {c.name: c.current if (len(c.current) > 1) else c.current[0] for c in Cutoff.cutoffs}
            print(' ', d)
            
            i += 1
    
class CutoffIterator(object):
    def __init__(self):
        pass
    
    def __iter__(self):
        return self
    
    def __next__(self):
        d = {c.name: c.current if (len(c.current) > 1) else c.current[0] for c in Cutoff.cutoffs}
        finished = Cutoff.loosen()
        if finished:
            raise StopIteration
        else:
            return d
    
    
        
