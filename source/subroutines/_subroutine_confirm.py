#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/subroutines/_subroutine_confirm.py

# Import standard packages
import os
import logging
import copy
import random
import math
from collections import namedtuple

# Import non-standard packages
import regex

# Import included AddTag-specific modules
from . import subroutine
from ..feature import Feature
from .. import utils
from .. import nucleotides
from .. import aligners
from .. import thermodynamics

class ConfirmParser(subroutine.Subroutine):
    def __init__(self, subparsers):
        self.subparsers = subparsers
        
        self.name = 'confirm'
        self.description = (
            "description:" "\n"
            "  Design primers for confirming whether each step of genome engineering is" " \n"
            "  successful or not. This does not design multiplexable primers." "\n"
        )
        self.help = "Design primers for confirming whether each step of genome engineering is successful or not."
        self.epilog = (
            "example:" "\n"
            "  You can list CRISPR/Cas motifs by running:" "\n"
            "   $ python3 {__program__} {__subroutine__}" "\n"
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
        
        aligner_choices = [x.name for x in aligners.aligners]
        self.parser.add_argument("--aligner", type=str, choices=aligner_choices, default='blastn',
            help="Program to calculate pairwise alignments. Please note that the 'addtag' internal aligner is very slow.")
        
        #self.parser.add_argument("--number_pcr_conditions", metavar="N", type=int, default=None,
        #    help="Number of PCR conditions to develop primers for. All amplicons \
        #    within each condition will have similar size (nt) and melting temperatures. \
        #    If unspecified, will default to the number of target features.")
        
        self.parser.add_argument("--max_number_designs_reported", metavar="N", type=int, default=5,
            help="The maximum number of final cPCR designs to report.")
        
        self.parser.add_argument("--primer_scan_limit", metavar="N", type=int, default=2*60,
            help="Number of seconds to limit each primer scan.")
        
        self.parser.add_argument("--primer_pair_limit", metavar="N", type=int, default=5*60,
            help="Amount of time (in seconds) to limit primer pairings.")
        
        self.parser.add_argument("--primers", nargs="+", metavar="*.fasta", type=str, default=[],
            help="A FASTA file for each round containing primer sequences that \
            must be used for that round. Usually, these correspond to flanktags.")
        
        oligo_choices = [x.name for x in thermodynamics.oligos]
        self.parser.add_argument("--oligo", type=str, choices=oligo_choices, default='UNAFold',
            help="Program to perform thermodynamic calculations.")
        
        self.parser.add_argument("--max_evalue", metavar="N", type=float, default=0.001,
            help="The maximum accepted alignment E-value to consider a recombination.")
        
        self.parser.add_argument("--min_length", metavar="N", type=int, default=35,
            help="The minimum accepted alignment length to consider a recombination.")
        
        self.parser.add_argument("--skip_round", metavar="N", nargs="+", type=int, default=[],
            help="Skip primer calculations for these rounds.")
        
        # Nucleotide matching stuff
        #  - number errors (for fuzzy regex)
        
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
        
        #Datum = namedtuple('Datum', ['r', 'qname', 'sname', 'us_seq', 'ds_seq', 'insert_seq', 'feature_seq', 'q_hih_seq', 'q_ush_seq', 's_ush_seq', 'q_dsh_seq', 's_dsh_seq'])
        #Datum = namedtuple('Datum', ['r', 'sname', 'ush_start', 'ush_end', 'dsh_start', 'dsh_end'])
        Datum = namedtuple('Datum', ['dDNA_r', 'dDNA_contig', 'genome_r', 'genome_contig', 'ush_start', 'ush_end', 'dsh_start', 'dsh_end', 'ins_start', 'ins_end'])
        
        # Create 'r0'
        logging.info('Working on round r{}'.format(0))
        # Parse the input FASTA
        genome_contigs_list.append(utils.load_multiple_fasta_files(args.fasta))
        dDNA_contigs_list.append(None) # The 'r0' round of genome engineering has no dDNA
        dDNA_fasta_file_list.append(None)
        dDNA_alignment_file_list.append(None)
        
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
            logging.info('Working on round r{}'.format(r))
            
            # Load input dDNA FASTA into the list of dicts
            dDNA_contigs_list.append(utils.old_load_fasta_file(dDNA_filename))
            
            # Build an index of the previous round's genome
            genome_index_file = args.selected_aligner.index(genome_fasta_file_list[r-1], os.path.basename(genome_fasta_file_list[r-1]), args.folder, args.processors)
            
            # Write current dDNA to FASTA
            dDNA_fasta_file_list.append(utils.write_merged_fasta(dDNA_contigs_list[r], os.path.join(args.folder, 'dDNA-r'+str(r)+'.fasta')))
            
            # We align dDNA against the genome+dDNA(prev)
            dDNA_alignment_file_list.append(args.selected_aligner.align(dDNA_fasta_file_list[r], genome_index_file, 'dDNA-r'+str(r)+'-alignment.'+args.selected_aligner.output, args.folder, args.processors))
            
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
                        logging.info('Working on {} vs {}'.format(qname, sname))
                        logging.info("crossover: " + str((i, qname, sname, record_list)))
                        
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
                            logging.info('  us_record is reverse complemented')
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
                                    m = regex.search('(?:'+q_overlap+'){e<'+str(len(q_overlap)//2)+'}$', s_ush_seq, flags=regex.IGNORECASE|regex.BESTMATCH)
                                    if m:
                                        sush_end = sush_start + m.start()
                                    
                                    # Set the query upstream homology end
                                    qush_end = qdsh_start
                                # If subject overlaps
                                elif (sdsh_start < sush_end):
                                    # Calculate the query trim
                                    s_overlap = scontig[sdsh_start:sush_end]
                                    m = regex.search('(?:'+s_overlap+'){e<'+str(len(s_overlap)//2)+'}$', q_ush_seq, flags=regex.IGNORECASE|regex.BESTMATCH)
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
                        #genome_contigs_list[r][sname] = us_seq + q_hih_seq + ds_seq # <---------------------- this should be done in the 'r' outer loop (earlier draft)
                        
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
                        # This code only works for a single-locus per chromosome! <------ need to update this so it works for multiple loci per chromosome
                        #updated_contigs[sname] = us_seq + q_hih_seq + ds_seq
                        ### end alternate ###
                        
                        # Add this parsed record data to the 'group_links'
                        # Key refers to dDNA file and contig name, value referes to gDNA file, contig name, and homology regions
                        #group_links.append(Datum(r, qname, r-1, sname, sush_start, sush_end, sdsh_start, sdsh_end, us_record.query_position[1], ds_record.query_position[0]))
                        group_links.append(Datum(r, qname, r-1, sname, sush_start, sush_end, sdsh_start, sdsh_end, qush_end, qdsh_start))
                        
                else:
                    logging.info(qname + ' vs ' + sname + ' has ' + str(len(record_list)) + ' regions of homology (2 needed)')
            
            # Rename the contigs based on the modifications
            # This code only works for a single-locus per chromosome! <------ need to update this so it works for multiple loci per chromosome
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
            logging.info('Working on round r{}'.format(r))
            
            # Align each dDNA to the genome it produced
            genome_index_file = args.selected_aligner.index(genome_fasta_file_list[r], os.path.basename(genome_fasta_file_list[r]), args.folder, args.processors)
            temp_alignment = args.selected_aligner.align(dDNA_fasta_file_list[r], genome_index_file, 'dDNA-temp-r'+str(r)+'-alignment.'+args.selected_aligner.output, args.folder, args.processors)
            
            # Parse the alignment to make an intermediary record data structure
            records = self.filter_alignment_records2(args, temp_alignment, dDNA_contigs_list[r])
            
            # Find the exact locations of each of the dDNAs within the engineered genomes
            # Iterate through the sorted records
            for i, qs_pair in enumerate(records):
                qname, sname = qs_pair
                record = records[qs_pair]
                
                logging.info('Working on {} vs {}'.format(qname, sname))
                
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
                    logging.info('  us_record is reverse complemented')
                    # We must reverse-complement these dDNA substrings so their orientation matches the genome
                    q_hih_seq = nucleotides.rc(q_hih_seq)
                else: # not reverse-complemented
                    pass
            
                #group_links.setdefault((r, qname), []).append(Datum(r, sname, sush_start, sush_end, sdsh_start, sdsh_end))
                group_links.append(Datum(r, qname, r, sname, sush_start, sush_end, sdsh_start, sdsh_end, ins_start, ins_end))
        
        ######## Begin for linking the loci ########
        
        logging.info('contig_groups:')
        for cg in contig_groups:
            logging.info('  ' + str(cg))
        
        def in_same_contig_group(contig1, contig2, cgroups=contig_groups):
            for g in cgroups:
                if ((contig1 in g) and (contig2 in g)):
                    return True
            else:
                return False
        
        # Make a copy of 'group_links' that we can pop stuff from
        datum_list = copy.deepcopy(group_links)
        
        logging.info('datum_list:')
        for datum in datum_list:
            logging.info('  ' + str(datum))
        
        
        # For each (round, locus), we need the (contig name, upstream_homology_length, downstream_homology_length)
        # loop through these:
        
        # For each cross-over event, there will be a datum
        # Typically, there should be 1 Datum object for haploid genomes,
        # and 2 Datum objects for diploid genomes
        
        
        # Populate initial groups based on if there is crossing over in the genome-r0
        datum_groups = []
        to_pop = []
        for di, datum in enumerate(datum_list):
            if datum.genome_contig in genome_contigs_list[0].keys():
                datum_groups.append([datum])
                to_pop.append(di)
        for di in sorted(to_pop, reverse=True):
            datum_list.pop(di)
        to_pop = []
        
        logging.info('datum_groups (r0):')
        for dg in datum_groups:
            logging.info('  ' + str(dg))
        
        logging.info('datum_list:')
        for datum in datum_list:
            logging.info('  ' + str(datum))
        
        # Populate the intermediary groups
        # New try:
        # If there is ANY OVERLAP AT ALL between
        #      round n,    AFTER: ush_start-dsh_end
        # and  round n+1, BEFORE: ush_start-dsh_end
        # Then these loci should be linked!
        no_more_after = False
        #no_more_before = False
        groups = []
        wcount = 0
        while (not no_more_after):
            # Specify an error if for some reason this doesn't work
            wcount += 1
            if (wcount > 10000):
                raise Exception('Too many iterations taken to link loci between engineered genomes.')
            
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
                    #else:
                    #    no_more_before = True
                
                else:
                    # If there is no "before" to pair with selected "after"
                    # then this is a terminal genome engineering
                    # Not sure if this code would work if a locus was "skipped" between rounds as follows:
                    #    genome-r0: locus A wt
                    #    genome-r1: locus A ko
                    #    genome-r2: locus A ko (no change)
                    #    genome-r3: locus A ki
                    
                    # Since this is a terminal, then we need to select the correct "before" using homology
                    pass
                if ((datum_after != None) and (datum_before != None)):
                    # Now that the two datum are selected, determine whether or not they overlap
                    # Feature.overlap_coverage(start1, end1, start2, end2, index_base=0)
                    overlap = Feature.overlap_coverage(datum_after.ush_start, datum_after.dsh_end, datum_before.ush_start, datum_before.dsh_end)
                    
                    # If they overlap
                    if (overlap > 0):
                        # Then add them to the same group
                        group = [datum_after, datum_before]
                        groups.append(group)
                        
                        logging.info('group:')
                        for g in group:
                            logging.info('  ' + str(g))
                        
                        for di in sorted(to_pop, reverse=True):
                            datum_list.pop(di)
                    # Otherwise, they are non-overlapping engineering events on the same genome contig
                    else:
                        # Do nothing
                        pass
        
        logging.info('datum_list:')
        for datum in datum_list:
            logging.info('  ' + str(datum))
        
        # Assuming we have 'group' with the [datum_after, datum_before]
        # Now we match homology regions with r0
        for dg in datum_groups:
            # group = [
            #     Datum(dDNA_r=1, dDNA_contig='exDonor-55', genome_r=1, genome_contig='Ca22chr3A_C_albicans_SC5314-r1[exDonor-55]', ush_start=943366, ush_end=None, dsh_start=None, dsh_end=943468)
            #     Datum(dDNA_r=2, dDNA_contig='reDonor-0', genome_r=1, genome_contig='Ca22chr3A_C_albicans_SC5314-r1[exDonor-55]', ush_start=942963, ush_end=943420, dsh_start=943414, dsh_end=943996)
            # ]
            ush_len = abs(dg[-1].ush_end - dg[-1].ush_start)
            dsh_len = abs(dg[-1].dsh_end - dg[-1].dsh_start)
            
            dg_ush_seq = genome_contigs_list[dg[-1].genome_r][dg[-1].genome_contig][dg[-1].ush_start:dg[-1].ush_start+ush_len]
            dg_dsh_seq = genome_contigs_list[dg[-1].genome_r][dg[-1].genome_contig][dg[-1].dsh_end-dsh_len:dg[-1].dsh_end]
            
            for group in groups:
                g_ush_seq = genome_contigs_list[group[0].genome_r][group[0].genome_contig][group[0].ush_start:group[0].ush_start+ush_len]
                g_dsh_seq = genome_contigs_list[group[0].genome_r][group[0].genome_contig][group[0].dsh_end-dsh_len:group[0].dsh_end]
                
                # If the rounds match, the dDNA_contigs match, and the genome_contigs match
                # and the homology regions match
                if ((dg[-1].dDNA_r == group[0].dDNA_r) and
                    (dg[-1].dDNA_contig == group[0].dDNA_contig) and
                    in_same_contig_group(group[0].genome_contig, dg[-1].genome_contig) and
                    (nucleotides.lcs(dg_ush_seq, g_ush_seq).size > 0.9*ush_len) and
                    (nucleotides.lcs(dg_dsh_seq, g_dsh_seq).size > 0.9*dsh_len)):
                    # Then populate their missing fields
                    #if (group[0].ush_end == None):
                    #    group[0].ush_end = group[-1].ush_start+ush_len
                    #if (group[0].dsh_start == None):
                    #    group[0].dsh_start = group[-1].dsh_end-dsh_len
                    
                    # Then add them to the same group
                    for g in group:
                        dg.append(g)
        
        logging.info('datum_groups:')
        for dg in datum_groups:
            logging.info('  ' + str(dg))
        
        # Populate the final/terminal group 'rN'
        # Using homology! (like we did with r0 above)
        for dg in datum_groups:
            # group = [
            #     Datum(dDNA_r=1, dDNA_contig='exDonor-55', genome_r=1, genome_contig='Ca22chr3A_C_albicans_SC5314-r1[exDonor-55]', ush_start=943366, ush_end=None, dsh_start=None, dsh_end=943468)
            #     Datum(dDNA_r=2, dDNA_contig='reDonor-0', genome_r=1, genome_contig='Ca22chr3A_C_albicans_SC5314-r1[exDonor-55]', ush_start=942963, ush_end=943420, dsh_start=943414, dsh_end=943996)
            # ]
            ush_len = abs(dg[-1].ush_end - dg[-1].ush_start)
            dsh_len = abs(dg[-1].dsh_end - dg[-1].dsh_start)
            
            dg_ush_seq = genome_contigs_list[dg[-1].genome_r][dg[-1].genome_contig][dg[-1].ush_start:dg[-1].ush_start+ush_len]
            dg_dsh_seq = genome_contigs_list[dg[-1].genome_r][dg[-1].genome_contig][dg[-1].dsh_end-dsh_len:dg[-1].dsh_end]
            
            to_pop = []
            for di, datum in enumerate(datum_list):
                d_ush_seq = genome_contigs_list[datum.genome_r][datum.genome_contig][datum.ush_start:datum.ush_start+ush_len]
                d_dsh_seq = genome_contigs_list[datum.genome_r][datum.genome_contig][datum.dsh_end-dsh_len:datum.dsh_end]
                
                # If the rounds match, the dDNA_contigs match, and the genome_contigs match
                # and the homology regions match
                if ((dg[-1].dDNA_r == datum.dDNA_r) and
                    (dg[-1].dDNA_contig == datum.dDNA_contig) and
                    in_same_contig_group(datum.genome_contig, dg[-1].genome_contig) and
                    #group[0].genome_contig.startswith(dg[-1].genome_contig) and
                    (nucleotides.lcs(dg_ush_seq, d_ush_seq).size > 0.9*ush_len) and
                    (nucleotides.lcs(dg_dsh_seq, d_dsh_seq).size > 0.9*dsh_len)):
                    
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
        logging.info('datum_groups (finished):')
        for dg in datum_groups:
            logging.info('  ' + str(dg))
        
        logging.info('datum_list (finished):')
        if (len(datum_list) == 0):
            logging.info('  EMPTY')
        else:
            for datum in datum_list:
                logging.info('  ' + str(datum))
        
        ######## End code for linking the loci ########
        
        # 'datum_groups' contains the linked loci
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
                        logging.info('too far')
                        us_done = True
                    fus = genome_contigs_list[datum.genome_r][datum.genome_contig][max(0, datum.ush_start-us_dist):datum.ush_start]
                    m = nucleotides.lcs(fus0, fus)
                    logging.info('us_dist = ' + str(us_dist))
                    logging.info('   fus0 = ' + fus0)
                    logging.info(('fus'+str(di+1)).rjust(7) +' = ' + fus)
                    logging.info('      m = ' + str(m))
                    if (m.size < 200):
                        us_dist += 100
                        break
                    else:
                        fus0 = fus0[m.a:m.a+m.size]
                else:
                    us_done = True
            
            logging.info('far_upstream_dist: ' + str(us_dist))
            logging.info('far_upstream_seq >= 200: ' + fus0)
            
            # Find the downstream LCS
            ds_dist = 600
            ds_done = False
            fds0 = None
            while (ds_done == False):
                fds0 = genome_contigs_list[dg[0].genome_r][dg[0].genome_contig][dg[0].dsh_end:dg[0].dsh_end+ds_dist]
                for di, datum in enumerate(dg[1:]):
                    if ((dg[0].dsh_end+ds_dist > len(genome_contigs_list[dg[0].genome_r][dg[0].genome_contig])) or (datum.dsh_end+ds_dist > len(genome_contigs_list[datum.genome_r][datum.genome_contig]))):
                        logging.info('too far')
                        ds_done = True
                    fds = genome_contigs_list[datum.genome_r][datum.genome_contig][datum.dsh_end:datum.dsh_end+ds_dist]
                    m = nucleotides.lcs(fds0, fds)
                    logging.info('us_dist = ' + str(ds_dist))
                    logging.info('   fus0 = ' + fds0)
                    logging.info(('fus'+str(di+1)).rjust(7) +' = ' + fds)
                    logging.info('      m = ' + str(m))
                    if (m.size < 200):
                        ds_dist += 100
                        break
                    else:
                        fds0 = fds0[m.a:m.a+m.size]
                else:
                    ds_done = True
            
            logging.info('far_downstream_dist: ' + str(ds_dist))
            logging.info('far_downstream_seq >= 200: ' + fds0)
            
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
            
            logging.info("checking 'sF' region:")
            logging.info("              fus0: " + fus0)
            for prpi2, prp2 in enumerate(pcr_region_positions[-1]):
                logging.info(str(prpi2) + ' ' + str(prp2[:2]).rjust(16) + ": " + genome_contigs_list[dg[prpi2].genome_r][dg[prpi2].genome_contig][prp2[0]:prp2[1]] + ' ' + str(dg[prpi2]))
            logging.info("checking 'sR' region:")
            logging.info("              fds0: " + fds0)
            for prpi2, prp2 in enumerate(pcr_region_positions[-1]):
                logging.info(str(prpi2) + ' ' + str(prp2[2:]).rjust(16) + ": " + genome_contigs_list[dg[prpi2].genome_r][dg[prpi2].genome_contig][prp2[2]:prp2[3]] + ' ' + str(dg[prpi2]))
        
        # Identify the feature/insert sequence where the 'rN-oF', 'rN-oR', 'rN-iF', 'rN-iR' primers should be located
        # and store them in 'insert_seqs'
        Insert = namedtuple('Insert', ['genome_r', 'genome_contig', 'seq', 'us_seq', 'ds_seq', 'fus_dist', 'fds_dist', 'type'])
        max_primer_length = 35
        
        def round_counter():
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
        
        for i, dg in enumerate(datum_groups):
            fus_seq, fds_seq = pcr_regions[i]
            
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
                        #datum.dDNA_r, # genome_r <----------------------- may need to change this to be the genome (not the dDNA)
                        #datum.dDNA_contig, # genome_contig <------------- may need to change this to be the genome (not the dDNA)
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
            
            #                                             shared_forward  shared_reverse  features/inserts
            pair_list, insert_pair_list, round_labels = self.calculate_them_primers(args, fus_seq,        fds_seq,        insert_seqs)
            
            # Filter 'pair_list' to get the top 10
            pair_list = sorted(pair_list, reverse=True)[:args.max_number_designs_reported]
            
            # Print the calculated weights of the sF/sR/oF/oR sets
            for ppli, (w, pp_list) in enumerate(pair_list):
                output0.append([ppli, i, w])
            
            # Print the primers for table output 1
            for ppli, (w, pp_list) in enumerate(pair_list):
                round_labels_iter = iter(round_labels)
                
                for pp_i, pp in enumerate(pp_list):
                    amp_name, round_n = next(round_labels_iter)
                    #amp_name = '-'
                    #round_n = '-'
                    
                    ra_counter = round_counter()
                    for ia in range(len(insert_pair_list)):
                        rac = next(ra_counter)
                        if (pp and (round_n == rac)):
                            output1.append([ppli, pp_i, i, amp_name, round_n, pp.forward_primer.name, pp.reverse_primer.name, pp.get_amplicon_size(), pp.get_tms()])
                        else:
                            if (amp_name == 'A'):
                                output1.append([ppli, pp_i, i, amp_name, rac, 'sF', 'sR', '-', '-'])
                            elif (amp_name == 'B'):
                                output1.append([ppli, pp_i, i, amp_name, rac, 'sF', 'r'+round_n+'-oR', '-', '-'])
                            elif (amp_name == 'C'):
                                output1.append([ppli, pp_i, i, amp_name, rac, 'r'+round_n+'-oF', 'sR', '-', '-'])
                    
                    #if pp: # primer_set_index, primer_pair_index, locus
                    #    #if (pp.forward_primer.name == 'sF'):
                    #    #    if (pp.reverse_primer.name == 'sR'):
                    #    #        amp_name = 'A'
                    #    #        round_n = next(A_counter)
                    #    #    else:
                    #    #        amp_name = 'B'
                    #    #        round_n = next(B_counter)
                    #    #else:
                    #    #    if (pp.reverse_primer.name == 'sR'):
                    #    #        amp_name = 'C'
                    #    #        round_n = next(C_counter)
                    #    #    else:
                    #    #        amp_name = 'D'
                    #    output1.append([ppli, pp_i, i, amp_name, round_n, pp.forward_primer.name, pp.reverse_primer.name, pp.get_amplicon_size(), pp.get_tms(), 'Template'])
                    #else:
                    #    output1.append([ppli, pp_i, i, amp_name, round_n, '-', '-', '-', '-', '-'])
            
            # Filter 'insert_pair_list' to get the top 10 (They are already sorted)
            insert_pair_list = [x[:args.max_number_designs_reported] for x in insert_pair_list]
            #D_counter = round_counter() # Hacky way to get the round
            #for ip_r, iF_iR_paired_primers in enumerate(insert_pair_list):
            #    amp_name = 'D'
            #    round_n = next(D_counter)
            #    if (len(iF_iR_paired_primers) == 0):
            #        output1.append(['-', ip_r, i, amp_name, round_n, '-', '-', '-', '-', '-'])
            #    else:
            #        for ppli, pp in enumerate(iF_iR_paired_primers):
            #            if pp:
            #                output1.append([ppli, ip_r, i, amp_name, round_n, pp.forward_primer.name, pp.reverse_primer.name, pp.get_amplicon_size(), pp.get_tms(), 'Template'])
            #            else:
            #                output1.append([ppli, ip_r, i, amp_name, round_n, '-', '-', '-', '-', '-'])
            
            if (len(insert_pair_list) > 0):
                number_D_sets = max([len(ipl) for ipl in insert_pair_list])
            else:
                number_D_sets = 0
            
            #D_sets = [] # [['0b', '0a', '1b', '1a'], ...]
            #for ppli in range(number_D_sets):
            #    amp_name = 'D'
            #    D_counter = round_counter() # Hacky way to get the round
            #    D_sets.append([])
            #    for ip_r, iF_iR_paired_primers in enumerate(insert_pair_list):
            #        round_n = next(D_counter)
            #        
            #        try:
            #            pp = iF_iR_paired_primers[ppli]
            #            D_sets[-1].append([ppli, ip_r, i, amp_name, round_n, pp.forward_primer.name, pp.reverse_primer.name, pp.get_amplicon_size(), pp.get_tms()])
            #        except IndexError, AttributeError:
            #            D_sets[-1].append([ppli, ip_r, i, amp_name, round_n, round_n+'-iF', round_n+'-iR', '-', '-'])
            #
            #for D_set in D_sets:
            #    for ds in D_set:
            #        output1.append(ds)
            
            amp_name = 'D'
            for ppli in range(number_D_sets):
                D_counter = round_counter() # Hacky way to get the round
                for ip_r, iF_iR_paired_primers in enumerate(insert_pair_list):
                    round_n = next(D_counter)
                    try:
                        pp = iF_iR_paired_primers[ppli]
                    except IndexError:
                        pp = None
                    ra_counter = round_counter()
                    for ia in range(len(insert_pair_list)):
                        rac = next(ra_counter)
                        if (pp and (round_n == rac)):
                            output1.append([ppli, ip_r, i, amp_name, round_n, pp.forward_primer.name, pp.reverse_primer.name, pp.get_amplicon_size(), pp.get_tms()])
                        else:
                            output1.append([ppli, ip_r, i, amp_name, rac, 'r'+round_n+'-iF', 'r'+round_n+'-iR', '-', '-'])
                    
            
            
            # Print the primers for table output 2
            for ppli, (w, pp_list) in enumerate(pair_list):
                for pp_i, pp in enumerate(pp_list):
                    if pp:
                        # 'f_locations' and 'r_locations' should have the same length
                        f_locations = self.get_primer_location(genome_fasta_file_list, genome_contigs_list, pp.forward_primer.sequence)
                        r_locations = self.get_primer_location(genome_fasta_file_list, genome_contigs_list, pp.reverse_primer.sequence)
                        for loc_i in range(len(f_locations)):
                            for loc in f_locations[loc_i]:
                                output2.append([ppli, pp_i, i, pp.forward_primer.name, pp.forward_primer.sequence, loc[0], loc[1], loc[2], loc[3], loc[4]])
                            for loc in r_locations[loc_i]:
                                output2.append([ppli, pp_i, i, pp.reverse_primer.name, pp.reverse_primer.sequence, loc[0], loc[1], loc[2], loc[3], loc[4]])
                        #output2.append([ppli, pp_i, i, pp.reverse_primer.name, pp.reverse_primer.sequence, 'File', 'Contig', 'Start', 'End', '-'])
                    else:
                        output2.append([ppli, pp_i, i, '-', '-', '-', '-', '-', '-', '+'])
                        output2.append([ppli, pp_i, i, '-', '-', '-', '-', '-', '-', '-'])
            
            for ip_r, iF_iR_paired_primers in enumerate(insert_pair_list):
                for ppli, pp in enumerate(iF_iR_paired_primers):
                    if pp:
                        # 'f_locations' and 'r_locations' should have the same length
                        f_locations = self.get_primer_location(genome_fasta_file_list, genome_contigs_list, pp.forward_primer.sequence)
                        r_locations = self.get_primer_location(genome_fasta_file_list, genome_contigs_list, pp.reverse_primer.sequence)
                        for loc_i in range(len(f_locations)):
                            for loc in f_locations[loc_i]:
                                output2.append([ppli, ip_r, i, pp.forward_primer.name, pp.forward_primer.sequence, loc[0], loc[1], loc[2], loc[3], loc[4]])
                            for loc in r_locations[loc_i]:
                                output2.append([ppli, ip_r, i, pp.reverse_primer.name, pp.reverse_primer.sequence, loc[0], loc[1], loc[2], loc[3], loc[4]])
                    else:
                        output2.append([ppli, ip_r, i, '-', '-', '-', '-', '-', '-', '+'])
                        output2.append([ppli, ip_r, i, '-', '-', '-', '-', '-', '-', '-'])
        
        # Output header information for first output table
        print('#                                          Genome')
        print('# Amplicon  ──upstream─┐┌─homology─┐┌──insert/feature──┐┌─homology─┐┌─downstream──')
        print('#        A   sF ===>····················································<=== sR')
        print('#        B   sF ===>···················<=== rN-oR')
        print('#        C                             rN-oF ===>·······················<=== sR')
        print('#        D                        rN-iF ===>······<=== rN-iR')
        print('#')
        print('# F=forward, R=reverse')
        print('# s=shared, o=outer, i=inner')
        print('# rN=round number')
        print('#')
        
        print('# '+ '\t'.join(['Set', 'Locus', 'Weight']))
        for line in output0:
            print('\t'.join(map(str, line)))
        
        print('# '+ '\t'.join(['Set', 'Index', 'Locus', 'Amplicon', 'Template', 'F', 'R', 'Size', 'Tm']))
        for line in output1:
            print('\t'.join(map(str, line)))
        
        # Output header information for second output table
        print('# ' + '\t'.join(['Set', 'Index', 'Locus', 'Primer', 'Sequence', 'File', 'Contig', 'Start', 'End', 'Strand']))
        for line in output2:
            print('\t'.join(map(str, line)))
        
                        #####################
                        # Cut code was here #
                        #####################
                        
                        # This code output GenBank '*.gb' files
        
        # End 'compute()'

    def calculate_them_primers(self, args, far_us_seq, far_ds_seq, insert_seqs):
        """
        """
        tm_max_difference = 4.0
        amplicon_size = (200, 900)
        tm_range = (52, 64)
        primer_length_range = (19,32)
        min_delta_g = -5.0
        
        subset_size = 1000 #35 #200 # Temporary size limit for the number of primers that go into the pair() function
        
        logging.info('Scanning the regions (shared upstream (sF), shared downstream (sR), feature/insert (rN-oF,rN-oR,rN-iF,rN-iR) for all decent primers')
        
        # Make the 'temp_folder' path
        temp_folder = os.path.join('/dev/shm/addtag', os.path.basename(args.folder))
        
        logging.info("Scanning far upstream for 'sF' primers:")
        sF_list = sorted(
            args.selected_oligo.scan(far_us_seq, 'left', primer_size=primer_length_range, tm_range=tm_range, min_delta_g=min_delta_g, folder=temp_folder, time_limit=args.primer_scan_limit),
            key=lambda x: x.weight,
            reverse=True
        )
        logging.info('  len(sF_list) = {}'.format(len(sF_list)))
        logging.info('  sF: skipping {}/{} calculated primers'.format(max(0, len(sF_list)-subset_size), len(sF_list)))
        sF_list = sF_list[:subset_size]
        # Need to add a condition that if (len(sF_list) < N), then it should re-evaluate the far-upstream region to try to get
        # more sequence to search.
        
        logging.info("Scanning far downstream for 'sR' primers:")
        sR_list = sorted(
            args.selected_oligo.scan(far_ds_seq, 'right', primer_size=primer_length_range, tm_range=tm_range, min_delta_g=min_delta_g, folder=temp_folder, time_limit=args.primer_scan_limit),
            key=lambda x: x.weight,
            reverse=True
        )
        logging.info('  len(sR_list) = {}'.format(len(sR_list)))
        logging.info('  sR: skipping {}/{} calculated primers'.format(max(0, len(sR_list)-subset_size), len(sR_list)))
        sR_list = sR_list[:subset_size]
        # Need to add condition that if there aren't enough putative 'sR' primers, then the sR region is expanded
        # Otherwise, the script can just end prematurely
        
        insert_list = []
        for ins in insert_seqs:
            logging.info("Scanning feature/insert for 'rN-oF', 'rN-oR', 'rN-iF', 'rN-iR' primers:")
            iF_list = sorted(
                args.selected_oligo.scan(ins.seq, 'left',  primer_size=primer_length_range, tm_range=tm_range, min_delta_g=min_delta_g, folder=temp_folder, us_seq=ins.us_seq, ds_seq=ins.ds_seq, time_limit=args.primer_scan_limit),
                key=lambda x: x.weight,
                reverse=True
            )
            logging.info('  len(iF_list) = {}'.format(len(iF_list)))
            logging.info('  insert_F: skipping {}/{} calculated primers'.format(max(0, len(iF_list)-subset_size), len(iF_list)))
            iF_list = iF_list[:subset_size]
            
            iR_list = sorted(
                args.selected_oligo.scan(ins.seq, 'right', primer_size=primer_length_range, tm_range=tm_range, min_delta_g=min_delta_g, folder=temp_folder, us_seq=ins.us_seq, ds_seq=ins.ds_seq, time_limit=args.primer_scan_limit),
                key=lambda x: x.weight,
                reverse=True
            )
            logging.info('  len(iR_list) = {}'.format(len(iR_list)))
            logging.info('  insert_R: skipping {}/{} calculated primers'.format(max(0, len(iR_list)-subset_size), len(iR_list)))
            iR_list = iR_list[:subset_size]
            
            insert_list.append([iF_list, iR_list])
        
        # Select the best subset of primers so the pairing doesn't take too long
        # This will cap the maximum number of primer pairs for each region span
        # to be 100 x 100 = 10,000
        
        # Do the pair calculations that involve 'sF' and 'sR'
        logging.info("Calculating: 'sF' 'sR' paired primers...")
        sF_sR_paired_primers = args.selected_oligo.pair(sF_list, sR_list, amplicon_size=(0, amplicon_size[1]), tm_max_difference=tm_max_difference, intervening=0, min_delta_g=min_delta_g, folder=temp_folder, time_limit=args.primer_pair_limit)
        
        for pp in sF_sR_paired_primers:
            # Re-calculate the PrimerPair weights to prefer the smallest amplicon sizes
            pp.weight = pp.get_weight(minimize=True)
            
            # Re-name the primers so when they are printed, they are easy to distinguish
            pp.forward_primer.name = 'sF'
            pp.reverse_primer.name = 'sR'
        
        # Sort by weight
        sF_sR_paired_primers = sorted(
            sF_sR_paired_primers,
            key=lambda x: x.get_joint_weight(),
            reverse=True
        )
        logging.info('  len(sF_sR_paired_primers) = {}'.format(len(sF_sR_paired_primers)))
        
        # pair_list[ins index] = [sF_oR_paired_primers, oF_sR_paired_primers, iF_iR_paired_primers]
        pair_list = []
        insert_pair_list = []
        
        #for i, (iF_list, iR_list) in enumerate(insert_list):
        for i in range(len(insert_seqs)):
            iF_list, iR_list = insert_list[i]
            ins = insert_seqs[i]
            
            logging.info('Pairing primers: ' + str(i))
        
            # If there is a hard constraint for primers that should be used
            # That is, if flanktag primers should be used
            if (len(args.primers) > 0):
                pass
            # Otherwise, no flanktags are specified
            else:
                # sF rN-oR
                logging.info("  Calculating: 'sF' 'r"+str(ins.genome_r)+ins.type+"-oR' paired_primers...")
                sF_oR_paired_primers = sorted(
                    args.selected_oligo.pair(sF_list, iR_list, amplicon_size=amplicon_size, tm_max_difference=tm_max_difference, intervening=ins.fus_dist, min_delta_g=min_delta_g, folder=temp_folder, time_limit=args.primer_pair_limit),
                    key=lambda x: x.get_joint_weight(),
                    reverse=True
                )
                for pp in sF_oR_paired_primers:
                    pp.forward_primer.name, pp.reverse_primer.name = 'sF', 'r'+str(ins.genome_r)+ins.type+'-oR'
                logging.info('  len(sF_oR_paired_primers) = {}'.format(len(sF_oR_paired_primers)))
            
                # rN-0F sR
                logging.info("  Calculating: 'r"+str(ins.genome_r)+ins.type+"-oF' 'sR' paired_primers...")
                oF_sR_paired_primers = sorted(
                    args.selected_oligo.pair(iF_list, sR_list, amplicon_size=amplicon_size, tm_max_difference=tm_max_difference, intervening=ins.fds_dist, min_delta_g=min_delta_g, folder=temp_folder, time_limit=args.primer_pair_limit),
                    key=lambda x: x.get_joint_weight(),
                    reverse=True
                )
                for pp in oF_sR_paired_primers:
                    pp.forward_primer.name, pp.reverse_primer.name = 'r'+str(ins.genome_r)+ins.type+'-oF', 'sR'
                logging.info('  len(oF_sR_paired_primers) = {}'.format(len(oF_sR_paired_primers)))
                
                # rN-iF rN-iR
                logging.info("  Calculating: 'r"+str(ins.genome_r)+ins.type+"-iF' 'r"+str(ins.genome_r)+ins.type+"-iR' paired_primers...")
                iF_iR_paired_primers = sorted(
                    args.selected_oligo.pair(iF_list, iR_list, amplicon_size=amplicon_size, tm_max_difference=tm_max_difference, intervening=0, same_template=True, min_delta_g=min_delta_g, folder=temp_folder, time_limit=args.primer_pair_limit),
                    key=lambda x: x.get_joint_weight(),
                    reverse=True
                )
                for pp in iF_iR_paired_primers:
                    pp.forward_primer.name, pp.reverse_primer.name = 'r'+str(ins.genome_r)+ins.type+'-iF', 'r'+str(ins.genome_r)+ins.type+'-iR'
                logging.info('  len(iF_iR_paired_primers) = {}'.format(len(iF_iR_paired_primers)))
                
                # Add calculated primer pairs to the list
                pair_list.append([sF_oR_paired_primers, oF_sR_paired_primers, ins])
                insert_pair_list.append(iF_iR_paired_primers)
        
        # If a primer in feature/insert of one round ALSO is identical
        # with a primer in another round,
        # Then it should only be weighted for a SINGLE entry
        starting_set, finished_set, round_labels = self.calculate_them_best_set(args, sF_sR_paired_primers, pair_list)
        
        return finished_set, insert_pair_list, round_labels
    
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
    
    def calculate_them_best_set(self, args, sF_sR_paired_primers, pair_list, max_iterations=1000):
        """
        Takes input primer sets and calculates their weights
        
        This starts with the best from a sorted list of 'sF' 'sR' primer pairs
        Then it optimizes by swapping the worst component a finite number of times
        (or until the delta drops below a certain amount)
        """
        # Define helper functions specific to only this method
        def nsum(n):
            return n*(n+1)/2
        
        def rank_order(x, reverse=False, shift=0):
            return [y+shift for y in sorted(range(len(x)), key=x.__getitem__, reverse=reverse)]
        
        def rank(x, reverse=False, shift=0):
            return [z[0]+shift for z in sorted(enumerate(sorted(enumerate(x), key=lambda w: w[1], reverse=reverse)), key=lambda y: y[1][0])]
        
        def rank_dist(x):
            s = nsum(len(x))
            return [y/s for y in rank(x, reverse=False, shift=1)]
        
        def product(x):
            z = 1
            for y in x:
                z *= y
            return z
        
        def random_choices(population, weights, k=1):
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
        
        def filter_primer_pairs(pairs, forward=None, reverse=None):
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
        
        def random_primer_by_weight_old(pairs, forward=None, reverse=None):
            ooo = filter_primer_pairs(pairs, forward, reverse)
            
            if (len(ooo) == 0):
                return None
            elif ('choices' in random.__all__):
                return random.choices(ooo, [x.get_joint_weight() for x in ooo])[0]
            else:
                return random_choices(ooo, [x.get_joint_weight() for x in ooo])[0]
        
        def random_primer_by_weight(pairs):
            if (len(pairs) == 0):
                return None
            elif ('choices' in random.__all__):
                return random.choices(pairs, [x.get_joint_weight() for x in pairs])[0]
            else:
                return random_choices(pairs, [x.get_joint_weight() for x in pairs])[0]
        
        # This should be in a loop...
        #sF_oR_paired_primers, oF_sR_paired_primers, iF_iR_paired_primers, ins = pair_list[0]
        
        # Set initial iteration count to zero
        iteration_count = 0
        
        # Create lists to hold 5-set of primer pairs, and their joint weight
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
            logging.info('sF_sR_paired_primers: skipping {}/{} calculated primer pairs'.format(max(0, len(sF_sR_paired_primers)-max_iterations), len(sF_sR_paired_primers)))
        
        while (iteration_count < max_iterations):
            # Increment the count
            iteration_count += 1
            
            # Create random 6-sets, with each pair's probability determined by its joint weight
            #uf_dr_pair = random_primer_by_weight(sF_sR_paired_primers) # comment this out to prevent random sampling of this one
            try:
                sF_sR_pair = next(sF_sR_iterator)
                #sF_sR_pair.forward_primer.name = 'sF'
                #sF_sR_pair.reverse_primer.name = 'sR'
            except StopIteration:
                break
            
            logging.info('loop {}:'.format(iteration_count))
            
            # Populate set with all primer pairs that have 'sF' and 'sR'
            pp_sources = []
            pp_setN = []
            p_setN = []
            for sF_oR_paired_primers, oF_sR_paired_primers, ins in pair_list:
                # Add the 'oR' from 'sF-oR' pair
                pp_sources.append(filter_primer_pairs(sF_oR_paired_primers, forward=sF_sR_pair.forward_primer))
                #for pp in pp_sources[-1]:
                #    if pp:
                #        pp.forward_primer.name, pp.reverse_primer.name = 'sF', 'r'+str(ins.genome_r)+ins.type+'-oR'
                pp = random_primer_by_weight(pp_sources[-1]) # Pick a random primer by weight
                pp_setN.append(pp) # Will be 'None' if 'pp_sources[-1]' is empty
                if pp:
                    p_setN.append(pp.reverse_primer)
                else:
                    p_setN.append(None)
                
                # Add the 'oF' from 'oF-sR' pair
                pp_sources.append(filter_primer_pairs(oF_sR_paired_primers, reverse=sF_sR_pair.reverse_primer))
                #for pp in pp_sources[-1]:
                #    if pp:
                #        pp.forward_primer.name, pp.reverse_primer.name = 'r'+str(ins.genome_r)+ins.type+'-oF', 'sR'
                pp = random_primer_by_weight(pp_sources[-1]) # Pick a random primer by weight
                pp_setN.append(pp) # Will be 'None' if 'pp_sources[-1]' is empty
                if pp:
                    p_setN.append(pp.forward_primer)
                else:
                    p_setN.append(None)
            
            logging.info('  length of sources: ' + str([len(x) for x in pp_sources]))
            
            # # Pick random primer by weight for each source
            # pp_setN = []
            # p_setN = []
            # for pp_list in pp_sources:
            #     pp = random_primer_by_weight(pp_list)
            #     pp_setN.append(pp) # Will be 'None' if 's' is empty
            # 
            # # Break down the 'pp_setN' which contains primer pairs into a set of primers
            # p_setN = []
            # for pp in pp_setN:
            #     if pp:
            #         p_setN.append(pp.reverse_primer)
            #         p_setN.append(pp.forward_primer)
            #     else:
            #         p_setN.append(None)
            
            group_weight = args.selected_oligo.group_weight([sF_sR_pair.forward_primer, sF_sR_pair.reverse_primer]+p_setN)
            joint_weight = group_weight * product(x.get_joint_weight() for x in [sF_sR_pair]+pp_setN if x)
            starting_set.append((joint_weight, [sF_sR_pair]+pp_setN))
            
            logging.info('  starting_set: ' + str(starting_set[-1]))
            
            # iteratively improve until a local maxima is found
            # By continually swapping out the least-weighted component primer
            min_weight_delta = 1e-20
            weight_delta = 1
            while(weight_delta > min_weight_delta):
                # Get list of indices from smallest weight to largest weight (excluding 'sF' and 'sR')
                wi_order = rank_order([x.weight if (x != None) else math.inf for x in p_setN])
                
                # If primer is 'None', then its weight will be 'math.inf', and the index
                # corresponding to it will be the last element in wi_order:
                #   rank_order([1, 2, 2.5, 0.001, math.inf]) # [3, 0, 1, 2, 4]
                
                # Start at the worst, and go to the next-worst, then next, then next
                for wi in wi_order:
                    # Copy the list of primer pairs
                    pp_setNN = pp_setN[:]
                    
                    # Swap the worst-performing primer with an alternative.
                    # If the alternative gives a better weight, then keep it
                    # and break out of the loop
                    for source_pp in pp_sources[wi]:
                        pp_setNN[wi] = source_pp
                        
                        p_setNN = []
                        for pp in pp_setNN:
                            if pp:
                                if (wi % 2 == 0):
                                    p_setNN.append(pp.reverse_primer)
                                else:
                                    p_setNN.append(pp.forward_primer)
                            else:
                                p_setNN.append(None)
                        
                        new_group_weight = args.selected_oligo.group_weight([sF_sR_pair.forward_primer, sF_sR_pair.reverse_primer]+p_setNN)
                        new_joint_weight = new_group_weight * product(x.get_joint_weight() for x in [sF_sR_pair]+pp_setNN if x)
                        
                        if (new_joint_weight > joint_weight):
                            logging.info('           set: ' + str((new_joint_weight, [sF_sR_pair]+pp_setNN)))
                            
                            weight_delta = new_joint_weight - joint_weight
                            
                            # replace values for next iteration
                            pp_setN = pp_setNN
                            joint_weight = new_joint_weight
                            break
                    else:
                        weight_delta = 0
                    
                    if (weight_delta > 0):
                        break
            
            # Does not re-calculate 'pp.weight', thus the amplicon_size will not be further penalized
            # However, the 'pp.get_amplicon_size()' will accurately reflect the updated 'pp.intervening' length
            #sF_insert_sR_pair = copy.deepcopy(sF_sR_pair)
            #sF_insert_sR_pair.intervening = len(q_hih_seq)
            #sF_feature_sR_pair = copy.deepcopy(sF_sR_pair)
            #sF_feature_sR_pair.intervening = len(s_ush_seq)+len(sseq_feature)+len(s_dsh_seq)
            
            # Duplicate 'sF_sR_pair' for each template, giving it the proper 'intervening' length (amplicon size)
            #sF_sR_pair_list = []
            #for sF_oR_paired_primers, oF_sR_paired_primers, ins in pair_list:
            #    sF_sR_pair_list.append(copy.deepcopy(sF_sR_pair))
            #    sF_sR_pair_list[-1].intervening = ins.fus_dist + len(ins.seq) + ins.fds_dist
            
            finished_set.append((joint_weight, [sF_sR_pair]+pp_setN)) # 1 copy of 'sF_sR_pair'
            #finished_set.append((joint_weight, sF_sR_pair_list+pp_setN)) # 4 copies of 'sF_sR_pair': ('0b', '0a', '1b', '1a')
        
        return starting_set, finished_set, round_labels
    
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
        with open(alignment_filename, 'r') as flo:
            record = None
            while True:
                record = args.selected_aligner.load_record(flo)
                if (record == None):
                    break
                else:
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
                                logging.info(str(qs_pair) + ' has too many valid alignments')
                                invalid_records.add(qs_pair)
                        
                        elif (record.query_position[1] > len(dDNA_contigs[record.query_name])-permitted_edge): # right
                            if (records[qs_pair][1] == None):
                                records[qs_pair][1] = record
                            else:
                                logging.info(str(qs_pair) + ' has too many valid alignments')
                                invalid_records.add(qs_pair)
        
        logging.info("before filtering:")
        for kkk, vvv in records.items():
            logging.info("  "+str(kkk) + " " + str(vvv))
        
        # Queue for removal the (query, subject) pairs that only have 1 of the 2 required alignments
        for qs_pair, record_pair in records.items():
            if (None in record_pair):
                invalid_records.add(qs_pair)
        
        # Remove the (query, subject) pairs that had too many valid alignments
        for bad_key in invalid_records:
            records.pop(bad_key)
            logging.info('Removing invalid record: ' + str(bad_key))
        
        logging.info("after filtering:")
        for kkk, vvv in records.items():
            logging.info("  "+str(kkk) + " " + str(vvv))
        
        return records
    
    def filter_alignment_records2(self, args, alignment_filename, dDNA_contigs):
        records = {}
        
        # Allow for alignments to omit up to 3 nt at edges (of homology regions)
        permitted_edge = 3
        
        # Make a set of invalid records that will be removed because they failed some quality threshold
        invalid_records = set()
        
        # Read the alignment file
        with open(alignment_filename, 'r') as flo:
            record = None
            while True:
                record = args.selected_aligner.load_record(flo)
                if (record == None):
                    break
                else:
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
                                logging.info(str(qs_pair) + ' has too many valid alignments')
                                invalid_records.add(qs_pair)
        for bad_qs_pair in invalid_records:
            records.pop(bad_qs_pair)
            logging.info('Removing invalid record: ' + str(bad_qs_pair))
        
        return records