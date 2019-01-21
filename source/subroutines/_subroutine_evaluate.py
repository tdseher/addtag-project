#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/subroutines/_subroutine_evaluate.py

# Import included AddTag-specific modules
from . import subroutine
from .. import utils
from ..feature import Feature
from ..targets import Target, ExcisionTarget, ReversionTarget
from ..donors import Donor, ExcisionDonor, ReversionDonor


class EvaluateParser(subroutine.Subroutine):
    def __init__(self, subparsers):
        self.subparsers = subparsers
        
        self.name = 'evaluate'
        self.description = (
            "description:" "\n"
            "  Evaluate pre-designed CRISPR/Cas oligonucleotide sequences." "\n"
        )
        self.help = "Evaluate pre-designed CRISPR/Cas oligonucleotide sequences."
        self.epilog = (
            "example:" "\n"
            "  No example (yet)." "\n"
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
        required_group.add_argument("--gff", required=True, metavar="*.gff", type=str,
            help="GFF file specifying chromosomal features that should be \
                 multiplexed together.")
        required_group.add_argument("--folder", required=True, metavar="FOLDER",
            type=str, help="Path of folder to store generated files.")
        
        # Add required, mutually-exclusive group
        me_group = self.parser.add_mutually_exclusive_group(required=True)
        me_group.add_argument("--spacers", metavar="*.fasta",
            help="Evaluate where spacers would target, and for each site, what \
                 its scores would be.")
        
        me_group.add_argument("--gRNAs", metavar="*.fasta",
            help="Evaluate where gRNAs would target, and for each site, what \
                 its scores would be.")
        
        me_group.add_argument("--dDNAs", nargs="+", metavar="*.fasta", type=str,
            help="Homology arms, presence of unique gRNA target site, and which \
                 features these dDNAs would disrupt, and whether-or-not the \
                 dDNA would introduce a mutation.")
        
        me_group.add_argument("--primers", metavar="*.fasta",
            help="Evaluate what amplicon sizes to expect for wt/ko/ki for \
                 input primers. Also evaluate the Tm/sensitivity/dimerization/etc \
                 for the input primers. Does not evaluate these in multiplex.")
        
        # Maybe: make --gRNAs, --dDNAs, and --primers all accept the same number of FASTA files
        #   Then, when it does calculations, it does all of them in serial steps.
        #  step 1: --gRNAs[0], --dDNAs[0], --primers[0]
        #  step 2: --gRNAs[1], --dDNAs[1], --primers[1]
        #  etc...
        
        # Add optional arguments
        
        # If evaluating knock-in gRNA, then include the ko-dDNA with the "--fasta" argument??
        #  addtag evaluate
        #   --ko-gRNA            
        #   --ko-dDNA            Evaluate knock-out dDNA for which DNA features it would knock-out, and evaluate whether it has a unique gRNA target
        #   --ki-gRNA            Evaulate if knock-in gRNA would target the input ko-dDNA
        #   --ki-dDNA            Evaluate if knock-in dDNA will restore wild type correctly, or where the mutations will be
        #   --primers FASTA      Evaluate what amplicon sizes to expect for wt/ko/ki for input primers.
        #                        Also evaluate the Tms/sensitivity/dimerization/etc for the input primers
        #                        (DOES NOT EVALUATE THESE IN MULTIPLEX)
        #                        >feature1 F
        #                        NNNNNNNNNNNNNNNNNNNN
        #                        >feature1 R
        #                        NNNNNNNNNNNNNNNNNNNN
        
        self.parser.add_argument("--motifs", metavar="MOTIF", nargs="+", type=str,
            default=["N{20}>NGG"],
            help="Find only targets with these 'SPACER>PAM' motifs, written from \
            5' to 3'. '>' points toward PAM. IUPAC ambiguities accepted. '{a,b}' \
            are quantifiers. '/' is a sense strand cut, '\\' is an antisense strand \
            cut, and '|' is a double-strand cut. '.' is a base used for positional \
            information, but not enzymatic recognition. Be sure to enclose each \
            motif in quotes so your shell does not interpret STDIN/STDOUT redirection.")
        self.parser.add_argument("--off_target_motifs", metavar="MOTIF", nargs="+", type=str,
            default=[],
            help="Defaults to the same as the on-target motif. Definition syntax is identical.")
        
    def compute(self, args):
        '''UNDER DEVELOPMENT'''
        print("Perform the evaluation here.")
        
        contig_sequences = utils.load_multiple_fasta_files(args.fasta)
        
        dDNA_file = ExcisionDonor.generate_fasta(os.path.join(args.folder, 'excision-dDNAs.fasta'))
        genome_fasta_file = utils.write_merged_fasta(contig_sequences, os.path.join(args.folder, 'genome.fasta'))
        
        # search spacer FASTA against genome+dDNA FASTA
        genome_index_file = args.selected_aligner.index(genome_fasta_file, os.path.basename(genome_fasta_file), args.folder, args.processors)
        dDNA_index_file = args.selected_aligner.index(dDNA_file, os.path.basename(dDNA_file), args.folder, args.processors)
        
        q2gDNA_align_file = args.selected_aligner.align(ex_query_file, genome_index_file, 'excision-query-2-gDNA.'+args.selected_aligner.output, args.folder, args.processors)
        q2exdDNA_align_file = args.selected_aligner.align(ex_query_file, ex_dDNA_index_file, 'excision-query-2-excision-dDNA.'+args.selected_aligner.output, args.folder, args.processors)
        
        ExcisionTarget.load_alignment(q2gDNA_align_file, args, contig_sequences)
        ExcisionTarget.load_alignment(q2exdDNA_align_file, args, ExcisionDonor.get_contig_dict())
        
        # Calculate off-target/guide scores for each algorithm
        logging.info("ExcisionTarget after SAM parsing and off-target scoring")
        ##### Add short-circuit/heuristic #####
        for et_seq, et_obj in ExcisionTarget.sequences.items():
            et_obj.score_off_targets(args, homologs)
            logging.info(et_obj)
            for a in et_obj.alignments:
                logging.info('  ' + str(a))
        
        # Batch calculate with new ExcisionTarget class
        ExcisionTarget.score_batch()
        
        logging.info("ExcisionTarget after Azimuth calculation")
        for et_seq, et_obj in ExcisionTarget.sequences.items():
            logging.info(et_obj)
        
        # End 'compute()'
        