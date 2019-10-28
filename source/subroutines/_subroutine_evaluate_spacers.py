#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/subroutines/_subroutine_evaluate_spacers.py

# Import standard packages
import os
import logging

# Import included AddTag-specific modules
from . import subroutine
from .. import utils
from .. import algorithms
from .. import aligners
from ..feature import Feature
from ..targets import ExcisionTarget

logger = logging.getLogger(__name__)

class EvaluateSpacersParser(subroutine.Subroutine):
    logger = logger.getChild('EvaluateParser')
    
    def __init__(self, subparsers):
        self.subparsers = subparsers
        
        self.name = 'evaluate_spacers'
        self.description = (
            "description:" "\n"
            "  Evaluate pre-designed CRISPR/Cas oligonucleotide sequences." "\n"
            "  Calculates the score for each Algorithm, and a final weight for each input Spacer" "\n"
        )
        self.help = "Evaluate scores of pre-designed spacers."
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

        required_group.add_argument("--gff", required=True, nargs="+", metavar="*.gff", type=str,
            help="GFF files specifying chromosomal features that should be \
            targetd (in multiplex) or excluded.")

        required_group.add_argument("--folder", required=True, metavar="FOLDER",
            type=str, help="Path of folder to store generated files.")

        required_group.add_argument("--spacers", required=True, nargs="+", metavar="*.fasta", type=str,
            help="Evaluate where spacers would target, and for each site, what \
            its scores would be.")

        # Add required, mutually-exclusive group
        # me_group = self.parser.add_mutually_exclusive_group(required=True)
        #
        #
        # me_group.add_argument("--dDNAs", nargs="+", metavar="*.fasta", type=str,
        #     help="Homology arms, presence of unique gRNA target site, and which \
        #          features these dDNAs would disrupt, and whether-or-not the \
        #          dDNA would introduce a mutation.")
        #
        # me_group.add_argument("--primers", metavar="*.fasta",
        #     help="Evaluate what amplicon sizes to expect for wt/ko/ki for \
        #          input primers. Also evaluate the Tm/sensitivity/dimerization/etc \
        #          for the input primers. Does not evaluate these in multiplex.")
        
        # Maybe: make --gRNAs, --dDNAs, and --primers all accept the same number of FASTA files
        #   Then, when it does calculations, it does all of them in serial steps.
        #  step 1: --gRNAs[0], --dDNAs[0], --primers[0]
        #  step 2: --gRNAs[1], --dDNAs[1], --primers[1]
        #  etc...
        
        # Add optional arguments

        self.parser.add_argument("--homologs", metavar="*.homologs", type=str, default=None,
            help="Path to text file containing homologous features on the same \
            line, separated by TAB characters")

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

        self.parser.add_argument("--tag", metavar='TAG', type=str, default='ID',
            help="GFF tag with feature names. Examples: 'ID', 'Name', 'Gene', 'Parent', or 'locus_tag'")

        self.parser.add_argument("--selection", metavar='FEATURE', nargs="+", type=str, default=None, #default=argparse.SUPPRESS, # '==SUPPRESS=='
            help="Select only certain features rather than all features in input GFF file.")

        self.parser.add_argument("--ambiguities", type=str,
            choices=["exclusive", "discard", "disambiguate", "keep"],
            default="discard",
            help="How generated gRNAs should treat ambiguous bases: \
            exclusive - gRNAs will only be created for ambiguous locations; \
            discard - no gRNAs will be created where the FASTA has an ambiguous base; \
            disambiguate - gRNAs containing ambiguous bases will be converted to a set of non-ambiguous gRNAs; \
            keep - gRNAs can have ambiguous bases.")

        self.parser.add_argument("--case", type=str, default="ignore",
            choices=["ignore", "upper-only", "lower-only", "mixed-lower", "mixed-upper", "mixed-only"],
            help="Restrict generation of gRNAs based on case of nucleotides in input FASTA: \
            ignore - keep all spacers; \
            upper-only - discard any potential spacer with a lower-case character in the genome; \
            lower-only - discard all potential spacers with upper-case characters; \
            mixed-lower - discard spacers that are all upper-case; \
            mixed-upper - discard spacers that are all lower-case; \
            mixed-only - only use spacers that have both lower- and upper-case characters.")

        self.parser.add_argument("--features", metavar="FEATURE", type=str, nargs="+", default=["gene"],
            help="Features to design gRNAs against. Must exist in GFF file. \
            Examples: 'CDS', 'gene', 'mRNA', 'exon', 'intron', 'tRNA', 'rRNA'.\
            The special 'all' feature type will include all listed features.")

        self.parser.add_argument("--dDNA_gDNA_ratio", metavar="N", type=int, default=1000,
            help="Ratio of donor DNA to genomic DNA for calculating off-target scores")

        self.parser.add_argument("--processors", metavar="N", type=int, default=(os.cpu_count() or 1),
            help="Number of processors to use when performing pairwise sequence alignments.")

        aligner_choices = [x.name for x in aligners.aligners]
        self.parser.add_argument("--aligner", type=str, choices=aligner_choices, default='bowtie2',
            help="Program to calculate pairwise alignments. Please note that the 'addtag' internal aligner is very slow.")

        prefilter_choices = [C.name for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.prefilter]
        self.parser.add_argument("--prefilters", nargs='+', type=str,
            choices=prefilter_choices, default=['GC', 'PolyT'],
            help="Specific algorithms for determining gRNA goodness.")

        off_target_choices = [C.name for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.off_target]
        self.parser.add_argument("--offtargetfilters", nargs='+', type=str,
            choices=off_target_choices, default=['CFD', 'Hsu-Zhang'],
            help="Specific algorithms for determining gRNA goodness.")

        on_target_choices = [C.name for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.on_target]
        self.parser.add_argument("--ontargetfilters", nargs='+', type=str,
            choices=on_target_choices, default=['Azimuth'],
            help="Specific algorithms for determining gRNA goodness.")

        postfilter_choices = [C.name for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.postfilter]
        self.parser.add_argument("--postfilters", nargs='+', type=str,
            choices=postfilter_choices, default=['Errors'],
            help="Specific algorithms for determining gRNA goodness.")

        self.parser.add_argument("--target_specificity", choices=['exclusive', 'all', 'any'], default='all',
            help="For 'exclusive', only spacer sequences that diagnostically target alleles will be designed. \
            With this option enabled, spacers will target polymorphisms in the feature. \
            For 'all', spacers will target invariant sites within the feature.\
            For 'any', targets will not be checked against homologous features.")

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
        

    def compute(self, args):
        '''UNDER DEVELOPMENT'''
        print("Perform the evaluation here.")

        # Load the FASTA file specified on the command line
        contig_sequences = utils.load_multiple_fasta_files(args.fasta)

        # Open and parse the GFF file specified on the command line
        #features = utils.load_gff_file(args.gff, args.features, args.tag)
        # Filter features by what is selected
        #features = self.filter_features(features, args.selection)
        for gff_file in args.gff:
            Feature.load_gff_file(gff_file, args.features, args.excluded_features, args.selection, args.tag)
        Feature.assert_features(args.selection, contig_sequences)

        # Load '--homologs' file
        # Make dict linking each feature to its gene
        # Make dict linking features to each other as homologs
        if args.homologs:
            homologs, feature2gene = utils.load_homologs(args.homologs)
        else:
            homologs, feature2gene = utils.dummy_homologs()

        # Assign default 'self.homologs' and 'self.gene'
        Feature.assign_homologs(homologs)
        Feature.assign_gene(feature2gene)

        self.logger.info('Feature.features')
        for f_name, f in sorted(Feature.features.items()):
            self.logger.info("  {}:{}:{}:{}..{}".format(f.name, f.contig, f.strand, f.start, f.end))
        self.logger.info('Feature.excluded_features')
        for exf_name, f in sorted(Feature.excluded_features.items()):
            self.logger.info("  {}:{}:{}:{}..{}".format(f.name, f.contig, f.strand, f.start, f.end))

        spacers = ExcisionTarget.load_spacers(args, args.spacers) # Maybe should be 'ExcisionTarget.load_fasta()'?

        # Re-format input spacers as their own FASTA
        ex_query_file = ExcisionTarget.generate_query_fasta(os.path.join(args.folder, 'excision-query.fasta'))

        # Re-format input genomes
        genome_fasta_file = utils.write_merged_fasta(contig_sequences, os.path.join(args.folder, 'genome.fasta'))

        # Index the genome
        genome_index_file = args.selected_aligner.index(genome_fasta_file, os.path.basename(genome_fasta_file), args.folder, args.processors)

        # Use selected alignment program to find all matches in the genome and dDNAs
        q2gDNA_align_file = args.selected_aligner.align(ex_query_file, genome_index_file, 'excision-query-2-gDNA.'+args.selected_aligner.output, args.folder, args.processors)

        # Load the SAM files and add Alignments to ExcisionTarget sequences
        ExcisionTarget.load_alignment(q2gDNA_align_file, args, contig_sequences)

        # Calculate off-target/guide scores for each algorithm
        self.logger.info("ExcisionTarget after SAM parsing and off-target scoring")
        ##### Add short-circuit/heuristic #####
        for et_seq, et_obj in ExcisionTarget.sequences.items():
            et_obj.score_off_targets(args, homologs)
            self.logger.info(et_obj)
            for a in et_obj.alignments:
                self.logger.info('  ' + str(a))

        # Batch calculate with new ExcisionTarget class
        # TODO: Ideally, all scoring Algorithms will be run in batch
        ExcisionTarget.score_batch(args)

        self.logger.info("ExcisionTarget after Azimuth calculation")
        for et_seq, et_obj in ExcisionTarget.sequences.items():
            self.logger.info(str(et_obj) + '\t' + '\t'.join(map(str, et_obj.locations)))

        # Generate the FASTA with the final scores
        excision_spacers_file = ExcisionTarget.generate_spacers_fasta(os.path.join(args.folder, 'excision-spacers.fasta'))

        # Write scores to STDOUT
        self.write_output()
        
        # End 'compute()'

    def write_output(self):
        pass
