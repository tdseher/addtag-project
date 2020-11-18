#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/subroutines/_subroutine_generate_all.py

# TODO: Add a new output that shows a graphical view of the candidate exTargets against the input Features
#       and make it output by default to STDOUT! (not just to log)
#       For example:
#         Feature_01
#         position            10        20        30        40        50        60        70        80        90        100
#                      GAGTAGCAGTTATAAATACTAATCTTTTTGTTGTTAGATATATGAAGCCATCTTGTGCTTGAATCGACACTTTTCTCTCATTAGGATTAGGATTTGTGGT
#         exTarget-1                                               3'-CCATCTTGTGCTTGAATCGA-5'
#         exTarget-2                                                               5'-TCGACACTTTTCTCTCATTAGG-3'

# TODO: Add Tm calculation as an Algorithm
# TODO: Add full and partial (N{20}>NGG, N{12}>NGG, N{8}>NGG) alignments for finding off-targets
#       Thus there will be 3 off-target scores: the full score, the 12nt score, and the 8nt score.
#       This approach is derived from 'CRISPRdirect' software
# TODO: Report NUMBER of off-target alignments in the STDOUT table, and not just the Algorithm scores
# TODO: Make it so that an Algorithm can return more than a single score
#       (for instance, a 'Composition' algoritm can return both %GC and %AT)
#       Each return value would be added as columns in STDOUT table

# TODO: Make it so that a star-tag can have any combination of these elements in its insert:
#         [x] addtag (RGN Target)
#         [x] bartag (guaranteed edit distance from others)
#         [x] sigtag (restriction site)
#         [x] flanktags 

# Import standard packages
import sys
import os
import logging

# Import included AddTag-specific modules
#from .__init__ import Feature, Donor, ExcisionDonor, ReversionDonor, Target, ExcisionTarget, ReversionTarget
#from . import __init__
#from source import Feature, Donor, ExcisionDonor, ReversionDonor, Target, ExcisionTarget, ReversionTarget
#from .__init__ import Donor, ExcisionDonor, ReversionDonor, Target, ExcisionTarget, ReversionTarget
from . import subroutine
from ..feature import Feature
from ..targets import ExcisionTarget, ReversionTarget
from ..donors import ExcisionDonor, ReversionDonor
from .. import utils
from .. import algorithms
from .. import aligners
from .. import thermodynamics

logger = logging.getLogger(__name__)

class GenerateAllParser(subroutine.Subroutine):
    logger = logger.getChild('GenerateParser')
    
    def __init__(self, subparsers):
        self.subparsers = subparsers
        
        self.name = 'generate_all'
        self.description = (
            "description:" "\n"
            "  Design full sets of oligonucleotide sequences for CRISPR/Cas genome" "\n"
            "  engineering experiment." "\n"
            "  " "\n"
            "  The objective of the 'generate' subroutine is to identify gDNA" "\n"
            "  Targets within or near to specific genomic Features that can be" "\n"
            "  cut with an experimenter’s RGN of choice, and then produce" "\n"
            "  dDNAs that will replace the Features." "\n"
        )
        self.help = "Design full sets of oligonucleotide sequences for CRISPR/Cas genome engineering experiment."
        self.epilog = (
            "example:" "\n"
            "  $ python3 {__program__} {__subroutine__} --fasta chromosomes.fasta --gff GENE.gff --folder GENEg" "\n"
            "    --homologs homologs.txt --ko-dDNA mintag --ko-gRNA --ki-gRNA --ki-dDNA --postfilters " "\n"
            "    --offtargetfilters CFD > GENEg.out 2> GENEg.err" "\n"
        ).format(**dict(list(subroutine.__dict__.items()) + list({"__subroutine__": self.name}.items()))) # key:value pairs in the latter dict will replace any instances in the prior dict
        
        self.define_parser()
        self.define_arguments()
    
    def define_arguments(self):
        # Add required arguments
        required_group = self.parser.add_argument_group('required arguments')
        required_group.add_argument("--fasta", required=True, nargs="+", metavar="*.fasta", type=str,
            help="FASTA files with contigs to find unique RGN target sites. Input FASTA \
            must not be compressed. User can decide whether ambiguous bases can \
            be chosen for RGN target sites. All FASTA sequences should have unique \
            primary headers (everything between the '>' symbol and the first \
            whitespace should be unique). You should include FASTA of the genome \
            and any plasmids.")
        required_group.add_argument("--gff", required=True, nargs="+", metavar="*.gff", type=str,
            help="GFF files specifying chromosomal features that should be \
            targeted (in multiplex) or excluded.")
        required_group.add_argument("--folder", required=True, metavar="FOLDER",
            type=str, help="Path of folder to store generated files.")
        
        # Add optional arguments
        self.parser.add_argument("--exhaustive", action="store_true", 
            help="Perform brute force search for optimal RGN target design. \
            This will significantly increase runtime.") # Default to on when not trimming???
        #parser.add_argument("--feature_homolog_regex", metavar="REGEX", type=str, default=None, help="regular expression with capturing group containing invariant feature. Example: '(.*)_[AB]' will treat features C2_10010C_A and C2_10010C_B as homologs")
        # okay idea, but needs more thought before implementation
        self.parser.add_argument("--homologs", metavar="*.homologs", type=str, default=None,
            help="Path to text file containing homologous features on the same \
            line, separated by TAB characters")
        #parser.add_argument("--pams", metavar="SEQ", nargs="+", type=str,
        #    default=["NGG"], help="Constrain finding only targets with these PAM sites")
        #parser.add_argument("--target_lengths", nargs=2, metavar=('MIN', 'MAX'),
        #    type=int, default=[17, 20],
        #    help="The length range of the 'target'/'spacer'/gRNA site")
        # Replacement for --pams and --target_lengths with this:
        self.parser.add_argument("--motifs", metavar="MOTIF", nargs="+", type=str,
            default=["N{17}|N{3}>NGG"],
            help="Find only targets with these 'SPACER>PAM' motifs, written from \
            5' to 3'. '>' points toward PAM. Ambiguous nucleotide characters are accepted. '{a,b}' \
            are quantifiers. '(a,b,…)' are permitted alternatives. \
            '/' is a sense strand cut, '\\' is an antisense strand \
            cut, and '|' is a double-strand cut. '.' is a base used for positional \
            information, but not enzymatic recognition. Be sure to enclose each \
            motif in quotes so your shell does not interpret STDIN/STDOUT redirection.")
        self.parser.add_argument("--off_target_motifs", metavar="MOTIF", nargs="+", type=str,
            default=[],
            help="Defaults to the same as the on-target motif. Definition syntax is identical.")
        self.parser.add_argument("--protruded_targets", action="store_true", default=False,
            help="Allow RGN targets to only partially overlap input Features.")
        # TODO: Add '--constructs FILE' option for automatically stitching together SPACER & SCAFFOLD
        # Need to decide if construct inputs should be TSV, or FASTA
        # And whether or not there should be an upstream parameter separate from
        # a downstream one. or if they are the same, then what?
        #self.parser.add_argument("--constructs", metavar="*.fasta", nargs="+", type=str,
        #    default=[], help="The first sequence will be prepended, and the second \
        #    sequence will be appended to the generated spacer sequences to form \
        #    the construct sequences. It is useful to put the gRNA promotor as the \
        #    first sequence, and the scaffold sequence and terminator as the \
        #    second. Specify one FASTA file for each motif.")
        self.parser.add_argument("--tag", metavar='TAG', type=str, default='ID',
            help="GFF tag with feature names. Examples: 'ID', 'Name', 'Gene', 'Parent', or 'locus_tag'")
        self.parser.add_argument("--selection", metavar='FEATURE', nargs="+", type=str, default=None, #default=argparse.SUPPRESS, # '==SUPPRESS=='
            help="Select only certain features rather than all features in input GFF file.")
        self.parser.add_argument("--ambiguities", type=str,
            choices=["exclusive", "discard", "disambiguate", "keep"],
            default="discard",
            help="How generated RGN targets should utilize ambiguous bases: \
            exclusive - targets will only be identified at ambiguous locations; \
            discard - targets cannot overlap ambiguous bases; \
            disambiguate - targets containing ambiguous bases will be converted to a set of non-ambiguous sequences; \
            keep - targets can have ambiguous bases.")
        self.parser.add_argument("--case", type=str, default="ignore",
            choices=["ignore", "upper-only", "lower-only", "mixed-lower", "mixed-upper", "mixed-only"],
            help="Restrict identification of RGN targets based on case of nucleotides in input FASTA: \
            ignore - keep all targets; \
            upper-only - discard any potential target with a lower-case character in the genome; \
            lower-only - discard all potential target with upper-case characters; \
            mixed-lower - discard targets that are all upper-case; \
            mixed-upper - discard targets that are all lower-case; \
            mixed-only - only use targets that have both lower- and upper-case characters.")
        #parser.add_argument("--min_contig_edge_distance", metavar="N", type=int, default=500,
        #    help="Minimum distance from contig edge a site can be found")
        self.parser.add_argument("--features", metavar="FEATURE", type=str, nargs="+", default=["gene"],
            help="Features to design RGN targets against. Must exist in GFF file. \
            Examples: 'CDS', 'gene', 'mRNA', 'exon', 'intron', 'tRNA', 'rRNA'.\
            The special 'all' feature type will include all listed features.")
        #self.parser.add_argument("--warning_features", metavar='FEATURE', nargs="+", type=str, default=['all'],
        #    help="GFF tags that will trigger a warning if they overlap with the \
        #    target feature. Examples: 'CDS', 'gene', 'mRNA', 'exon', 'intron', 'tRNA', 'rRNA'")
        self.parser.add_argument("--excluded_features", metavar='FEATURE', nargs='+', type=str, default=[],
            help="Prevent feature expansion if it would overlap one of these features present in the GFF. \
            Examples: 'CDS', 'gene', 'mRNA', 'exon', 'intron', 'tRNA', 'rRNA'.\
            The special 'all' feature type will exclude all listed features \
            except those specified by the '--features' option.")
        self.parser.add_argument("--dDNA_gDNA_ratio", metavar="N", type=int, default=1000,
            help="Ratio of donor DNA to genomic DNA for calculating off-target scores")
        #self.parser.add_argument("--target_gc", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[25, 75],
        #    help="Generated gRNA spacers must have %%GC content between these values (excludes PAM motif)") # Moved to prefilter
        
        #  center_feature: --------HHHH[...............FEATURE.........TARGET]HHHH------------------------
        #                  -----HHHH[..................FEATURE.........TARGET...]HHHH--------------------- pad=3
        #   center_target: -----------------------HHHH[FEATURE.........TARGET................]HHHH--------
        #                  --------------------HHHH[...FEATURE.........TARGET...................]HHHH----- pad=3
        #     center_both: -----------------------HHHH[FEATURE.........TARGET]HHHH------------------------
        #                  --------------------HHHH[...FEATURE.........TARGET...]HHHH--------------------- pad=3
        # justify_feature: -----------------------HHHH[FEATURE.........TARGET]HHHH------------------------
        #                  -----------------------HHHH[FEATURE.........TARGET...]HHHH--------------------- pad=3
        #  justify_target: -----------------------HHHH[FEATURE.........TARGET]HHHH------------------------
        #                  --------------------HHHH[...FEATURE.........TARGET]HHHH------------------------ pad=3
        self.parser.add_argument("--feature_expansion_format", type=str, default=None,
            choices=['center_feature', 'center_target', 'center_both', 'justify_feature', 'justify_target'],
            help="If a feature needs to be expanded to contain an RGN target, \
            expand the feature such that either the feature, the target, or \
            both the feature and the target are in the center.")
        self.parser.add_argument("--feature_expansion_lengths", nargs=2, metavar=("MIN", "MAX"), type=int, default=[100, 4000],
            help="If a feature needs to be expanded to contain an RGN target (SPACER>PAM site), \
            the length of the expanded feature must be within this range, inclusive. \
            If a derived feature is smaller than MIN, then it is padded appropriately \
            until it reaches the minimum size.")
        self.parser.add_argument("--feature_expansion_pad", metavar="N", type=int, default=None,
            help="If a feature needs to be expanded to contain an RGN target, \
            and the expanded feature is within the permitted lengths, then \
            the expanded feature will have N number of nucleotides padded.")
        #    This is useful when the expanded feature would be small, and should \
        #    be artificially enlarged to ensure decent iF/iR primer amplification. \
        #    If left unspecified, feature expansion will be automatically applied \
        #    to small features (<300 nt).")
        
        # TODO: Consider removing '--excise_upstream_homology' and '--excise_downstream_homology' and replacing them
        #       entirely with '--excise_donor_lengths'
        self.parser.add_argument("--excise_upstream_homology", nargs=2, metavar=("MIN", "MAX"), type=int, default=[50,50],
            help="Range of homology lengths acceptable for knock-out dDNAs, inclusive.")
        self.parser.add_argument("--excise_downstream_homology", nargs=2, metavar=("MIN", "MAX"), type=int, default=[47,50],
            help="Range of homology lengths acceptable for knock-out dDNAs, inclusive.")
        self.parser.add_argument("--homology_distance", nargs=3, metavar=("M", "G", "E"), type=int, default=[0, 0, 0], # formerly '--max_homology_errors'
            help="Maximum number of mismatches (M), gaps (G), and total errors (E) to allow in \
            each dDNA homology arm across all homologous Features. \
            If a number of M/G/E greater than the maximum exists, then the feature will be \
            expanded until a homology region with appropriate M/G/E is found.")
        self.parser.add_argument("--target_distance", nargs=3, metavar=("M", "G", "E"), type=int, default=[0, 0, 0],
            help="Maximum number of mismatches (M), gaps (G), and total errors (E) to consider \
            Targets in homologous Features to be equivalent.")
        
        # TODO: Change '--excise_donor_lengths MIN MAX' so it accepts only a single length, and becomes
        #       '--excise_donor_length N'
        self.parser.add_argument("--excise_donor_lengths", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[100, 100],
            help="Range of lengths acceptable for knock-out dDNAs, inclusive.")
        self.parser.add_argument("--excise_insert_lengths", nargs=2, metavar=("MIN", "MAX"), type=int, default=[0, 3],
            help="Range for inserted DNA lengths, inclusive (mintag). \
            If MIN < 0, then regions of dDNA homology (outside the feature) will be removed.")
        self.parser.add_argument("--excise_feature_edge_distance", metavar="N", type=int, default=0,
            help="If positive, gRNAs won't target any nucleotides within this distance \
            from the edge of the feature. If negative, gRNAs will target nucleotides \
            this distance outside the feature.")
        #
        # May need to remove '--excise_upstream_feature_trim' and '--excise_downstream_feature_trim'
        # And replace with '--expand_feature {upstream,downstream,both}'
        self.parser.add_argument("--excise_upstream_feature_trim", nargs=2, metavar=('MIN', 'MAX'),
            type=int, default=[0, 0], help="Between MIN and MAX number of nucleotides \
            upstream of the feature will be considered for knock-out when designing \
            donor DNA.")
        self.parser.add_argument("--excise_downstream_feature_trim", nargs=2, metavar=("MIN", "MAX"),
            type=int, default=[0, 0], help="Between MIN and MAX number of nucleotides \
            downstream of the feature will be considered for knock-out when designing \
            donor DNA.")
        #self.parser.add_argument("--expand_feature", type=str, default="both",
        #    choices=["none", "upstream", "downstream", "both"],
        #    help="Sequence immediately up/down-stream of intended feature may \
        #    be knocked-out in order to generate good mintag sites.")
        #
        #
        #parser.add_argument("--min_donor_insertions", metavar="N", type=int, default=2,
        #    help="The uniqueness of final donor DNA compared to the rest of the genome")
        #parser.add_argument("--min_donor_deletions", metavar="N", type=int, default=2,
        #    help="The uniqueness of final donor DNA compared to the rest of the genome")
        #parser.add_argument("--min_donor_substitutions", metavar="N", type=int, default=2,
        #    help="The uniqueness of final donor DNA compared to the rest of the genome")
        #parser.add_argument("--min_donor_errors", metavar="N", type=int, default=3,
        #    help="The uniqueness of final donor DNA compared to the rest of the genome")
        
        #self.parser.add_argument("--revert_upstream_homology", nargs=2, metavar=("MIN", "MAX"), type=int, default=[100,400],
        #    help="Range of homology lengths acceptable for knock-in dDNAs, inclusive.")
        #self.parser.add_argument("--revert_downstream_homology", nargs=2, metavar=("MIN", "MAX"), type=int, default=[100,400],
        #    help="Range of homology lengths acceptable for knock-in dDNAs, inclusive.")
        self.parser.add_argument("--revert_homology_length", nargs=2, metavar=("MIN", "MAX"), type=int, default=[100,400],
            help="Range of homology lengths acceptable for knock-in dDNAs, inclusive.")
        self.parser.add_argument("--revert_amplification_primers", action="store_true", default=False,
            help="Do PCR calculations for amplifying the wild-type allele.")
        
        #self.parser.add_argument("--revert_donor_lengths", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[0, 100000],
        #    help="Range of lengths acceptable for knock-in dDNAs.")
    #    parser.add_argument("--min_donor_distance", metavar="N", type=int, default=36,
    #        help="The minimum distance in bp a difference can exist from the edge of donor DNA") # homology with genome
        #self.parser.add_argument("--max_consecutive_ts", metavar="N", type=int, default=4, # Moved to prefilter algorithm
        #    help="The maximum number of Ts allowed in generated gRNA sequences.")
        self.parser.add_argument("--max_number_sequences_reported", metavar="N", type=int, default=5,
            help="The maximum number of sequences to report for each step.")
        self.parser.add_argument("--min_weight_reported", metavar="N", type=float, default=0.01,
            help="Only gRNA-dDNA pairs with at least this much weight will be reported.")
        # program currently will only search 'both' strands
        #parser.add_argument("--strands", type=str, choices=["+", "-", "both"], default="both",
        #    help="Strands to search for gRNAs")
        
        # Add command line arguments for the additional hard constraints:
        #  Only report potential targets that have no off targets with mismatches within 8, 12, N nt from 3' end
        self.parser.add_argument("--processors", metavar="N", type=int, default=(os.cpu_count() or 1),
            help="Number of processors to use when performing pairwise sequence alignments.")
        
        aligner_choices = [x.name for x in aligners.pw_aligners]
        self.parser.add_argument("--aligner", type=str, choices=aligner_choices, default='bowtie2',
            help="Program to calculate pairwise alignments. Please note that the 'addtag' internal aligner is very slow.")
        
        oligo_choices = [x.name for x in thermodynamics.oligos]
        self.parser.add_argument("--oligo", type=str, choices=oligo_choices, default='Primer3',
            help="Program to perform thermodynamic calculations.")
        
        prefilter_choices = [C.name for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.prefilter]
        self.parser.add_argument("--prefilters", nargs='+', type=str,
            choices=prefilter_choices, default=['GC', 'PolyT'],
            help="Specific algorithms for determining RGN target goodness.")
        
        off_target_choices = [C.name for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.off_target]
        self.parser.add_argument("--offtargetfilters", nargs='+', type=str,
            choices=off_target_choices, default=['CFD', 'Hsu-Zhang'],
            help="Specific algorithms for determining RGN target goodness.")
        
        on_target_choices = [C.name for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.on_target]
        self.parser.add_argument("--ontargetfilters", nargs='+', type=str,
            choices=on_target_choices, default=['Azimuth'],
            help="Specific algorithms for determining RGN target goodness.")
        
        postfilter_choices = [C.name for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.postfilter]
        self.parser.add_argument("--postfilters", nargs='*', type=str,
            choices=postfilter_choices, default=['Errors'],
            help="Specific algorithms for determining RGN target goodness.")
        
        self.parser.add_argument("--weights", nargs='+', type=str, default=[], metavar='ALGORITHM:INFLECTION±SLOPE,*',
            help="Specify any number of non-default Algorithm weights with the following syntax. \
            If SLOPE>0, then weight increases with score. If SLOPE<0, then weight decreases as score increases. \
            The ALGORITHM name followed by a colon (':'). Then a comma-separated list of INFLECTION point, \
            sign ('+' or '-'), and SLOPE for the sigmoidal curve. For example, '--weights GC:40+1.7,60-1.7' \
            will resemble a bell-shaped curve with peak weight at scores of 50.")
        
        self.parser.add_argument("--primer_scan_limit", metavar="N", type=int, default=2*60,
            help="Amount of time (in seconds) to limit each primer scan.")
        
        self.parser.add_argument("--primer_pair_limit", metavar="N", type=int, default=5*60,
            help="Amount of time (in seconds) to limit primer pairings.")
        
        # Temporary stopgap to make sure the calculations don't take too long
        #self.parser.add_argument("--max_time", metavar="SECONDS", type=float, default=60,
        #    help="Maximum amount of time, in seconds, for each feature to spend calculating dDNAs.")
        
        self.parser.add_argument("--bartag_number", metavar="N", type=int, default=1,
            help="Number of bartags per locus to generate. \
            For efficiency's sake, a maximum 'features × bartags = 200' is enforced.")
        
        self.parser.add_argument("--bartag_motif", metavar="MOTIF", type=str, default='N{10}',
            help="Structure of nucleotides that should be generated. \
            Longer oligonucleotide lengths take longer to calculate.")
        
        self.parser.add_argument("--bartag_distance", metavar="N", type=int, default=3,
            help="Minimum required edit distance between each bartag.")
        
        self.parser.add_argument("--flanktags", nargs='+', action=subroutine.ValidateFlanktags, type=str,
            metavar=('{*.fasta,uniform,specific}', '{single,pair}'), default=None,
            help="If 'uniform', then all features will share the same flanktags. \
            If 'specific', then each feature will have its own flanktags. \
            If '*.fasta', then flanktags will be taken from an input FASTA file with either a single or paired flanktags. \
            If 'single', then design a single oligo to serve as both the uptag and dntag. \
            If 'pair', then design uptag different to dntag.")
        
        # Splitting up "generate" into specific sub-tasks
        
        # TODO: Running 'addtag generate_all --ko-dDNA mintag' without '--ko-gRNA' prints nothing to STDOUT.
        #       This needs to be fixed.
        self.parser.add_argument("--ko-gRNA", action='store_true', default=False,
            help="Design gRNAs to target features in genome")
        
        #self.parser.add_argument("--ko-dDNA", type=str, default=None,
        #    choices=['mintag', 'addtag', 'unitag', 'bartag'],
        #    help="'mintag' are unique us/i/ds junction targets specific to each feature. \
        #    'addtag' are unique targets for each feature. \
        #    'unitag' is a single, invariant target for ALL features. \
        #    'bartag' are unique barcodes for each feature (does not guarantee targets).")
        
        # TODO: If '--ko-dDNA X' is not specified on the commandline, then 'args.ko_dDNA' takes the value 'None'
        #       If it is 'None', then errors will happen. Need to fix this
        self.parser.add_argument("--ko-dDNA", type=str, action=subroutine.ValidateKodDNA,
            metavar='{*.fasta,mintag,addtag,unitag,bartag,sigtag}', default=None,
            help="'*.fasta' is a FASTA file containing user-specified sequences. \
            'mintag' are unique us/i/ds junction targets specific to each feature. \
            'addtag' are unique targets for each feature. \
            'unitag' is a single, invariant target for ALL features. \
            'bartag' are unique barcodes for each feature (does not guarantee targets). \
            'sigtag' are sequences containing restriction sites (does not guarantee targets).")
        
        self.parser.add_argument("--ki-gRNA", action='store_true', default=False,
            help="Design gRNAs to target the ko-dDNA. \
            Defaults to True if '--ko-dDNA mintag' is specified.")
        
        self.parser.add_argument("--ki-dDNA", nargs='?', type=str, default=None,
            metavar='*.fasta', const=True, # If no command-line argument follows, the value of const will be assumed instead.
            help="If no file is specified, then design wild type dDNA. \
            If a FASTA file is specified, then design dDNA to replace features \
            with the sequences in this FASTA file. The primary sequence header \
            should either be the TAG identifier for the feature to replace \
            (specified by the '--tag' option) or the gene name (specified \
            within the '--homologs' file).")
        #    help="If no file is specified, then design wild type dDNA. \
        #    If a FASTA file is specified, then design dDNA to replace features \
        #    with the sequences in this FASTA file. The sequence header should \
        #    have a 'TAG=FEATURE' field where 'TAG' is the tag specified with \
        #    '--tag' option, and 'FEATURE' corresponds to the value of that \
        #    'TAG' within the input GFF.")
        
        # subparsers
        #  addtag generate
        #   --ko-gRNA                              Design gRNAs to target features in genome
        #   --ko-dDNA mintag|addtag|unitag|bartag  Design knock-out dDNAs for target features in genome
        #                                          (these will be constrained to dDNAs targetable by ki-gRNAs)
        #             mintag                          Designs unique us/i/ds junction target specific to each feature
        #                    addtag                   Designs unique target for each feature
        #                           unitag            Designs a single/uniform, unique target for ALL features
        #                                  bartag     Designs unique barcode for each feature
        #   --bartag_length INT                    Specify length of sigtag barcode in nt (Or should this be --sigtag_motif NNNNNN?)
        #   --flanktags generate                   Design identical uptag/dntag primer sequences to addtag/sigtag/unitag for each site/feature
        #   --flanktag_length INT INT              Specify length of generated uptag and dntag sequences in nt
        #   --flanktags input FASTA                Use file specifying uptag/dntag primers sequences for addtag/sigtag for all sites
        #                                          >seq1 ID=feature1 flanktag=uptag
        #                                          NNNNNNNNNNNNNNNNNNNN
        #                                          >seq2 ID=feature1 dntag
        #                                          NNNNNNNNNNNNNNNNNNNN
        #   
        #   --ki-gRNA     Design knock-in gRNAs that target the ko-dDNA
        #   --ki-dDNA wt        Design the dDNA for knocking the wild type feature back in
        #   --ki-dDNA FASTA     If a FASTA file is specified, then use mutants for these indicated features
        #                       >seq1 ID=feature1
        #                       NNNNNNNNNNNNNNNNNNN
        #                       The entire feature will be replaced by this sequence
        
        self.parser.add_argument("--target_specificity", choices=['exclusive', 'all', 'any'], default='all',
            help="For 'exclusive', only spacer sequences that diagnostically target alleles will be designed. \
            With this option enabled, spacers will target polymorphisms in the feature. \
            For 'all', spacers will target invariant sites within the feature.\
            For 'any', targets will not be checked against homologous features.")
        
        # ^^^^vvvv These two options, ideally, should be only one option ^^^^vvvv
        #self.parser.add_argument("--allele-specific-primers", action="store_true", default=False,
        #    help="If primers are calculated for amplifying the wild type, then they should be allele-specific.")
        
        self.parser.add_argument("--donor_specificity", choices=['exclusive', 'all', 'any'], default='all',
            help="For 'exclusive' (allele-specific), homology arms of dDNAs should be unique for each homologous feature. \
            Primer amplicons will either be diagnostially-different sizes, or primer sequences themselves will be different. \
            For 'all' (multi-allelic), dDNA homology arms will minimize polymorphisms so they target all homologous features, and primer pairs will amplify all homologs. \
            For 'any' (allele-agnostic), dDNA homology arms and primer pairs will be calculated for each feature independently.")
        
        # Future implementation: allow for user to input a choice for each dDNA
        #   --donor_specificity all exclusive      --donor_specificity exclusive all
        #     A0   B0   (g0: wild type)              A0   B0   (g0: wild type)
        #      \   /      1 ko-dDNA                  |     |     2 ko-dDNAs
        #       AB1     (g1: intermediary)           A1   B1   (g1: intermediary)
        #      /   \      2 ki-dDNAs                  \   /      1 ki-dDNA
        #     A2   B2   (g2: knock-in)                 AB2     (g2: knock-in)
        #
        # But the current implementation will do only the following:
        #   --donor_specificity all all            --donor_specificity exclusive exclusive
        #     A0   B0   (g0: wild type)              A0   B0   (g0: wild type)
        #      \   /      1 ko-dDNA                  |     |     2 ko-dDNAs
        #       AB1     (g1: intermediary)           A1   B1   (g1: intermediary)
        #        |        1 ki-dDNA                  |     |     2 ki-dDNAs
        #       AB2     (g2: knock-in)               A2   B2   (g2: knock-in)
        
        # Get TEMP directory
        if sys.platform.startswith('win'):
            from tempfile import gettempdir
            temp_default = gettempdir()
        else: # 'linux' or 'darwin'
            temp_default = '/dev/shm'
        self.parser.add_argument("--temp_folder", type=str, default=temp_default,
            help="Directory to store temporary files. RAMdisk recommended.")
        
    def compute(self, args):
        """Perform complete CRISPR/Cas analysis for input"""
        #contigs = utils.load_fasta_file(args.fasta[0]) # To do: Do this for all FASTA files, then merge them into the same dictionary
        
        # Load the FASTA file specified on the command line
        # Merge all sequence information into the same dictionary
        #fasta_index, contig_index, contig_sequences = utils.load_indexed_fasta_files(args.fasta)
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
            self.logger.info("    homologs={}".format(f.homologs))
        self.logger.info('Feature.excluded_features')
        for exf_name, f in sorted(Feature.excluded_features.items()):
            self.logger.info("  {}:{}:{}:{}..{}".format(f.name, f.contig, f.strand, f.start, f.end))
        
        # Run Features through a check. If homologous features are on the same contig, then stop
        # and print an error/warning that they won't be handled. And the user should split them
        # into separate contigs.
        h_groups = Feature.group_features_by_gene(feature2gene)
        #for k, v in h_groups.items():
        #    print(k, v, file=sys.stderr)
        # EFG1 [Feature(name=CR_07890W_A, gene=EFG1, location=Ca22chrRA_C_albicans_SC5314:+:1723488..1725146), Feature(name=CR_07890W_B, gene=EFG1, location=Ca22chrRB_C_albicans_SC5314:+:1722891..1724556)]
        h_exit = False
        for h_gene, h_features in h_groups.items():
            if (len(set([f.contig for f in h_features])) != len(h_features)):
                print('Homologous Features for gene {} lie on the same contig: {}'.format(h_gene, h_features), file=sys.stderr)
                h_exit = True
        if h_exit:
            raise Exception('Some homologs are on the same contig, and thus cannot be processed. \
            All homologous Features must be on separate contigs.')
        
        # Merge features?
        #features = merge_features(features)
        
        #### Some checks that need to be added ####
        # The feature being targeted MUST not go up to the edge of the contigs
        # Otherwise there will be no junction. i.e. the upstream and downstream
        # regions won't exist.
        # Actually, these regions must be at minimum, 50 nt
        ###########################################
        
        # This code is performed in 'ExcisionDonor.generate_donors()'
        # if (args.ko_dDNA): # (mintag/addtag/unitag/bartag)
        #     if (args.ko_dDNA == 'mintag'):
        #         pass
        #     elif (args.ko_dDNA == 'addtag'):
        #         pass
        #     elif (args.ko_dDNA == 'unitag'):
        #         pass
        #     elif (args.ko_dDNA == 'bartag'):
        #         pass
        
        # Search for good targets within specified features
        # If no good target is found, then expand the feature
        # TODO: Need to de-couple the Target search and the homology search.
        #       if '--ko-dDNA', then homology search should happen
        #       if '--ko-gRNA', then Target search should happen
        if (args.ko_dDNA or args.ko_gRNA):
            # Only expand features if '--feature_expansion_format' is not 'None'
            if (args.feature_expansion_format != None):
                #Feature.expand_all_features(args, contig_sequences) # OLD code
                Feature.new_new_expand_all_features(args, contig_sequences, h_groups)
                
                # Print the set of new features
                self.logger.info('Feature.features')
                for f_name, f in sorted(Feature.features.items()):
                    self.logger.info("  {}:{}:{}:{}..{} PARENT={}".format(f.name, f.contig, f.strand, f.start, f.end, f.get_parent().name))
        
        # TODO: Will need to rewrite this filtering step BECAUSE some original (non-derived)
        #       Features will need to be filtered out
        # Get list of features that have compatible homology regions (Takes WAAAAYYYY too long)
        #sim_features = Feature.match_features_by_homology(args, contig_sequences) # Old code
        
        # Remove any features whose homology regions violate command-line argument '--homology_distance'
        #Feature.filter_features([x.name for x in sim_features])
        
        
        # Print to log what the filtered Features looks like
        self.logger.info('Feature.features after filtering')
        self.logger.info("  {}:{}:{}:{}..{} {} {} {}".format('Name', 'Contig', 'Strand', 'Start', 'End', 'Gene', 'Parent', 'Homologs'))
        for i, (f_name, f) in enumerate(sorted(Feature.features.items())):
            self.logger.info("  {} {}:{}:{}:{}..{} {} {} {}".format(i, f.name, f.contig, f.strand, f.start, f.end, f.get_gene(), f.get_parent().name, f.homologs))
        
        if (len(Feature.features) == 0):
            print("No Features match homology similarity requirements", file=sys.stderr)
            sys.exit(1)
        
        if (args.ko_gRNA):
            # Create 'Target' objects by searching for motif matches within 'Feature' sequences
            ExcisionTarget.search_all_features(args, contig_sequences)
            
            # Write the query list to FASTA
            # TODO: Make a method 'ExcisionTarget.generate_query(args, filename)' that automatically chooses
            #       if it will create a 'fasta' or 'fastq' file, depending on which Aligner is selected.
            if (args.selected_aligner.input == 'fasta'):
                ex_query_file = ExcisionTarget.generate_query_fasta(os.path.join(args.folder, 'excision-query.fasta'))
            elif (args.selected_aligner.input == 'fastq'):
                ex_query_file = ExcisionTarget.generate_query_fastq(os.path.join(args.folder, 'excision-query.fastq'))
        
        # Generate excision dDNAs and their associated reversion gRNA spacers
        # TODO: This should only be run if '--ko-dDNA' is specified
        ExcisionDonor.generate_donors(args, contig_sequences, feature2gene)
        
        # Design gRNAs to target the ko-dDNA.
        if (args.ki_gRNA):
            ReversionTarget.create_target_objects()
        ex_dDNA_file = ExcisionDonor.generate_fasta(os.path.join(args.folder, 'excision-dDNAs.fasta')) # TODO: Program will fail with error if this file is empty...
        
        if (args.ki_gRNA):
            # TODO: See previous TODO--make automatic choice depending on which Aligner is used.
            if (args.selected_aligner.input == 'fasta'):
                re_query_file = ReversionTarget.generate_query_fasta(os.path.join(args.folder, 'reversion-query.fasta'))
            elif (args.selected_aligner.input == 'fastq'):
                re_query_file = ReversionTarget.generate_query_fastq(os.path.join(args.folder, 'reversion-query.fastq'))
        
        if (args.ki_dDNA != None): # args.ki_dDNA can take one of 3 possible values: (None/True/'*.fasta')
            # Generate reversion dDNAs and write them to FASTA
            # TODO: Homology search for US/DS regions of excision-dDNA currently works.
            #       Now, you need to add a homology check for the reversion-dDNA homology regions
            ReversionDonor.generate_donors(args, contig_sequences, feature2gene)
            re_dDNA_file = ReversionDonor.generate_fasta(os.path.join(args.folder, 'reversion-dDNAs.fasta'))
        
        # Merge input FASTA files into a single one
        genome_fasta_file = utils.write_merged_fasta(contig_sequences, os.path.join(args.folder, 'genome.fasta'))
        
        # Index args.fasta for alignment
        #index_file = index_reference(args)
        genome_index_file = args.selected_aligner.index(genome_fasta_file, self.pathstrip(genome_fasta_file), args.folder, args.processors)
        ex_dDNA_index_file = args.selected_aligner.index(ex_dDNA_file, self.pathstrip(ex_dDNA_file), args.folder, args.processors)
        if (args.ki_dDNA == True):
            re_dDNA_index_file = args.selected_aligner.index(re_dDNA_file, self.pathstrip(re_dDNA_file), args.folder, args.processors)
        
        if (args.ko_gRNA):
            # Use selected alignment program to find all matches in the genome and dDNAs
            #ex_genome_align_file = align(ex_query_file, genome_index_file, args)
            exq2gDNA_align_file = args.selected_aligner.align(ex_query_file, genome_index_file, 'excision-query-2-gDNA', args.folder, args.processors)
            exq2exdDNA_align_file = args.selected_aligner.align(ex_query_file, ex_dDNA_index_file, 'excision-query-2-excision-dDNA', args.folder, args.processors)
            
            #print("ExcisionTarget before SAM parsing")
            #for et_seq, et_obj in ExcisionTarget.sequences.items():
            #    print(et_obj)
            
            # Load the SAM files and add Alignments to ExcisionTarget sequences
            ExcisionTarget.load_alignment(exq2gDNA_align_file, args, contig_sequences)
            ExcisionTarget.load_alignment(exq2exdDNA_align_file, args, ExcisionDonor.get_contig_dict())
            
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
        excision_spacers_file = ExcisionTarget.generate_spacers_fasta(os.path.join(args.folder, 'excision-targets.fasta'))
        
        # Use selected alignment program to find all matches in the genome and dDNAs
        #re_align_file = align(re_query_file, genome_index_file, args)
        if (args.ki_gRNA):
            req2gDNA_align_file = args.selected_aligner.align(re_query_file, genome_index_file, 'reversion-query-2-gDNA', args.folder, args.processors)
            req2exdDNA_align_file = args.selected_aligner.align(re_query_file, ex_dDNA_index_file, 'reversion-query-2-excision-dDNA', args.folder, args.processors)
        if (args.ki_dDNA == True):
            req2redDNA_align_file = args.selected_aligner.align(re_query_file, re_dDNA_index_file, 'reversion-query-2-reversion-dDNA', args.folder, args.processors)
        
        # Load the SAM files and add Alignments to ReversionTarget sequences
        if (args.ki_gRNA):
            ReversionTarget.load_alignment(req2gDNA_align_file, args, contig_sequences)
            ReversionTarget.load_alignment(req2exdDNA_align_file, args, ExcisionDonor.get_contig_dict())
        if (args.ki_dDNA == True):
            ReversionTarget.load_alignment(req2redDNA_align_file, args, ReversionDonor.get_contig_dict())
        
        
        if (args.ki_gRNA):
            # Calculate off-target/guide scores for each algorithm
            self.logger.info("ReversionTarget after SAM parsing and off-target scoring")
            # Somehow need to prioritize sequences for scoring.
            # Sequences with best diversity should be scored first.
            # If calculation time runs out, then just stop.
            # Should time each feature separately?
            ##### Subroutine #####
            #start_time = time.time()
            #if (time.time() - start_time >= args.max_time):
            #    self.logger.info('Site search terminated due to time constraints.')
            #    return
            ##### Subroutine #####
            for re_seq, re_obj in ReversionTarget.sequences.items():
                re_obj.score_off_targets(args, homologs)
                self.logger.info(re_obj)
                for a in re_obj.alignments:
                    self.logger.info('  ' + str(a))
            
            # Batch calculate with new ReversionTarget class
            # TODO: Ideally, all scoring Algorithms will be run in batch
            ReversionTarget.score_batch(args)
            
            self.logger.info("ReversionTarget after Azimuth calculation")
            for rt_seq, rt_obj in ReversionTarget.sequences.items():
                self.logger.info(rt_obj)
            
            # Generate the FASTA with the final scores
            reversion_spacers_file = ReversionTarget.generate_spacers_fasta(os.path.join(args.folder, 'reversion-targets.fasta'))
        
        # Test code to generate alignments
        ExcisionDonor.generate_alignments()
        
        # Pick out the best ones and print them out
        self.log_results(args, homologs, n=5)
        
        # If '--ki-gRNA' was on command line, then print to STDOUT
        if args.ki_gRNA:
            self.print_reTarget_results(args, homologs, feature2gene)
        
        # If '--ko-gRNA' was on command line, then print to STDOUT
        if args.ko_gRNA:
            self.print_exTarget_results(args, homologs, feature2gene)
        
        self.print_reDonor_results(args)
        self.print_exDonor_results(args)
        
        # If reversion AmpF/AmpR primers were desired, then print them to STDOUT
        if args.revert_amplification_primers:
            self.print_AmpFR_results(args)
        
        #print('======')
        #self.get_best_table(args, homologs, feature2gene)
        
        # Existing pipeline - Want to mutate feature, and cut site within feature
        # ,,,,,,,,,,,,,,,,NNNNNNNNNNNNNNNNNNNN111NNN................... Feature to be mutated contains cut site 1
        # ,,,,,,,,,,,,,,,,aaa222aaa.................................... Feature is knocked-out, and unique cut site 2 added
        # ,,,,,,,,,,,,,,,,mmmmmmmmmmmmmmmmmmmmmmmmmm................... Cut site 2 targeted and replaced by mutant
        
        # New pipeline - Want to mutate feature, but cut site outside of feature
        # ,,,,,,,,,,,,,,,,NNNNNNNNNNNNNN---------------111---.......... Feature to be mutated does not contain closest cut site 1
        # ,,,,,,,,,,,,,,,,nnnnnnnnnnnnnnnnnnnnnnnnnnnnn111nnn.......... Feature expanded to include cut site
        # ,,,,,,,,,,,,,,,,aaa222aaa.................................... Feature is knocked-out, and unique cut site 2 added
        # ,,,,,,,,,,,,,,,,mmmmmmmmmmmmmm---------------111---.......... Cut site 2 targeted and replaced by mutant (input) plus wild type cut site 1
        
        # 1) input: feature start and end
        # 2) if closest/best cut site is outside the feature
        #    Then expand the feature to include the closest/best cut site
        #    Also, concatenate the wild-type version of the expanded site to the mutant
        #    to create the new ki-dDNA
        
        # What if multiple features are adjacent to each other?
        # Then they should be combined into a single transformation.
        # ,,,,,,,,,NNNNNNNNNNN111NNN........NNNNNNNNN111NNN............ Since features 1 and 2 are independent, they can be kept as separate reactions
        # ,,,,,,,,,NNNNNNNNNN--------111-------NNNNNNNNNNNN............ Since features 1 and 2 both have the closest/best cut site, they are combined
        # ,,,,,,,,,nnnnnnnnnnnnnnnnnn111nnnnnnnnnnnnnnnnnnn............
        # ,,,,,,,,,aaa222aaa...........................................
        # ,,,,,,,,,mmmmmmmmmm--------111-------mmmmmmmmmmmm............ the dDNA is expanded to include the wt region in between both of these features
        
        # End 'compute()'
    
    @classmethod
    def get_reTarget_homologs(cls, homologs):
        """Get ReversionTarget objects for each homologous feature group"""
        ret_dict2 = {}
        
        if homologs:
            groups = set() # A set of all the parent (non-derived) feature names
            for feature_name, f in Feature.features.items():
                #contig, start, end, strand = features[feature]
                # Convert to tuple
                #groups.add(tuple(sorted(homologs[feature_name])))
                groups.add(tuple(sorted(homologs[f.get_parent().name])))
            # groups = {('F1_A', 'F1_B'), ('F2_A', 'F2_B')} # A set of tuples
            
            for g in groups:
                ret_dict2[g] = set()
                # Get all ReversionTargets that have all feature of this group in its location information
                
                for name, obj in ReversionTarget.indices.items():
                    #obj_features = set(x[0] for x in obj.locations) # these are comma-separated
                    obj_features = set(utils.flatten(x[0].split(',') for x in obj.locations))
                    
                    if (len(obj_features.intersection(g)) == len(g)):
                        ret_dict2[g].add(obj) # Store with the tuple form of the group as key
                
                #print('ret_dict2', len(g), g)
                #for obj in sorted(ret_dict2[g], key=lambda x: int(x.name.split('-')[1])):
                #    print(' ', obj)
        
        # key = ('F1_A', 'F1_B') # tuple of homologous features
        # value = [ReversionTarget(), ReversionTarget(), ...] # list of ReversionTarget objects
        return ret_dict2
    
    @classmethod
    def rank_donors(cls, donor_list):
        """Uses a heuristic to find the best Donor in the list"""
        # If the best ReversionTarget has multiple locations, then
        # choose the best ExcisionDonor location:
        #    1) minimize the distance between x and y: w..x:mAT:y..z
        #    2) minimize the length of mAT
        #    3) report all ties (don't break ties)
        
        rank_list = []
        for d in donor_list: # these should all be ExcisionDonor objects
            gaps = []
            for l in d.locations:
                gaps.append(len(l[4]) + l[5][0]-l[3][1])
            rank_list.append((min(gaps), d))
        
        return sorted(rank_list, key=lambda x: x[0]) # smallest insert size/gap length will be first
    
    @classmethod
    def rank_targets(cls, args, target_list):
        """Uses a heuristic to find the best Target in the list"""
        
        # Hsu-Zhang off-target score should be >95
        #   if Hsu-Zhange < 95, score should drop quickly
        # Azimuth on-target score should be >60
        #   if Azimuth < 60, score should drop quickly
        # CFD off-target score should be >50
        #   if CFD < 50, score should drop quickly
        
        #rank_azimuth = lambda x: 1/(1+1.17**(50-x))
        #rank_hsuzhang = lambda x: 1/(1+1.8**(90-x))
        #rank_cfd = lambda x: 1/(1+1.2**(40-x))
        #rank = lambda x, y, z: rank_azimuth(x)*rank_hsuzhang(y)*rank_cfd(z)
        
        # Returns list of targets, sorted such that the one with the highest
        # aggregate score is the 0th index
        rank_list = []
        for t in target_list:
            #rank_list.append((rank(t.score['Azimuth'], t.off_targets['Hsu-Zhang'], t.off_targets['CFD']), t))
            
            rank = 1.0
            for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms:
                if C.name in args.offtargetfilters + args.ontargetfilters:
                    if C.off_target:
                        rank *= C.weight(t.off_targets[C.name])
                    else:
                        rank *= C.weight(t.score[C.name])
            rank_list.append((rank, t))
        
        #return sorted(target_list, key=lambda x: rank(x.score['Azimuth'], x.off_targets['Hsu-Zhang'], x.off_targets['CFD']), reverse=True)
        return sorted(rank_list, key=lambda x: x[0], reverse=True)
    
    def print_reTarget_results(self, args, homologs, feature2gene):
        """
        Print to STDOUT the final output for the '_generate()' function.
        Lists the best 'reTarget' and their corresponding 'exDonor' objects.
        """
        # TODO: To my knowledge, this doesn't respect EITHER Feature.homologs when printing OR Target equivalents
        #       So this needs to be modified so that it DOES (somehow)...
        
        print('# reTarget results')
        header = ['# gene', 'features', 'weight', 'reTarget name', 'reTarget sequence'] + ['OT:{}'.format(x) for x in args.offtargetfilters] + [x for x in args.ontargetfilters] + ['exDonors', 'us-trim:mAT:ds-trim', 'feature:contig:strand:start..end', 'warnings']
        print('\t'.join(header))
        
        results = []
        # for weight, obj in self.rank_targets(ret_dict2[feature_homologs]):
        #for rt in ReversionTarget.indices.items()
        for weight, rt in self.rank_targets(args, ReversionTarget.indices.values()):
            
            genes = set()
            features = set()
            positions = []
            
            for l in rt.locations:
                # (feature, contig, orientation, start, end, upstream, downstream)
                for ll in l[0].split(','):
                    genes.add(Feature.get_gene_from_feature(ll, feature2gene))
                    features.add(Feature.features[ll].get_parent().name)
                    positions.append(rt.format_location(l))
            
            
            # Get the ExcisionDonor objects for this ki-spacer, and weigh them
            ex_donors = self.rank_donors(rt.get_donors())
            join_code = sorted(set(utils.flatten(map(lambda x: x[1].get_inserts_and_trims(), ex_donors))))
            
            genes = ','.join(sorted(genes))
            features = ','.join(sorted(features))
            positions = ','.join(positions)
            ex_donors = ','.join([ x[1].name for x in sorted(ex_donors, key=lambda y: int(y[1].name.split('-')[1])) ])
            join_code = ','.join('{}:{}:{}'.format(x[1], x[0], x[2]) for x in join_code)
            
            #othz = round(rt.off_targets['Hsu-Zhang'], 2)
            #otcfd = round(rt.off_targets['CFD'], 2)
            #azimuth = round(rt.score['Azimuth'], 2)
            
            off_scores = [round(rt.off_targets[x], 2) for x in args.offtargetfilters]
            on_scores = [round(rt.score[x], 2) for x in args.ontargetfilters]
            
            warnings = 'None'
            
            results.append([genes, features, weight, rt.name, rt.format_sequence()] + off_scores + on_scores + [ex_donors, join_code, positions, warnings])
        
        for r in results:
            print('\t'.join(map(str, r)))
    
    def print_exTarget_results(self, args, homologs, feature2gene):
        """
        Print to STDOUT the final output for the '_generate()' function.
        Lists the best 'exTarget' and their corresponding 'reDonor' objects.
        """
        
        print('# exTarget results')
        
        # If '--ki-dDNA' is specified, then cPCR searches will be performed
        # for each feature (parent and derived).
        # Based on the primer pair weights, AND the 'exTarget' weights, the
        # possible 'reDonor's can be given weights, which can be used to
        # rank them
        
        # TODO: Make it so the 'reDonors' and 'feature:contig:strand:start..end' just print the SMALLEST
        #       features that contain this Target (if the Features are derived)
        #       As it is now, the huge list of derived features isn't useful.
        
        header = ['# gene', 'features', 'weight', 'exTarget name', 'exTarget sequence'] + ['OT:{}'.format(x) for x in args.offtargetfilters] + [x for x in args.ontargetfilters] + ['reDonors', 'None', 'feature:contig:strand:start..end', 'warnings']
        print('\t'.join(header))
        
        results = []
        
        for weight, et in self.rank_targets(args, ExcisionTarget.indices.values()):
        
        #for feature in sorted(Feature.features):
        #    et_list = []
        #    for name, et in ExcisionTarget.indices.items():
        #        if feature in et.get_location_features():
        #            et_list.append(et)
        #    
        #    #rd_list = []
        #    #for name, rd in ReversionDonor.indices.items():
        #    #    if feature in rd.get_location_features():
        #    #        rd_list.append(rd)
        #    self.logger.info('---feature: {}'.format(feature))
        #    self.logger.info('---et_list: {}'.format(et_list))
        #    #self.logger.info('---rd_list: {}'.format(rd_list))
        #    
        #    # Print the top 20
        #    for weight, et in self.rank_targets(et_list)[:20]:
        #    self.logger.info('---weight, et: {}, {}'.format(weight, et))
            
            genes = set()
            features = set()
            positions = []
            
            for l in et.locations:
                # (feature, contig, orientation, start, end, upstream, downstream)
                for ll in l[0].split(','):
                    genes.add(Feature.get_gene_from_feature(ll, feature2gene))
                    features.add(Feature.features[ll].get_parent().name)
                    positions.append(et.format_location(l))
            
            
            # Get the ReversionDonor objects for this ko-gRNA, and weigh them
            #re_donors = self.rank_donors(et.get_donors())
            #re_donors = rd_list
            re_donors = set()
            for name, rd in ReversionDonor.indices.items():
                obj_features = set(rd.get_features())
                if (len(obj_features.intersection(et.get_features())) > 0):
                    re_donors.add(rd)
            #re_donors = sorted(re_donors, key=lambda x: int(x.name.split('-')[1]))
            #re_donors = ','.join(x.name for x in re_donors)
            re_donors = ','.join([ x.name for x in sorted(re_donors, key=lambda y: int(y.name.split('-')[1])) ])
            
            
            genes = ','.join(sorted(genes))
            features = ','.join(sorted(features))
            positions = ','.join(positions)
            #re_donors = ','.join([ x[1].name for x in sorted(re_donors, key=int(x[1].name.split('-')[1])) ])
            
            join_code = 'None'
            
            #othz = round(et.off_targets['Hsu-Zhang'], 2)
            #otcfd = round(et.off_targets['CFD'], 2)
            #azimuth = round(et.score['Azimuth'], 2)
            
            off_scores = [round(et.off_targets[x], 2) for x in args.offtargetfilters]
            on_scores = [round(et.score[x], 2) for x in args.ontargetfilters]
            
            warnings = 'None'
            
            results.append([genes, features, weight, et.name, et.format_sequence()] + off_scores + on_scores + [re_donors, join_code, positions, warnings])
        
        for r in results:
            print('\t'.join(map(str, r)))
    
    def print_reDonor_results(self, args):
        '''
        Print to STDOUT the top-weighted reDonor objects
        :param args: Argparse Namespace object containing commandline parameters
        :return: None
        '''
        # TODO: Finish this function that prints reDonors to STDOUT
        print('# reDonor results')
        headers = ['gene', 'features', 'weight', 'reDonor name', 'exTargets', 'feature:contig:strand:start..end', 'sequence']
        
        # After struggling with the current, less-than-clear output AddTag returns, here are a few ideas of how they can be improved:
        #   # Table of exTargets
        #   gene	features	weight	exTarget name	exTarget sequence	SCORE	SCORE	SCORE	best reDonor	best (e)Feature	best exDonor	warnings
        #   
        #   # Table of designs
        #   gene	features	reTarget	exTarget	exDonor	reDonor
    
    def print_exDonor_results(self, args):
        '''
        Print to STDOUT the top-weighted exDonor objects
        :param args: Argparse Namespace object containing commandline parameters
        :return: None
        '''
        # TODO: Finish this function that prints exDonors to STDOUT
        print('# exDonor results')
        # TODO: Need to invent a weight for a dDNA. Components:
        #        * weight of best reTarget
        #        * amount of extra DNA expanded from original Feature
        #        * polymorphisms in homology regions (should be boolean)
        headers = ['gene', 'features', 'weight', 'exDonor name', 'reTargets', 'feature:contig:strand:start..end', 'sequence']
    
    def print_AmpFR_results(self, args):
        '''
        Print to STDOUT the top-weighted AmpF/AmpR PrimerPairs
        :param args: Argparse Namespace object containing commandline parameters
        :return: None
        '''
        header = ['gene', 'features', 'weight', 'reDonor name', 'AmpF', 'AmpR', 'amplicon sizes', 'Tms', 'GCs', 'min(delta-G)']
        print('\t'.join(header))
        # TODO: Finish this function that prints AmpF/AmpR results to STDOUT
        #       For previous implementation, see lines around 1050~1100 in 'donors.py'
    
    @classmethod
    def log_results(cls, args, homologs, n=None):
        """Function that prints to the log file the best spacers and dDNAs for each feature"""
        
        cls.logger.info('Log of best results...')
        
        # Print best ReversionTargets calculated and their corresponding ExcisionDonors
        cls.logger.info("Best 'ReversionTarget's calculated and their corresponding 'ExcisionDonor's...")
        ret_dict2 = cls.get_reTarget_homologs(homologs)
        for k in ret_dict2:
            cls.logger.info(str(k) + ' ' + str(len(ret_dict2[k])))
            
            # Get the value for N
            if (n == None):
                display_num = len(ret_dict2[k])
            else:
                display_num = n
            
            # Print the top N
            for rank, obj in cls.rank_targets(args, ret_dict2[k])[:display_num]:
                cls.logger.info(' ' + str(rank) + ' ' + str(obj))
                # Get the ExcisionDonor objects for this ki-spacer, and rank them
                rds = cls.rank_donors(obj.get_donors())
                # filter out all but the top-ranked ones
                rds = [x for x in rds if (x[0] == rds[0][0])]
                for gap, exd_obj in rds:
                    cls.logger.info('   ' + str(gap) + ' ' + str(exd_obj.get_trims()) + ' ' + str(exd_obj))
        
        # Print best ExcisionTargets (not necessarily homozygous) for each feature
        cls.logger.info("Best 'ExcisionTarget's for each feature...")
        # and the ReversionDonor
        for feature in sorted(Feature.features):
            cls.logger.info(feature)
            et_list = []
            for name, obj in ExcisionTarget.indices.items():
                if feature in obj.get_location_features():
                    et_list.append(obj)
            
            # Get the value for N
            if (n == None):
                display_num = len(et_list[k])
            else:
                display_num = n
            
            # Print the top N
            for rank, obj in cls.rank_targets(args, et_list)[:display_num]:
                cls.logger.info('  ' + str(rank) + ' ' + str(obj))
            
            red_list = []
            for name, obj in ReversionDonor.indices.items():
                if feature in obj.get_location_features():
                    red_list.append(obj)
            for obj in red_list: # In no particular order (ampF/ampR primers for reDonors DO have weights, however)
                cls.logger.info('  ' + str(obj))
    
    def rolling_picker(self):
        """
        Rolling picker algorithm for multiplex gRNA/dDNA.
        
        Once all candidate target sequences are fully annotated and ranked, the
        sgRNA designer cycles through the list of candidates, attempting to pick
        sequences in order to achieve the best final set of sequences.
        
        We select all best-ranked gRNA, and if there is a conflict, then the
        worse one is replaced with a different candidate.
        
        After each round of picking, constraints are somewhat relaxed until
        a usable set of gRNA are found.
        
        """
        pass
    
    ###### DEPRECATED ######
    def merge_features(self, features):
        """Combine overlapping features?"""
        return features
    ###### DEPRECATED ######
    
    ######################## Start DEPRECATED ########################

    def get_exTarget_homologs(self, homologs):
        """Get ExcisionTarget objects for each homologous feature group"""
        ext_dict = {}
        if homologs:
            groups = set()
            for feature_name, f in Feature.features.items():
                #groups.add(tuple(sorted(homologs[feature_name])))
                groups.add(tuple(sorted(homologs[f.get_parent().name])))
            for g in groups:
                ext_dict[g] = set()
                for name, obj in ExcisionTarget.indices.items():
                    obj_features = set(utils.flatten(x[0].split(',') for x in obj.locations))
                    if (len(obj_features.intersection(g)) == len(g)):
                        ext_dict[g].add(obj)
        return ext_dict
    
    def get_exTarget_allele_specific(self):
        """Gets allele-specific ExcisionTarget objects"""
        ext_dict = {}
        for feature_name, f in Feature.features.items():
        #for f in features:
            ext_dict[feature_name] = set()
            for name, obj in ExcisionTarget.indices.items():
                obj_features = set(utils.flatten(x[0].split(',') for x in obj.locations))
                if ((len(obj_features) == 1) and (feature_name in obj_features)):
                    ext_dict[feature_name].add(obj) # Store with the feature as key
        return ext_dict
    
    def get_reTarget_allele(self, features):
        """Gets all ReversionTarget objects for each feature"""
        ret_dict = {}
        for f in features:
            ret_dict[f] = set()
            for name, obj in ReversionTarget.indices.items():
                obj_features = set(utils.flatten(x[0].split(',') for x in obj.locations))
                if (f in obj_features):
                    ret_dict[f].add(obj) # Store with the feature as key
        return ret_dict
    
    def get_reTarget_allele_specific(self):
        """Gets allele-specific ReversionTarget objects"""
        ret_dict = {}
        for feature_name, f in Feature.features.items():
        #for f in features:
            ret_dict[feature_name] = set()
            for name, obj in ReversionTarget.indices.items():
                obj_features = set(utils.flatten(x[0].split(',') for x in obj.locations))
                if ((len(obj_features) == 1) and (feature_name in obj_features)):
                    ret_dict[feature_name].add(obj) # Store with the feature as key
        return ret_dict
    
    def get_best_table(self, args, homologs, feature2gene):
        """
        Identify and print the best spacers and dDNAs for each feature, given
        each mAT insert size, and us/ds trim length.
        Thus, the user can see the best spacer for any given combination of these.
        """
        
        ########## Results table for the knock-out dDNAs and the efficiency of their gRNA for cutting them (ExcisionDonor+ReversionTarget or ko-dDNA+ki-gRNA) ##########
        # Add a column for "WARNING", that tells the user if the target feature overlaps with another feature
        #  it will just list the feature IDs that overlap:
        #   C4_03620C_A-T,C4_03620C_A-T-E1
        
        # Only include columns for algorithms whose weight does not equal 1 uniformly (alternatively, whose classes override the default 'weight()' method)
        #algorithms.weighted_algorithms
        
        #header = ['gene', 'features', 'insert', 'mAT', 'translations', '(us, ds) trim', 'weight', 'OT:Hsu-Zhang', 'OT:CFD', 'Azimuth', 'reTarget name', 'reTarget sequence', 'ExDonors']
        # Add columns for:
        #  contig:start..end
        #   contig_A:40221..40241,contig_B:40155..40175
        #  feature:start..end
        #   CR_02630C_A:221..241,CR_02630C_B:224..244
        #  feature,feature:contig:strand:start..end
        #   C1_06280C_A,C1_06280C_B:exDonor-222:-:43..66
        #  strand
        #   +     or     -
        #
        # 'hairpin ΔG', 'homodimer ΔG', 'heterodimer ΔG'
        # locations0s: 'features:contig:strand:start..end'
        header = ['gene', 'features', 'contig:strand:start..end', 'us-trim:mAT:ds-trim', 'translations', 'weight', 'OT:Hsu-Zhang', 'OT:CFD', 'Azimuth', 'reTarget name', 'reTarget sequence', 'exDonors', 'warning']
        print('\t'.join(header))
        
        # Get best ReversionTargets by calculating their weights, and also getting
        # the ExcisionDonors that correspond to the ReversionTarget
        ret_dict2 = self.get_reTarget_homologs(homologs)
        for feature_homologs in sorted(ret_dict2):
            # Get these two things which are common to all the top hits
            gene = feature2gene[feature_homologs[0]] # Get the gene name
            csfeatures = ','.join(feature_homologs) # Add the features as comma-separated list
            
            features_pos = '' #features[gene]
            
            outputs = {}
            
            # Print the top N for each insert size and trim
            for weight, obj in self.rank_targets(args, ret_dict2[feature_homologs]):
                othz = round(obj.off_targets['Hsu-Zhang'], 2)
                otcfd = round(obj.off_targets['CFD'], 2)
                azimuth = round(obj.score['Azimuth'], 2)
                
                # Get the ExcisionDonor objects for this ki-spacer, and weigh them
                rds = self.rank_donors(obj.get_donors())
                # filter out all but the top-weighted ones
                rds = [x for x in rds if (x[0] == rds[0][0])]
                
                exdonors = ','.join(map(lambda x: x[1].name, rds))
                
                key0 = sorted(set(utils.flatten(map(lambda x: x[1].get_inserts_and_trims(), rds))))
                key0s = ','.join('{}:{}:{}'.format(x[1], x[0], x[2]) for x in key0)
                key1 = sorted(set([(len(x[0]), x[1], x[2]) for x in key0])) # replace mAT with length
                key1s = ','.join(map(str, key1))
                #locations0 = sorted(set(utils.flatten([xo[1].format_location(x) for x in xo[1].locations] for xo in rds)))
                #locations0s= ','.join(locations0)
                translations = None
                
                sline = [gene, csfeatures, features_pos, key0s, translations, weight, othz, otcfd, azimuth, obj.name, obj.spacer+'|'+obj.pam, exdonors]
                
                if (weight >= args.min_weight_reported):
                    if (len(outputs.get(key1s, [])) < args.max_number_sequences_reported):
                        outputs.setdefault(key1s, []).append(sline)
            
            for k in sorted(outputs): # sort by insert/trims
            #for k in sorted(outputs, key=lambda x: outputs[x][4], reverse=True): # sort by weight
                for sline in outputs[k]:
                    print('\t'.join(map(str, sline)))
            
            if (len(outputs) == 0):
                self.logger.info('No spacers for targeting knocked-out ' + csfeatures)
        
        # Print a table of allele-specific knock-in spacer and dDNA pairs
        ret_dict = self.get_reTarget_allele_specific()
        for feature_name in sorted(ret_dict):
            #gene = feature2gene[feature] # Get the gene name
            gene = self.get_gene_from_feature(feature_name, feature2gene)
            
            features_pos = ''
            
            outputs = {}
            
            for weight, obj in self.rank_targets(args, ret_dict[feature_name]):
                othz = round(obj.off_targets['Hsu-Zhang'], 2)
                otcfd = round(obj.off_targets['CFD'], 2)
                azimuth = round(obj.score['Azimuth'], 2)
                
                rds = self.rank_donors(obj.get_donors())
                rds = [x for x in rds if (x[0] == rds[0][0])]
                
                exdonors = ','.join(map(lambda x: x[1].name, rds))
                
                key0 = sorted(set(utils.flatten(map(lambda x: x[1].get_inserts_and_trims(), rds))))
                key0s = ','.join('{}:{}:{}'.format(x[1], x[0], x[2]) for x in key0)
                key1 = sorted(set([(len(x[0]), x[1], x[2]) for x in key0])) # replace mAT with length
                key1s = ','.join(map(str, key1))
                translations = None
                
                sline = [gene, feature_name, features_pos, key0s, translations, weight, othz, otcfd, azimuth, obj.name, obj.spacer+'>'+obj.pam, exdonors] # should automatically determine the direction of the SPACER>PAM based on the motif
                
                if (weight >= args.min_weight_reported):
                    if (len(outputs.get(key1s, [])) < args.max_number_sequences_reported):
                        outputs.setdefault(key1s, []).append(sline)
            
            for k in sorted(outputs):
                for sline in outputs[k]:
                    print('\t'.join(map(str, sline)))
            
            # If there are no allele-specific records, then nothing is printed
            if (len(outputs) == 0):
                self.logger.info('No allele-specific spacers for targeting knocked-out ' + feature_name)
        
        ########## Results table for cutting the wild-type genome (ExcisionTarget or ko-gRNA) ##########
        self.logger.info("Results table for cutting the wild-type genome ('ExcisionTarget' or 'ko-gRNA')")
        
        # to add to header: contig:strand:start..end
        header = ['gene', 'features', 'contig:strand:start..end', 'weight', 'OT:Hsu-Zhang', 'OT:CFD', 'Azimuth', 'exTarget name', 'exTarget sequence', 'reDonors']
        print('\t'.join(header))
        
        # Print the best homozygous ExcisionTargets for each feature set
        ext_dict2 = self.get_exTarget_homologs(homologs)
        self.logger.info("len(ext_dict2) = {}".format(len(ext_dict2)))
        
        for feature_homologs in sorted(ext_dict2):
            self.logger.info('feature_homologs = {}'.format(feature_homologs))
            gene = feature2gene[feature_homologs[0]] # Get the gene name
            csfeatures = ','.join(feature_homologs) # Add the features as comma-separated list
            
            red_list = set()
            for name, obj in ReversionDonor.indices.items():
                obj_features = set(obj.get_location_features())
                if (len(obj_features.intersection(feature_homologs)) > 0):
                    red_list.add(obj)
            red_list = sorted(red_list, key=lambda x: int(x.name.split('-')[1]))
            
            outputs = []
            for weight, obj in self.rank_targets(args, ext_dict2[feature_homologs]):
                othz = round(obj.off_targets['Hsu-Zhang'], 2)
                otcfd = round(obj.off_targets['CFD'], 2)
                azimuth = round(obj.score['Azimuth'], 2)
                
                redonors = ','.join(x.name for x in red_list)
                
                sline = [gene, csfeatures, '', weight, othz, otcfd, azimuth, obj.name, obj.spacer+'>'+obj.pam, redonors] # should automatically determine the direction of the SPACER>PAM based on the motif
                
                if (weight >= args.min_weight_reported):
                    if (len(outputs) < args.max_number_sequences_reported):
                        outputs.append(sline)
            
            for sline in outputs:
                print('\t'.join(map(str, sline)))
            
            if (len(outputs) == 0):
                self.logger.info('No spacers for targeting ' + csfeatures)
        
        # Print a table of allele-specific knock-out spacer and dDNA pairs
        ext_dict = self.get_exTarget_allele_specific()
        for feature_name in sorted(ext_dict):
            #gene = feature2gene[feature] # Get the gene name
            gene = self.get_gene_from_feature(feature_name, feature2gene)
            
            red_list = set()
            for name, obj in ReversionDonor.indices.items():
                if feature_name in obj.get_location_features():
                    red_list.add(obj)
            red_list = sorted(red_list, key=lambda x: int(x.name.split('-')[1]))
            
            outputs = []
            
            for weight, obj in self.rank_targets(args, ext_dict[feature_name]):
                othz = round(obj.off_targets['Hsu-Zhang'], 2)
                otcfd = round(obj.off_targets['CFD'], 2)
                azimuth = round(obj.score['Azimuth'], 2)
                
                redonors = ','.join(x.name for x in red_list)
                
                sline = [gene, feature_name, '', weight, othz, otcfd, azimuth, obj.name, obj.spacer+'|'+obj.pam, redonors]
                
                if (weight >= args.min_weight_reported):
                    if (len(outputs) < args.max_number_sequences_reported):
                        outputs.append(sline)
            
            for sline in outputs:
                print('\t'.join(map(str, sline)))
            
            # If there are no allele-specific records, then nothing is printed
            if (len(outputs) == 0):
                self.logger.info('No allele-specific spacers for targeting wild-type ' + feature_name)
    
    def pathstrip(self, path):
        return os.path.splitext(os.path.basename(path))[0]
    
    def old_get_best(self, args, features, contigs):
        #exd_dict = {}
        red_dict = {}
        ext_dict = {}
        ret_dict = {}
        
        # Non-exclusively separate all instances into feature groups
        for feature in features:
            red_dict[feature] = []
            for name, obj in ReversionDonor.indices.items():
                if feature in obj.get_location_features():
                    red_dict[feature].append(obj)
            
            ext_dict[feature] = []
            for name, obj in ExcisionTarget.indices.items():
                if feature in obj.get_location_features():
                    ext_dict[feature].append(obj)
            
            ret_dict[feature] = set()
            for name, obj in ReversionTarget.indices.items():
                for c in obj.get_contigs():
                    if feature in ExcisionDonor.indices[c].get_location_features():
                        ret_dict[feature].add(obj)
        
        # Find the best instances for ko/ki
        for feature in sorted(features):
            # Find the best ReversionTarget
            on_target_sorted = sorted(ret_dict[feature], key=lambda x: x.score['Azimuth'], reverse=True)
            ret_best = None
            for i in range(len(on_target_sorted)):
                if (on_target_sorted[i].off_targets['Hsu-Zhang'] >= 90):
                    ret_best = on_target_sorted[i]
                    break
            #if not ret_best:
            #    print(len(on_target_sorted))
            #    for tmp in on_target_sorted:
            #        print(' ', tmp)
            #if not ret_best: # this could fail if there are no sites
            #    off_target_sorted = sorted(ret_dict[feature], key=lambda x: x.score['Hsu-Zhang'], reverse=True)
            #    ret_best = off_target_sorted[0]
            
            # Find the best ExcisionTarget
            on_target_sorted = sorted(ext_dict[feature], key=lambda x: x.score['Azimuth'], reverse=True)
            ext_best = None
            for i in range(len(on_target_sorted)):
                if (on_target_sorted[i].off_targets['Hsu-Zhang'] >= 95):
                    ext_best = on_target_sorted[i]
                    break
            
            # Find the ExcisionDonors that correspond with the best ReversionTarget
            exd_best = []
            if ret_best:
                exd_best = [ExcisionDonor.indices[x] for x in ret_best.get_contigs()]
            
            # Find the ReversionDonor
            red_best = red_dict[feature]
            
            self.logger.info("")
            self.logger.info("  feature = " + feature)
            self.logger.info("ko spacer = " + str(ext_best))
            for d in exd_best:
                self.logger.info("  ko dDNA = " + str(d))
            self.logger.info("ki spacer = " + str(ret_best))
            for d in red_best:
                self.logger.info("  ki dDNA = " + str(d))
    
    ######################## End DEPRECATED ########################
