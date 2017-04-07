#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/__init__.py

# Import standard packages
import sys
import os
import argparse
import textwrap
import time

# Import non-standard packages
import regex

# Import included AddTag-specific modules
from . import utils
from . import nucleotides
from . import scores
from . import hsuzhang
from . import housden
from . import morenomateos
from . import doench
from . import bowtie2

# Define meta variables
__author__ = "Thaddeus D. Seher (@tdseher) & Aaron Hernday"
__date__ = utils.load_git_date()
__fullversion__ = utils.load_git_version()
__version__ = __fullversion__[:7]
__commits__ = utils.load_git_commits()
__program__ = os.path.basename(sys.argv[0])
__description__ = """\
description:
  Program for identifying unique endogenous gRNA sites 
  and creating unique synthetic gRNA sites.
  
  Copyright (c) 2017 {__author__}.
  All rights reserved.

version:
  short   {__version__}
  full    {__fullversion__}
  commits {__commits__}
  date    {__date__}

""".format(**locals())
__epilog__ = """\
example:
 $ python3 {__program__} genome.fasta genome.gff > results.txt
""".format(**locals())

class Sequence(object):
    """Data structure defining a sequence"""
    
    substitution_threshold = 4
    insertion_threshold = 2
    deletion_threshold = 2
    error_threshold = 5
    
    doench2014_threshold = 1.0
    doench2016_threshold = 1.0
    hsuzhang_threshold = 1.0
    linear_threshold = 80.0
    morenomateos_threshold = 1.0
    
    def __init__(self, feature, contig_sequence, pams, contig=None, contig_start=None, contig_end=None, contig_orientation=None, feature_orientation=None):
        """Create a structure for holding individual sequence information"""
        self.feature = feature
        self.feature_orientation = feature_orientation
        
        self.contig_sequence = contig_sequence
        self.contig_target, self.contig_pam = nucleotides.split_target_sequence(self.contig_sequence, pams, force=True)
        #self.disambiguated_sequences = disambiguate_iupac(self.contig_sequence, kind="dna") # need to apply pre-filters
        self.contig = contig
        self.contig_start = contig_start
        self.contig_end = contig_end
        self.contig_orientation = contig_orientation
        
        # List to store alignments
        self.alignments = []
        
        # query sequence only
        self.gc = scores.gc_score(self.contig_target)
        self.doench2014 = doench.on_target_score_2014(self.contig_target, self.contig_pam, upstream='', downstream='')
        self.doench2016 = 100 # assume perfect match (should NOT assume if ambiguities in reference)
        self.hsuzhang = 100 # assume perfect match (should NOT assume if ambiguities in reference)
        self.linear = 100 # assume perfect match (should NOT assume if ambiguities in reference)
        self.housden = housden.housden_score(self.contig_target)
        self.morenomateos = morenomateos.morenomateos_score(self.contig_target, self.contig_pam, upstream='', downstream='')
        
        # query x subject score
        self.off_target_doench2014 = None
        self.off_target_doench2016 = None
        self.off_target_hsuzhang = None
        self.off_target_linear = None
        self.off_target_morenomateos = None
    
    def add_alignment(self, aligned_sequence, pams, aligned_contig, aligned_start, aligned_end, aligned_orientation):
        """Add a genomic position to the list of alignments"""
        aligned_target, aligned_pam = nucleotides.split_target_sequence(aligned_sequence, pams, force=True)
        substitutions, insertions, deletions = nucleotides.count_errors(self.contig_sequence, aligned_sequence)
        seq = (
            aligned_sequence,
            aligned_target,
            aligned_pam,
            aligned_contig,
            aligned_start,
            aligned_end,
            aligned_orientation,
            substitutions,
            insertions,
            deletions,
            doench.on_target_score_2014(aligned_target, aligned_pam),
            doench.on_target_score_2016(self.contig_target, aligned_target, aligned_pam),
            hsuzhang.hsuzhang_score(self.contig_target, aligned_target, iupac=False),
            scores.linear_score(self.contig_target, aligned_target),
            housden.housden_score(aligned_target),
            morenomateos.morenomateos_score(aligned_target, aligned_pam),
            nucleotides.ridentities(self.contig_target, aligned_target),
            scores.r_score(self.contig_target, aligned_target, 4),
            scores.r_score(self.contig_target, aligned_target, 8),
            scores.r_score(self.contig_target, aligned_target, 12),
            scores.r_score(self.contig_target, aligned_target, 16),
        )
        self.alignments.append(seq)
    
    def score(self):
        """Calculate Guide scores for each algorithm"""
        doench2014_list = []
        doench2016_list = []
        hsuzhang_list = []
        linear_list = []
        morenomateos_list = []
        for a in self.alignments:
            if (a[3:6] != (self.contig, self.contig_start, self.contig_end)):
                # lambda a, b: all([a[0]<=b[0], a[1]<=b[1], a[2]<b[2], sum(a)<=b[3]])
                if ((a[7] <= self.substitution_threshold) and
                    (a[8] <= self.insertion_threshold) and
                    (a[9] <= self.deletion_threshold) and 
                    (sum(a[7:10]) <= self.error_threshold)
                ):
                    if (a[10] >= self.doench2014_threshold):
                        doench2014_list.append(a[10])
                    if (a[11] >= self.doench2016_threshold):
                        doench2016_list.append(a[11])
                    if (a[12] >= self.hsuzhang_threshold):
                        hsuzhang_list.append(a[12])
                    if (a[13] >= self.linear_threshold):
                        linear_list.append(a[13])
                    if (a[15] >= self.morenomateos_threshold):
                        morenomateos_list.append(a[15])
        self.off_target_doench2014 = scores.off_target_score(doench2014_list, (self.doench2014,))
        self.off_target_doench2016 = scores.off_target_score(doench2016_list, (self.doench2016,))
        self.off_target_hsuzhang = scores.off_target_score(hsuzhang_list, (self.hsuzhang,))
        self.off_target_linear = scores.off_target_score(linear_list, (self.linear,))
        self.off_target_morenomateos = scores.off_target_score(morenomateos_list, (self.morenomateos,))
    
    def __repr__(self):
        """Return a string containing a printable representation of the Sequence object."""
        return 'Sequence(feature=' + self.feature + ', ' + \
            self.contig + ':' + str(self.contig_start) + ':' + \
            str(self.contig_end) + ', ' + self.contig_target + '|' + \
            self.contig_pam + ', alignments=' + str(len(self.alignments)) + ')'

class CustomHelpFormatter(argparse.HelpFormatter):
    """Help message formatter which retains any formatting in descriptions
    and adds default values to argument help.
    
    Only the name of this class is considered a public API. All the methods
    provided by the class are considered an implementation detail.
    """
    # This class combines:
    #   argparse.ArgumentDefaultsHelpFormatter
    #   argparse.RawDescriptionHelpFormatter
    
    def _fill_text(self, text, width, indent):
        return ''.join([indent + line for line in text.splitlines(True)])
    
    def _get_help_string(self, action):
        help = action.help
        if '%(default)' not in action.help:
            if action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += ' (default: %(default)s)'
        return help

def load_sam_file_test(filename, pams, contigs, sep=':'):
    """Read in SAM file.
    sep is the separator for the header. Positions are converted to 0-index
    Creates a list of Sequence objects
    """
    
    alignments = {}
    with open(filename, 'r') as flo:
        for line in flo:
            if not line.startswith('@'):
                sline = line.rstrip().split("\t")
                if (len(sline) > 5):
                    feature, source_contig, source_start, source_end = sline[0].split(sep)
                    source = (feature, source_contig, int(source_start), int(source_end))
                    
                    # Get orientation
                    alignment_orientation = utils.sam_orientation(int(sline[1]))
                    
                    # Get alignment position
                    alignment_contig = sline[2]
                    alignment_start = int(sline[3])-1
                    alignment_end = int(sline[3])-1+utils.cigar_length(sline[5])
                    
                    # Reverse-complement if needed
                    alignment_sequence = contigs[alignment_contig][alignment_start:alignment_end]
                    actual_sequence = sline[9]
                    if (alignment_orientation == '-'):
                        alignment_sequence = nucleotides.rc(alignment_sequence)
                        actual_sequence = nucleotides.rc(actual_sequence)
                    
                    # if source not in alignments:
                    #    alignments[source] = s
                    # alignments[sournce].add_alignment(...)
                    
                    # Assuming creating an instance of Sequence() is cheaper
                    # than traversing alignments dict()
                    s = Sequence(
                        feature,
                        actual_sequence,
                        # contigs[source_contig][int(source_start):int(source_end)], # contig_sequence
                        pams, # ['NGG']
                        contig=source_contig,
                        contig_start=int(source_start),
                        contig_end=int(source_end),
                        contig_orientation=None,
                        feature_orientation=None
                    )
                    
                    alignments.setdefault(source, s).add_alignment(
                        alignment_sequence, # aligned_sequence (as when matched with reference, thus may be revcomp of initial query)
                        pams,
                        alignment_contig, # aligned_contig
                        alignment_start, # aligned_start
                        alignment_end, # aligned_end
                        alignment_orientation # aligned_orientation (+/-)
                    )
    
    #return list(alignments.values()) # unsorted list
    return list(map(lambda x: alignments[x], sorted(alignments))) # sorted

def parse_arguments():
    # Create the argument parser
    parser = argparse.ArgumentParser(
        description=__description__,
        epilog=__epilog__,
        formatter_class=CustomHelpFormatter
    )
    
    # Add required arguments
    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument("--fasta", required=True, metavar="*.fasta", type=str,
        help="FASTA file with contigs to find unique gRNA sites. Input FASTA \
            must not be compressed. Ambiguous bases within the FASTA will not \
            be chosen for gRNA sites. All FASTA sequences should have unique \
            primary headers (everything between the '>' and the first ' ' \
            should be unique).")
    required_group.add_argument("--gff", required=True, metavar="*.gff", type=str,
        help="GFF file specifying chromosomal features")
    required_group.add_argument("--folder", required=True, metavar="FOLDER",
        type=str, help="Path of folder to store generated files")
    
    # Add optional arguments
    parser.add_argument("--test", action="store_true",
        help="Perform tests only")
    #parser.add_argument("--pams", metavar="SEQ", nargs="+", type=str,
    #    default=["NGG"], help="Constrain finding only targets with these PAM sites")
    #parser.add_argument("--target_lengths", nargs=2, metavar=('MIN', 'MAX'),
    #    type=int, default=[17, 20],
    #    help="The length range of the 'target'/'spacer'/gRNA site")
    # Replacement for --pams and --target_lengths with this:
    parser.add_argument("--motifs", metavar="MOTIF", nargs="+", type=str,
        default=["N{17,20}>NGG"],
        help="Find only targets with these 'SPACER>PAM' motifs, written from \
        5' to 3'. '>' points toward PAM. IUPAC ambiguities accepted. '{a,b}' \
        are quantifiers. Examples: 'G{,2}N{19,20}>NGG', 'N{8,10}AN{10}>NGA', \
        'N{20,21}>NNGRRT', 'TTTN<N{20,23}'")
    #    maybe this: 'GGN<N{20}|N{4-20}|N{20}>NGG' where '|' separates spacer from non-targets?
    #       or this: 'GGN<N{20}.{4-20}N{20}>NGG' where '.' is any non-spacer nucleotide?
    parser.add_argument("--tag", metavar='TAG', type=str, default='ID',
        help="GFF tag with feature names. Examples: 'ID', 'Name', 'Parent', or 'locus_tag'")
    #parser.add_argument("--feature_homolog_regex", metavar="REGEX", type=str, default=None, help="regular expression with capturing group containing invariant feature. Example: '(.*)_[AB]' will treat features C2_10010C_A and C2_10010C_B as homologs")
    # okay idea, but needs more thought before implementation
    parser.add_argument("--feature_homologs", metavar="*.homologs", type=str, default=None,
        help="Path to text file containing homologous features on the same \
            line, separated by TAB characters")
    #parser.add_argument("--min_contig_edge_distance", metavar="N", type=int, default=500,
    #    help="Minimum distance from contig edge a site can be found")
    parser.add_argument("--features", metavar="FEATURE", type=str, nargs="+", default=["gene"],
        help="Features to design gRNA sites against. Must exist in GFF file. Examples: 'CDS', 'gene', 'mRNA', 'exon'")
    parser.add_argument("--target_gc", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[25, 75],
        help="Generated gRNAs must have %%GC content between these values (excluding PAM motif)")
    parser.add_argument("--excise_donor_homology", nargs=2, metavar=("MIN", "MAX"), type=int, default=[40,80],
        help="Range of homology lengths acceptable for knock-out dDNAs, inclusive.")
    parser.add_argument("--excise_donor_lengths", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[90, 100],
        help="Range of lengths acceptable for knock-out dDNAs, inclusive.")
    parser.add_argument("--excise_insert_lengths", nargs=2, metavar=("MIN", "MAX"), type=int, default=[0,7],
        help="Range for inserted DNA lengths, inclusive (mini-AddTag, mAT). If MIN < 0, then regions of dDNA homology (outside the feature) will be removed.")
    parser.add_argument("--revert_donor_lengths", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[300, 600],
        help="Range of lengths acceptable for knock-in dDNAs.")
    parser.add_argument("--min_feature_edge_distance", metavar="MIN", type=int, default=23,
        help="The minimum distance a gRNA site can be from the edge of the \
             feature. If negative, the maximum distance a gRNA site can be \
             outside the feature.")
    #parser.add_argument("--min_donor_insertions", metavar="N", type=int, default=2,
    #    help="The uniqueness of final donor DNA compared to the rest of the genome")
    #parser.add_argument("--min_donor_deletions", metavar="N", type=int, default=2,
    #    help="The uniqueness of final donor DNA compared to the rest of the genome")
    #parser.add_argument("--min_donor_substitutions", metavar="N", type=int, default=2,
    #    help="The uniqueness of final donor DNA compared to the rest of the genome")
    #parser.add_argument("--min_donor_errors", metavar="N", type=int, default=3,
    #    help="The uniqueness of final donor DNA compared to the rest of the genome")
    parser.add_argument("--min_donor_distance", metavar="N", type=int, default=36,
        help="The minimum distance in bp a difference can exist from the edge of donor DNA") # homology with genome
    parser.add_argument("--max_consecutive_ts", metavar="N", type=int, default=4,
        help="The maximum number of Ts allowed in generated gRNA sequences")
    # program currently will only search 'both' strands
    #parser.add_argument("--strands", type=str, choices=["+", "-", "both"], default="both",
    #    help="Strands to search for gRNAs")
    parser.add_argument("--ambiguities", type=str, choices=["discard", "disambiguate", "keep"], default="discard",
        help="How generated gRNAs should treat ambiguous bases: \
        discard - no gRNAs will be created where the FASTA has an ambiguous base; \
        disambiguate - gRNAs containing ambiguous bases will be converted to a set of non-ambiguous gRNAs; \
        keep - gRNAs can have ambiguous bases")
    parser.add_argument("--case", type=str, default="ignore",
        choices=["ignore", "discard-lower", "discard-upper", "invariant-lower", "invariant-upper"],
        help="Restrict generation of gRNAs based on case of nucleotides in input FASTA")
    # Add command line arguments for the additional hard constraints:
    #  Only report potential targets that have no off targets with mismatches within 8, 12, N nt from 3' end
    parser.add_argument("--processors", metavar="N", type=int, default=(os.cpu_count() or 1),
        help="Number of processors to use when performing pairwise sequence alignments")
    parser.add_argument("--aligner", type=str, choices=['addtag', 'blast+',
        'blat', 'bowtie', 'bowtie2', 'bwa', 'cas-offinder'], default='bowtie2',
        help="Program to calculate pairwise alignments. Please note that the 'addtag' internal aligner is very slow.")
    # Other aligners to consider: 'rmap', 'maq', 'shrimp2', 'soap2', 'star', 'rhat', 'mrsfast', 'stampy'
    parser.add_argument("--python2_path", type=str, default="python",
        help="Path to the Python 2.7+ program")
    parser.add_argument("--bowtie_path", type=str, default="bowtie",
        help="Path to the 'bowtie' executable")
    parser.add_argument("--bowtie-build_path", type=str, default="bowtie-build",
        help="Path to the 'bowtie-build' executable")
    parser.add_argument("--bowtie2_path", type=str, default="bowtie2",
        help="Path to the 'bowtie2' executable")
    parser.add_argument("--bowtie2-build_path", type=str, default="bowtie2-build",
        help="Path to the 'bowtie2-build' executable")
    parser.add_argument("--bwa_path", type=str, default="bwa",
        help="Path to the 'bwa' executable")
    parser.add_argument("--blastn_path", type=str, default="blastn",
        help="Path to the 'blastn' executable")
    parser.add_argument("--blat_path", type=str, default="blat",
        help="Path to the 'blat' executable")
    
    
    # Special version action
    parser.add_argument("-v", "--version", action='version', version='{__program__} {__version__}'.format(**globals()))
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Parse the motifs
    args.parsed_motifs = []
    for motif in args.motifs:
        args.parsed_motifs.append(parse_motif(motif))
    
    # Return the parsed arguments
    return args

def _parse_motif_helper(submotif):
    # Keep track of expanded sequences
    sequences = ['']
    
    # Keep track if a quantifier is being parsed
    quantifier = None
    
    # Iterate through the characters
    for c in submotif:
        if (c == '{'): # Start the quantifier
            quantifier = ''
        elif (c == '}'): # End the quantifier
            quantifier_list = quantifier.split(',')
            if (len(quantifier_list) == 1):
                min_length = int(quantifier)
                max_length = min_length
                
            elif (len(quantifier_list) == 2):
                if (quantifier_list[0] == ''):
                    min_length = 0
                else:
                    min_length = int(quantifier_list[0])
                if (quantifier_list[1] == ''):
                    raise Exception("Motif quantifier '{" + quantifier + "}' contains no maximum value")
                else:
                    max_length = int(quantifier_list[1])
                if (min_length > max_length):
                    raise Exception("Motif quantifier '{" + quantifier + "}' minimum and maximum lengths are invalid")
            
            last_chars = [ x[-1] for x in sequences ]
            
            sequences = [ x[:-1] for x in sequences ]
            
            new_sequences = []
            for i, s in enumerate(sequences):
                for length in range(min_length, max_length+1):
                    new_sequences.append(s + last_chars[i]*length)
            sequences = new_sequences
            quantifier = None
        elif (quantifier != None): # add current character to quantifier if it is open
            quantifier += c
        else: # add the current character to the expanded sequences
            for i in range(len(sequences)):
                sequences[i] = sequences[i] + c
    
    return sequences

def parse_motif(motif):
    # Eventually, replace --pams and --target_lengths with this:
    #parser.add_argument("--motifs", metavar="SEQ", nargs="+", type=str,
    #    default=["N{17,20}>NGG"], help="Find only targets with these \
    #    'SPACER>PAM' motifs, written from 5' to 3'. '>' points toward PAM. \
    #    IUPAC ambiguities accepted. '{a,b}' are quantifiers. \
    #    Examples: 'G{,2}N{19,20}>NGG', 'N{17,20}>NGA', 'N{20,21}>NNGRRT', 'TTTN<N{20,23}'")
    
    gt_count = motif.count('>')
    lt_count = motif.count('<')
    lb_count = motif.count('{')
    rb_count = motif.count('}')
    # Make sure motif does not violate basic rules
    if (gt_count + lt_count < 1):
        raise Exception("Motif lacks distinction between spacer and PAM sequences ('>' or '<' character)")
    elif (gt_count + lt_count > 1):
        raise Exception("Motif has too many '>' or '<' characters")
    if (lb_count != rb_count):
        raise Exception("Motif braces '{' and '}' do not match")
    if (motif.count(' ') > 0):
        raise Exception("Motif contains invalid space ' ' characters")
    if (motif.count('{}') > 0):
        raise Exception("Motif contains invalid quantifier '{}'")
    if (motif.count('{,}') > 0):
        raise Exception("Motif contains invalid quantifier '{,}'")
    
    if (gt_count == 1):
        spacer_motif, pam_motif = motif.split('>')
        side = '>'
    elif (lt_count == 1):
        pam_motif, spacer_motif = motif.split('<')
        side = '<'
    
    return _parse_motif_helper(spacer_motif), _parse_motif_helper(pam_motif), side

def list_pam_sites():
    # Code taken from
    #  https://github.com/maximilianh/crisporWebsite/crispor.py
    pamDesc = [
        # Motif      sgRNA+motif origin system
        ('NGG',      '20bp-NGG - Cas9 Streptococcus Pyogenes and Cas9-HF1'),
        ('TTTN',     'TTTN-23bp - Cpf1 Acidaminococcus / Lachnospiraceae'),
        #('TTN',      'TTN-23bp - Cpf1 F. Novicida'), # Jean-Paul: various people have shown that it's not usable yet
        ('NGA',      '20bp-NGA - Cas9 S. Pyogenes mutant VQR'),
        #('NAG',       '???'),
        ('NGCG',     '20bp-NGCG - Cas9 S. Pyogenes mutant VRER'),
        ('NNAGAA',   '20bp-NNAGAA - Cas9 S. Thermophilus'),
        ('NGGNG',    '20bp-NGGNG - Cas9 S. Thermophilus'),
        ('NNGRRT',   '21bp-NNG(A/G)(A/G)T - Cas9 S. Aureus'),
        ('NNNNGMTT', '20bp-NNNNG(A/C)TT - Cas9 N. Meningitidis'),
        ('NNNNACA',  '20bp-NNNNACA - Cas9 Campylobacter jejuni'),
    ]
    
    # from (https://benchling.com/pub/cpf1):
    # Cpf1 is an RNA-guided nuclease, similar to Cas9. It recognizes a T-rich
    # PAM, TTTN, but on the 5' side of the guide. This makes it distinct from
    # Cas9, which uses an NGG PAM on the 3' side. The cut Cpf1 makes is
    # staggered. In AsCpf1 and LbCpf1, it occurs 19 bp after the PAM on the
    # targeted (+) strand and 23 bp on the other strand, as shown here:
    #                        Cpf1
    #   TTTC GAGAAGTCATCTAATAAGG|CCAC TGTTA
    #   AAAG CTCTTCAGTAGATTATTCC GGTG|ACAAT
    #   -PAM =========gRNA====== =
    #
    # Benchling suggests 20 nt guides can be used for now, as there is not real
    # suggestion for optimal guide length. Robust guide scores for Cpf1 are
    # still in development, but simple scoring based on the number of off-target
    # sites is available on Benchling.
    # 
    # Cpf1 requires only a crRNA for activity and does not need a tracrRNA to
    # also be present.
    #
    # Two Cp1-family proteins, AsCpf1 (from Acidaminococcus)
    # and LbCpf1 (from Lachnospiraceae), have been shown to perform efficient
    # genome editing in human cells.
    #
    # Why use Cpf1 over Cas9?
    #  see https://benchling.com/pub/cpf1

#def score(sequence, algorithms):
#    """Code that scores a gRNA sequence
#    Returns scores"""
#    # two types of off-target scores
#    #  CFD off-target score
#    #  MIT off-target score
#    
#    # Histogram of off-targets:
#    #  For each number of mismatches, the number of off-targets is indicated.
#    #  Example:
#    #   1-3-20-50-60    This means 1 off-target with 0 mismatches, 3 off-targets with 1 mismatch, 20 off-targets with 2 mismatches, etc.
#    #   0-2-5-10-20     These are the off-targets that have no mismatches in the 12 bp adjacent to the PAM. These are the most likely off-targets.
#    #   
#    #   Off-targets are considered if they are flanked by one of the motifs NGG, NAG or NGA.

# Template 1
def generate_excise_target(args, feature):
    """Finds the gRNA sequence to cut within the specified places
    on the feature. May target either the + or - strand.
    """
    pass

def generate_revert_target(args, feature):
    """Creates the gRNA sequence that targets the excise donor DNA oligo.
    May target either the + or - strand.
    """
    pass

def generate_excise_donor(args, feature, contigs, revert_target):
    """
    Creates the DNA oligo with the structure:
    [upstream homology][unique gRNA][downstream homology]
    that excises the target feature
    """
    
    # Generate the full set of potential dDNAs
    
    # First, get the homology blocks up- and down-stream of the feature
    contig, start, end, strand = feature
    # start & end are 0-based indices, inclusive
    
    # assumes start < end
    # DNA 5' of feature is upstream
    upstream = contigs[contig][start-args.excise_donor_homology[1]:start]
    
    # DNA 3' of feature is downstream
    downstream = contigs[contig][end+1:end+1+args.excise_donor_homology[1]]
    
    # mini_addtags = itertools.something()
    # dDNAs = [upstream + x + downstream for x in itertools.something()]
    
    dDNAs = []
    targets = []
    
    # For each potential dDNA, evaluate how good it is
    for insert_length in range(args.excise_insert_lengths[0], args.excise_insert_lengths[1]+1):
        if insert_length >= 0:
            for mAT in nucleotides.kmers(insert_length):
                # Add this candidate dDNA to the list of all candidate dDNAs
                dDNAs.append(upstream + mAT + downstream)
                
                # Use code to generate targets for this revised region
                # between [maximum 5' distance] upstream mAT downstream [maximum 3' distance]
                targets.extend(get_targets(...))
        else:
            pass
    
    # Write the targets to a FASTA file
    query_file = utils.generate_query(os.path.join(args.folder, 'reversion-query.fasta'), targets)
    
    # Use selected alignment program to find all potential off-targets in the genome
    sam_file = align(query_file, index_file, args.folder, args.processors)
    
    # Calculate scores for each target
    

def generate_revert_donor(args, feature):
    """Use template DNA sequence to create oligo that will be used for fixing
    the DSB
    """
    # Optionally provide a list of restriction enzyme targeting sites
    # on either side for easier cloning
    pass

def format_output():
    # output should take the following format
    # contig    gRNA-start-pos   gRNA-end-pos    gRNA-strand   features    gRNA-sequence    PAM     on-target-score    off-target-score
    pass

def merge_features(features):
    """Combine overlapping features?
    """
    return features

def process(args):
    # The Cas9 cuts 3-4bp upstream of the PAM sequence:
    #  ======= gRNA ==== === PAM
    #  CGATGCATCGACTTTAC CGA AGG
    #                   ^ Cut
    
    features = []
    for f in features:
        # identify_pam_positions()
        et = generate_excise_target(args, f)
        rt = generate_revert_target(args, f)
        
        ed = generate_excise_donor(args, f, rt)
        rd = generate_revert_donor(args, f)
    
    # Cut site can be anywhere within target feature
    # genome          ACGGATTAGAGAGAGGCCTCCTCCCGGAT-GAATGGAAGACTAAACGGTAGATATAAG
    # feature         -------|ORF->                                <-ORF|-------
    # PAM                                              PAM         PAM
    # excise target                         CCCGGAT^GAATGG
    # excise genome   ACGGATTAGAGAGAGGCCTCCTCCCGGAT^GAATGGAAGACTAAACGGTAGATATAAG
    # excise donor    ACGGATT-----------------------------------AAACGGTAGATATAAG
    # revert target      GATT-----------------------------------AAACGG
    # revert donor     CGGATTAGAGAGAGGCCTCCTCCCGGAT-GAATGGAAGACTAAACGGTAGATATA
    
    pass

# Template 2
def find_candidate_targets(args):
    pass

def align_targets():
    pass

def find_targets():
    pass

def make_primers():
    pass

def make_donors():
    pass

def excise():
    candidate_targets = find_targets()
    
    excise_primers = make_primers()
    
    excise_donors, revert_targets = make_donors()

def revert(revert_targets):
    revert_donors = make_donors()
    
    revert_primers = make_primers()

def target_filter(seq, args):
    '''
    Filters the candidate gRNA sequence based on the following criteria:
     1) case: ignore, discard-lower, discard-upper (does not process invariant-lower/invariant-upper)
     2) ambiguous characters: discard, keep, disambiguate
     3) maximum consecutive Ts
     4) PAM site
     * %GC
     
    Assumes the sequence is of the appropriate length
    '''
    # check for ambiguity expansion
    # Check for PAM motif
    # check for poly-T
    
    # Check the case of the potential gRNA sequence
    if (args.case == "discard-lower"):
        if regex.search('[a-z]', seq):
            return [] # Reject this sequence because it has lower-case characters
    elif (args.case == "discard-upper"):
        if regex.search('[A-Z]', seq):
            return [] # Reject this sequence because it has upper-case characters
    # if (args.case == "ignore"), then do nothing
    
    # Convert sequences to upper-case so it can be evaluated by the scoring algorithms
    seq = seq.upper()
    
    # Check if target sequence has any ambiguities
    if (args.ambiguities == 'discard'):
        if regex.search('[^ATCGatcg]', seq):
            # Reject this sequence
            return []
        seqs = [seq]
    # Disambiguate sequences if necessary
    elif (args.ambiguities == 'disambiguate'):
        seqs = nucleotides.disambiguate_iupac(seq)
    # Do nothing if just 'keep'
    else:
        seqs = [seq]
    
    # Remove targets with T{5,}
    seqs = [ nt for nt in seqs if ('T'*(args.max_consecutive_ts+1) not in nt) ]
    
    # Remove targets that do not end with intended PAM sites
    temp_nts = []
    temp_targets = []
    temp_pams = []
    for nt in seqs:
        m = nucleotides.split_target_sequence(nt, args.pams)
        if m:
            temp_nts.append(nt)
            temp_targets.append(m[0])
            temp_pams.append(m[1])
    seqs = temp_nts
    #nts = [ nt for nt in nts if re_compiled.search(nt) ]
    
    # Remove targets whose %GC is outside the chosen bounds
    # hard-coded PAM length at 3
    temp_nts = []
    for i in range(len(nts)):
        if (args.target_gc[0] <= scores.gc_score(temp_targets) <= args.target_gc[1]):
            temp_nts.append(nts[i])
    seqs = temp_nts
    #nts = [ nt for nt in nts if (target_gc[0] <= scores.gc_score(nt[:-3]) <= target_gc[1]) ]
    
    # Add all targets for this feature to the targets list
    return seqs

def get_targets(args, contigs, features):
    ## Build a regex to only match strings with PAM sites specified in args.pams
    ## Assumes PAM site is at 3' end
    #re_pattern = '|'.join([ nucleotides.build_regex_pattern(p)+'$' for p in args.pams ])
    ##re_pams = []
    ##for p in pams:
    ##    re_pams.append(nucleotides.build_regex_pattern(p) + '$')
    ##re_pattern = '|'.join(re_pams)
    re_flags = regex.ENHANCEMATCH | regex.IGNORECASE
    #re_compiled = regex.compile(re_pattern, flags=re_flags)
    
    # Ideally, this code would procedurally write to the query file
    # as new sequences were being added, thus limiting the amount of memory
    # used to store sequences.
    
    # Find unique gRNA sites within each feature
    targets = []
    # Use a sliding window to make a list of queries
    for feature in features:
        
        contig, start, end, strand = features[feature]
        if (end == None):
            end = len(contigs[contig])
        
        # Make sure the contig the feature is on is present in the FASTA
        if contig in contigs:
            # print(feature, features[feature], file=sys.stderr)
            # Find a site within this feature that will serve as a unique gRNA
            
            # for each orientation:
            # if (args.strands in ['+', 'both']):
            #  targets.extend...
            # if (args.strands in ['-', 'both']):
            #  targets.extend...
            
            # Search both the '+' and '-' strands
            for contig_seq in [contigs[contig], utils.rc(contigs[contig])]:
                for pam in args.pams:
                    #pam_seq = 'NGG'
                    #pam_regex = '([ACGT]GG)'
                    pam_length = len(pam)
                    #pam_regex_compiled = regex.compile(pam_regex+'$', flags=re_flags
                    
                    # Problem: what if user wants lengths 17 & 20, but not 18 & 19?
                    # solution...use that regex-stype motif command line argument
                    for target_length in range(target_lengths[0], target_lengths[1]+1):
                        for pos in range(start, end+1-target_length+pam_length, 1):
                            nt = contig_seq[pos:pos+target_length]
                            
                            # Skip this sequence if the mask is on
                            # and it has lower-case characters
                            if (args.lower_case_mask):
                                if regex.search('[a-z]', nt):
                                    continue
                            
                            # Convert sequences to upper-case
                            nt = nt.upper()
                            
                            # Disambiguate sequences if necessary
                            if (args.ambiguities == 'discard'):
                                if regex.search('[^ATCGatcg]', nt):
                                    continue
                                nts = [nt]
                            elif (args.ambiguities == 'disambiguate'):
                                nts = nucleotides.disambiguate_iupac(nt)
                            
                            # Remove targets with T{5,}
                            nts = [ nt for nt in nts if ('T'*(args.max_consecutive_ts+1) not in nt) ]
                            
                            # Remove targets that do not end with intended PAM sites
                            temp_nts = []
                            temp_targets = []
                            temp_pams = []
                            for nt in nts:
                                m = nucleotides.split_target_sequence(nt, args.pams)
                                if m:
                                    temp_nts.append(nt)
                                    temp_targets.append(m[0])
                                    temp_pams.append(m[1])
                            nts = temp_nts
                            #nts = [ nt for nt in nts if re_compiled.search(nt) ]
                            
                            # Remove targets whose %GC is outside the chosen bounds
                            # hard-coded PAM length at 3
                            temp_nts = []
                            for i in range(len(nts)):
                                if (args.target_gc[0] <= scores.gc_score(temp_targets) <= args.target_gc[1]):
                                    temp_nts.append(nts[i])
                            nts = temp_nts
                            #nts = [ nt for nt in nts if (target_gc[0] <= scores.gc_score(nt[:-3]) <= target_gc[1]) ]
                            
                            # Add all targets for this feature to the targets list
                            #targets.extend(map(lambda x: (feature, contig,)+x, utils.sliding_window(contigs[contig], window=target_length, start=start, stop=end)))
                            #targets.extend(map(lambda x: (feature, contig, pos, pos+target_length, x), nts))
                            for nt in nts:
                                targets.append((feature, contig, pos, pos+target_length, nt))
    
    return targets
    
    
            # # Get all potential gRNAs from feature
            # ts = nucleotides.SlidingWindow(contigs[contig], window=target_length+3, start=start, stop=end)
            # 
            # # Convert sequences to upper-case
            # ts = [ (s, e, seq.upper()) for s, e, seq in ts ]
            # 
            # # Disambiguate sequences if necessary
            # if (ambiguities == 'discard'):
            #     ts = [ item for item in ts if not regex.search('[^ATCGatcg]', item[2]) ]
            # 
            # # This code takes way too much memory
            # elif (ambiguities == 'disambiguate'):
            #     #tts = []
            #     #for s in ts:
            #     #    for ds in utils.disambiguate_iupac(s[2]):
            #     #        tts.append((s[0], s[1], ds))
            #     #ts = tts
            #     #ts = [ map(lambda x: (s, e, x), utils.disambiguate_iupac(seq)) for s, e, seq in ts ]
            #     #ts = utils.flatten(ts)
            #     ts = [(a[0], a[1], x) for a in ts for x in nucleotides.disambiguate_iupac(a[2])]
            # 
            # # Remove targets with T{5,}
            # # ts = utils.filter_polyt(ts, args.max_consecutive_ts)
            # ts = [ item for item in ts if ('T'*(max_consecutive_ts+1) not in item[2]) ]
            # 
            # # Remove targets that do not end with intended PAM sites
            # ts = [ item for item in ts if re_compiled.search(item[2]) ]
            # 
            # # Remove targets whose %GC is outside the chosen bounds
            # # hard-coded PAM length at 3
            # ts = [ item for item in ts if (target_gc[0] <= scores.gc_score(item[2][:-3]) <= target_gc[1]) ]
            # 
            # # Add all targets for this feature to the targets list
            # #targets.extend(map(lambda x: (feature, contig,)+x, utils.sliding_window(contigs[contig], window=target_length, start=start, stop=end)))
            # targets.extend(map(lambda x: (feature, contig,)+x, ts))
            # 
            # #for target in utils.sliding_window(contigs[contig][start:end], target_length):
            # #    Use regex (slow) to find all matches in the genome
            # #    regex = utils.build_regex(target, max_substitutions=2)
            # #    matches = utils.find_target_matches(regex, contigs, overlap=True)
            # #    for m in matches:
            # #        print(m)
    

def index_reference(fasta, folder, processors):
    if (args.aligner == 'addtag'):
        index_file = fasta
    elif (args.aligner == 'blast+'):
        pass
    elif (args.aligner == 'blat'):
        pass
    elif (args.aligner == 'bowtie'):
        pass
    elif (args.aligner == 'bowtie2'):
        index_file = bowtie2.index_reference(fasta, tempdir=folder, threads=processors)
    elif (args.aligner == 'bwa'):
        pass
    elif (args.aligner == 'cas-offinder'):
        pass
    return index_file

def align(query_file, index_file, folder, processors):
    if (args.aligner == 'addtag'):
        sam_file = None
    elif (args.aligner == 'blast+'):
        pass
    elif (args.aligner == 'blat'):
        pass
    elif (args.aligner == 'bowtie'):
        pass
    elif (args.aligner == 'bowtie2'):
        sam_file = bowtie2.align(query_file, index_file, folder=folder, threads=processors)
    elif (args.aligner == 'bwa'):
        pass
    elif (args.aligner == 'cas-offinder'):
        pass
    return sam_file

def main():
    """Function to run complete AddTag analysis"""
    
    # Obtain command line arguments and parse them
    args = parse_arguments()
    
    # outputs:
    #  folder/                            folder holding output
    #  folder/log.txt                     log
    #  folder/bowtie2-index/*             bowtie2 index for input FASTA
    #  folder/excision-query.fasta             FASTA file of all candidate gRNAs
    #  folder/excision-query.sam               SAM file for alignment of candidate gRNAs
    #  folder/excision-query.err               STDOUT/STDERR from alignment
    #  folder/excision-gRNAs.fasta        sequences to be synthesized/amplified/transcribed (header contains scores)
    #  folder/excision-dDNAs.fasta        sequences to be synthesized (paired with reversion-gRNAs)
    #  folder/reversion-query.fasta
    #  folder/reversion-query.sam
    #  folder/reversion-query.err
    #  folder/reversion-gRNAs.fasta       sequences to be synthesized/amplified/transcribed (paired with excision-dDNAs)
    #  folder/reversion-dDNAs.fasta       sequence to be amplified
    #  folder/reversion-primers.fasta     primers for amplifying reversion dDNAs
    #  folder/off-target-dDNAs.fasta      wt dDNAs to prevent mutations at off-target Cas9 binding sites (contains edit distance in header s/i/d)
    #  folder/primers.fasta               primers to check for KO/KI (contains expected amplicon sizes in header)
    
    if args.test:
        # Perform test code
        test(args)
    else:
        # Create the project directory if it doesn't exist
        os.makedirs(args.folder, exist_ok=True)
        
        # Load the FASTA file specified on the command line
        contigs = utils.load_fasta_file(args.fasta)
        
        # Open and parse the GFF file specified on the command line
        features = utils.load_gff_file(args.gff, args.features, args.tag)
        
        # Merge features?
        #features = merge_features(features)
        
        # Generate the query list: [(feature, contig, start, end, sequence), ...]
        targets = get_targets(args, contigs, features)
        
        # Index the reference FASTA
        index_file = index_reference(args.fasta, args.folder, args.processors)
        
        # Write the query list to FASTA
        query_file = utils.generate_query(os.path.join(args.folder, 'excision-query.fasta'), targets)
        
        # Use selected alignment program to find all matches in the genome
        sam_file = align(query_file, index_file, args.folder, args.processors)
        
        # Open the SAM file
        alignments = load_sam_file_test(os.path.join(args.folder, 'excision-query.sam'), args.pams, contigs)
        for s in alignments:
            print(s)
            for a in s.alignments:
                print('  ', a)
        
        # Discard potential gRNAs that have mismatches with their target site
        #if (args.case == "invariant-lower"):
        #    pass
        #elif (args.case == "invariant-upper"):
        #    pass

def test(args):
    """Code to test the classes and functions in 'source/__init__.py'"""
    # Echo the command line parameters
    print(args, file=sys.stderr)
    
    sys.exit()
    
    # Get timestamp
    start = time.time()
    
    # Load the FASTA file specified on the command line
    print("=== FASTA ===")
    contigs = utils.load_fasta_file(args.fasta)
    print(list(contigs.keys())[:5])
    
    # Test SAM file parsing
    print("=== SAM ===")
    try:
        alignments = load_sam_file_test(os.path.join(args.folder, 'excision-query.sam'), args.pams, contigs)
        for s in alignments:
            print(s)
            for a in s.alignments:
                print('  ', a)
    except FileNotFoundError:
        print('Skipping...')
    
    # Open and parse the GFF file specified on the command line
    # returns a dictionary:
    #  features[ID] = (contig, start(bp), end(bp), frame)
    print("=== GFF ===")
    features = utils.load_gff_file(args.gff, args.features, args.tag)
    for k in list(features.keys())[:10]:
        print(k, features[k])
    
    # Test code to find all similar oligonucleotides in the FASTA
    print("=== Align ===")
    target = 'TCCGGTACAKTGAKTTGTAC'
    regex = nucleotides.build_regex(target, max_errors=2)
    matches = nucleotides.find_target_matches(regex, contigs, overlap=True)
    for m in matches:
        print(m)
        for seq in nucleotides.disambiguate_iupac(m[4]):
            print(seq, len(seq), hsuzhang.hsuzhang_score(target, seq))
    
    # Test Hsu score
    print("=== Hsu 2013 ===")
    a = 'CGATGGCTWGGATCGATTGAC'
    b = 'AAGTGCTCTTAAGAGAAATTC'
    c = 'ATGSCTCGGATCGATTGAC'
    print(hsuzhang.calcHitScore(a, a), hsuzhang.hsuzhang_score(a, a), hsuzhang.hsuzhang_score(a, a, True))
    print(hsuzhang.calcHitScore(a, b), hsuzhang.hsuzhang_score(a, b), hsuzhang.hsuzhang_score(a, b, True))
    print(hsuzhang.calcHitScore(a, c), hsuzhang.hsuzhang_score(a, c), hsuzhang.hsuzhang_score(a, c, True))
    
    # Test Doench 2014 score:
    print("=== Doench 2014 ===")
    pam = 'AGG'
    gRNAa = 'CGATGGCTTGGATCGATTGA'
    gRNAb = 'CGTTGGCTTGGATCGATTGA'
    gRNAc = 'CGATGGCTTCGATCGATTGA'
    gRNAd = 'CGATGGCTTCGAGCGATTGA'
    print(doench.on_target_score_2014(gRNAa, pam))
    print(doench.on_target_score_2014(gRNAb, pam))
    print(doench.on_target_score_2014(gRNAc, pam))
    print(doench.on_target_score_2014(gRNAd, pam))
    
    # Test Doench 2016 score
    print("=== Doench 2016 ===")
    pam = 'AGG'
    gRNAa = 'CGATGGCTTGGATCGATTGA'
    gRNAb = 'CGTTGGCTTGGATCGATTGA'
    gRNAc = 'CGATGGCTTCGATCGATTGA'
    gRNAd = 'CGATGGCTTCGAGCGATTGA'
    print(doench.on_target_score_2016(gRNAa, gRNAa, pam))
    print(doench.on_target_score_2016(gRNAa, gRNAb, pam))
    print(doench.on_target_score_2016(gRNAa, gRNAc, pam))
    print(doench.on_target_score_2016(gRNAa, gRNAd, pam))
    
    # Print time taken for program to complete
    print('Runtime: {}s'.format(time.time()-start), file=sys.stderr)