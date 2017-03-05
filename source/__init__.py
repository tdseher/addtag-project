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
    def __init__(self, feature, contig_sequence, pams, contig=None, contig_start=None, contig_end=None, contig_orientation=None, feature_orientation=None):
        """Create a structure for holding individual sequence information"""
        self.feature = feature
        self.feature_orientation = feature_orientation
        
        self.contig_sequence = contig_sequence
        self.contig_target, self.contig_pam = nucleotides.split_target_sequence(self.contig_sequence, pams)
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
        
        # query x subject score
        self.off_target_doench2014 = None
        self.off_target_doench2016 = None
        self.off_target_hsuzhang = None
        self.off_target_linear = None
    
    def add_alignment(self, aligned_sequence, pams, aligned_contig, aligned_start, aligned_end, aligned_orientation):
        """Add a genomic position to the list of alignments"""
        aligned_target, aligned_pam = nucleotides.split_target_sequence(aligned_sequence, pams)
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
        )
        self.alignments.append(seq)
    
    def __repr__(self):
        return 'Sequence(feature=' + self.feature + ', ' + \
            self.contig + ':' + str(self.contig_start) + '..' + \
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
    parser.add_argument("--pams", metavar="SEQ", nargs="+", type=str,
        default=["NGG"], help="Constrain finding only targets with these PAM sites")
    parser.add_argument("--target_lengths", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[17, 20],
        help="The length range of the 'target'/'spacer'/gRNA site")
    # Eventually, replace --pams and --target_lengths with this:
    #parser.add_argument("--motifs", metavar="SEQ", nargs="+", type=str,
    #    default=["N{17,20}>NGG"], help="Find only targets with these \
    #    'SPACER>PAM' motifs, written from 5' to 3'. '>' points toward PAM. \
    #    IUPAC ambiguities accepted. '{a,b}' are quantifiers. \
    #    Examples: 'G{,2}N{19,20}>NGG', 'N{17,20}>NGA', 'N{20,21}>NNGRRT', 'TTTN<N{20,23}'")
    parser.add_argument("--tag", metavar='TAG', type=str, default='ID',
        help="GFF tag with feature names. Examples: 'ID', 'Name', 'Parent', or 'locus_tag'")
    #parser.add_argument("--feature_homolog_regex", metavar="REGEX", type=str, default=None, help="regular expression with capturing group containing invariant feature. Example: '(.*)_[AB]' will treat features C2_10010C_A and C2_10010C_B as homologs")
    # okay idea, but needs more thought before implementation
    parser.add_argument("--feature_homologs", metavar="*.homologs", type=str, default=None,
        help="Path to text file containing homologous features on the same \
            line, separated by TAB characters")
    parser.add_argument("--min_contig_edge_distance", metavar="N", type=int, default=500,
        help="Minimum distance from contig edge a site can be found")
    parser.add_argument("--features", metavar="FEATURE", type=str, nargs="+", default=["gene"],
        help="Features to design gRNA sites against. Must exist in GFF file. Examples: 'CDS', 'gene', 'mRNA', 'exon'")
    parser.add_argument("--target_gc", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[25, 75],
        help="Generated gRNAs must have %%GC content between these values (excluding PAM motif)")
    parser.add_argument("--excise_donor_lengths", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[90, 100],
        help="Range of lengths acceptable for knock-out dDNAs.")
    parser.add_argument("--revert_donor_lengths", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[300, 600],
        help="Range of lengths acceptable for knock-in dDNAs.")
    parser.add_argument("--min_feature_edge_distance", metavar="N", type=int, default=23,
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
    parser.add_argument("--strands", type=str, choices=["+", "-", "both"], default="both",
        help="Strands to search for gRNAs")
    parser.add_argument("--ambiguities", type=str, choices=["discard", "disambiguate", "keep"], default="discard",
        help="How generated gRNAs should treat ambiguous bases: \
        discard - no gRNAs will be created where the FASTA has an ambiguous base; \
        disambiguate - gRNAs containing ambiguous bases will be converted to a set of non-ambiguous gRNAs; \
        keep - gRNAs can have ambiguous bases")
    parser.add_argument("--lower_case_mask", action="store_true",
        help="Do not generated gRNAs targeting lower case nucleotides in input FASTA")
    parser.add_argument("--overlap", action="store_true",
        help="Include exhaustive search for overlapping sites. May increase computation time.")
    parser.add_argument("--processors", metavar="N", type=int, default=(os.cpu_count() or 1),
        help="Number of processors to use when performing pairwise sequence alignments")
    parser.add_argument("--aligner", type=str, choices=['bowtie', 'bowtie2', 'bwa'], default='bowtie2',
        help="Program to calculate pairwise alignments")
    # Other aligners to consider: 'bowtie', 'bwa', 'blastn', 'blat', 'rmap', 'maq', 'shrimp2', 'soap2', 'star', 'rhat', 'mrsfast', 'stampy'
    
    # Special version action
    parser.add_argument("-v", "--version", action='version', version='{__program__} {__version__}'.format(**globals()))
    
    # Parse the arguments, and return
    return parser.parse_args()

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

def generate_excise_donor(args, feature, revert_target):
    """Creates the DNA oligo with the structure:
    homology--unique gRNA--homology
    that excises the target feature
    """
    pass

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
def find_candidate_targets():
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

def main():
    """Function to run complete AddTag analysis"""
    
    # Obtain command line arguments and parse them
    args = parse_arguments()
    
    # outputs:
    #  folder/                            folder holding output
    #  folder/log.txt                     log
    #  folder/bowtie2-index/*             bowtie2 index for input FASTA
    #  folder/alignment.fasta             FASTA file of all candidate gRNAs
    #  folder/alignment.sam               SAM file for alignment of candidate gRNAs
    #  folder/alignment.err               STDOUT/STDERR from alignment
    #  folder/excission-gRNAs.fasta       sequences to be synthesized/amplified/transcribed (header contains scores)
    #  folder/excission-dDNAs.fasta       sequences to be synthesized (paired with reversion-gRNAs)
    #  folder/reversion-gRNAs.fasta       sequences to be synthesized/amplified/transcribed (paired with escission-dDNAs)
    #  folder/reversion-dDNAs.fasta       sequence to be amplified
    #  folder/reversion-primers.fasta     primers for amplifying reversion dDNAs
    #  folder/off-target-dDNAs.fasta      wt dDNAs to prevent mutations at off-target Cas9 binding sites (contains edit distance in header s/i/d)
    #  folder/primers.fasta               primers to check for KO/KI (contains expected amplicon sizes in header)
    
    if args.test:
        # Perform test code
        test(args)
    else:
        # Load the FASTA file specified on the command line
        contigs = utils.load_fasta_file(args.fasta)
        
        # Open and parse the GFF file specified on the command line
        features = utils.load_gff_file(args.gff, args.features, args.tag)
        
        # Merge features?
        #features = merge_features(features)
        
        # Set gRNA target_length to be 20 nt, despite what the command line specifies
        # for simple testing purposes
        target_length = 20
        
        # Create the project directory if it doesn't exist
        os.makedirs(args.folder, exist_ok=True)
        
        # Index the reference FASTA
        if (args.aligner == 'bowtie2'):
            index_file = bowtie2.index_reference(args.fasta, tempdir=args.folder, threads=args.processors)
        
        # Build a regex to only match strings with PAM sites specified in args.pams
        re_pams = []
        for p in args.pams:
            re_pams.append(nucleotides.build_regex_pattern(p) + '$')
        re_pattern = '|'.join(re_pams)
        re_flags = regex.ENHANCEMATCH | regex.IGNORECASE
        re_compiled = regex.compile(re_pattern, flags=re_flags)
        
        # Ideally, this code would procedurally write to the query file
        # as new sequences were being added, thus limiting the amount of memory
        # used to store sequences.
        
        # Find unique gRNA sites within each feature
        targets = []
        # Use a sliding window to make a list of queries
        for feature in features:
            
            # for each orientation:
            # if (args.strands in ['+', 'both']):
            #  targets.extend...
            # if (args.strands in ['-', 'both']):
            #  targets.extend...
            
            contig, start, end, strand = features[feature]
            # Make sure the contig the feature is on is present in the FASTA
            if contig in contigs:
                # print(feature, features[feature], file=sys.stderr)
                # Find a site within this feature that will serve as a unique gRNA
                
                # Get all potential gRNAs from feature
                ts = nucleotides.SlidingWindow(contigs[contig], window=target_length+3, start=start, stop=end)
                
                # Convert sequences to upper-case
                ts = [ (s, e, seq.upper()) for s, e, seq in ts ]
                
                # Disambiguate sequences if necessary
                if (args.ambiguities == 'discard'):
                    #ts = [ item for item in ts if item[2] ]
                    pass
                # This code takes way too much memory
                if (args.ambiguities == 'disambiguate'):
                    #tts = []
                    #for s in ts:
                    #    for ds in utils.disambiguate_iupac(s[2]):
                    #        tts.append((s[0], s[1], ds))
                    #ts = tts
                    #ts = [ map(lambda x: (s, e, x), utils.disambiguate_iupac(seq)) for s, e, seq in ts ]
                    #ts = utils.flatten(ts)
                    ts = [(a[0], a[1], x) for a in ts for x in nucleotides.disambiguate_iupac(a[2])]
                
                # Remove targets with T{5,}
                # ts = utils.filter_polyt(ts, args.max_consecutive_ts)
                ts = [ item for item in ts if ('T'*(args.max_consecutive_ts+1) not in item[2]) ]
                
                # Remove targets that do not end with intended PAM sites
                ts = [ item for item in ts if re_compiled.search(item[2]) ]
                
                # Remove targets whose %GC is outside the chosen bounds
                # hard-coded PAM length at 3
                ts = [ item for item in ts if (args.target_gc[0] <= scores.gc_score(item[2][:-3]) <= args.target_gc[1]) ]
                
                # Add all targets for this feature to the targets list
                #targets.extend(map(lambda x: (feature, contig,)+x, utils.sliding_window(contigs[contig], window=target_length, start=start, stop=end)))
                targets.extend(map(lambda x: (feature, contig,)+x, ts))
                
                #for target in utils.sliding_window(contigs[contig][start:end], target_length):
                #    Use regex (slow) to find all matches in the genome
                #    regex = utils.build_regex(target, max_substitutions=2)
                #    matches = utils.find_target_matches(regex, contigs, overlap=True)
                #    for m in matches:
                #        print(m)
        
        name = 'alignment'
        query_file = bowtie2.generate_query(os.path.join(args.folder, name+'.fasta'), targets)
        
        # Use bowtie2 to find all matches in the genome
        sam_file = bowtie2.align(query_file, index_file, folder=args.folder)
        
def test(args):
    """Code to test the classes and functions in 'source/__init__.py'"""
    # Echo the command line parameters
    print(args, file=sys.stderr)
    
    # Get timestamp
    start = time.time()
    
    # Load the FASTA file specified on the command line
    print("=== FASTA ===")
    contigs = utils.load_fasta_file(args.fasta)
    print(list(contigs.keys())[:5])
    
    # Test SAM file parsing
    print("=== SAM ===")
    try:
        alignments = load_sam_file_test(os.path.join(args.folder, 'alignment.sam'), args.pams, contigs)
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