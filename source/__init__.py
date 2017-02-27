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
from . import hsu
from . import doench
from . import bowtie2

# Define meta variables
__author__ = "Thaddeus D. Seher (@tdseher) & Aaron Hernday"
__date__ = "2017-02-03"
__fullversion__ = utils.load_git_version()
__version__ = __fullversion__[:7]
__program__ = os.path.basename(sys.argv[0])
__description__ = """\
description:
  Program for identifying unique endogenous gRNA sites 
  and creating unique synthetic gRNA sites.
  
  Copyright (c) {__date__} {__author__}.
  All rights reserved.

version:
  short {__version__}
  full  {__fullversion__}

""".format(**locals())
__epilog__ = """\
example:
 $ python3 {__program__} genome.fasta genome.gff > results.txt
""".format(**locals())

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
    
    # Add optional arguments
    parser.add_argument("--test", action="store_true",
        help="Perform tests only")
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
    parser.add_argument("--target_lengths", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[19, 21],
        help="The length range of the 'target'/'spacer'/gRNA site")
    parser.add_argument("--donor_lengths", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[80, 100],
        help="The length range of the final computed donor DNA for each site")
    parser.add_argument("--min_feature_edge_distance", metavar="N", type=int, default=24,
        help="The minimum distance a gRNA site can be from the edge of the \
             feature. If negative, the maximum distance a gRNA site can be \
             outside the feature.")
    parser.add_argument("--min_donor_insertions", metavar="N", type=int, default=2,
        help="The uniqueness of final donor DNA compared to the rest of the genome")
    parser.add_argument("--min_donor_deletions", metavar="N", type=int, default=2,
        help="The uniqueness of final donor DNA compared to the rest of the genome")
    parser.add_argument("--min_donor_substitutions", metavar="N", type=int, default=2,
        help="The uniqueness of final donor DNA compared to the rest of the genome")
    parser.add_argument("--min_donor_errors", metavar="N", type=int, default=3,
        help="The uniqueness of final donor DNA compared to the rest of the genome")
    parser.add_argument("--min_donor_distance", metavar="N", type=int, default=36,
        help="The minimum distance in bp a difference can exist from the edge of donor DNA") # homology with genome
    parser.add_argument("--strands", type=str, choices=["+", "-", "both"], default="both",
        help="Strands to search for gRNAs")
    parser.add_argument("--overlap", action="store_true",
        help="Include exhaustive search for overlapping sites. May increase computation time.")
    parser.add_argument("--processors", type=int, default=(os.cpu_count() or 1),
        help="Number of processors to use when performing pairwise sequence alignments")
    parser.add_argument("--aligner", type=str, choices=['bowtie2'], default='bowtie2',
        help="Program to calculate pairwise alignments")
    # Other aligners to consider: 'bowtie', 'bwa', 'blastn', 'blat', 'rmap', 'maq', 'shrimp2', 'soap2', 'star', 'rhat', 'mrsfast', 'stampy'
    parser.add_argument("--temp", type=str, default=os.getcwd(),
        help="Path of folder to store temporary files")
    
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
        ('NGCG',     '20bp-NGCG - Cas9 S. Pyogenes mutant VRER'),
        ('NNAGAA',   '20bp-NNAGAA - Cas9 S. Thermophilus'),
        ('NGGNG',    '20bp-NGGNG - Cas9 S. Thermophilus'),
        ('NNGRRT',   '21bp-NNG(A/G)(A/G)T - Cas9 S. Aureus'),
        ('NNNNGMTT', '20bp-NNNNG(A/C)TT - Cas9 N. Meningitidis'),
        ('NNNNACA',  '20bp-NNNNACA - Cas9 Campylobacter jejuni'),
    ]

def scores():
    pass
    
    # two types of off-target scores
    #  CFD off-target score
    #  MIT off-target score
    
    # Histogram of off-targets:
    #  For each number of mismatches, the number of off-targets is indicated.
    #  Example:
    #   1-3-20-50-60    This means 1 off-target with 0 mismatches, 3 off-targets with 1 mismatch, 20 off-targets with 2 mismatches, etc.
    #   0-2-5-10-20     These are the off-targets that have no mismatches in the 12 bp adjacent to the PAM. These are the most likely off-targets.
    #   
    #   Off-targets are considered if they are flanked by one of the motifs NGG, NAG or NGA.
    

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
    # The Cas9 cuts 3-4bp upstream of the PAM sequence.
    
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

def test(args):
    """Run test code"""
    # Echo the command line parameters
    print(args, file=sys.stderr)
    
    # Get timestamp
    start = time.time()
    
    # Load the FASTA file specified on the command line
    print("=== FASTA ===")
    contigs = utils.load_fasta_file(args.fasta)
    print(list(contigs.keys())[:5])
    
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
    regex = utils.build_regex(target, max_errors=2)
    matches = utils.find_target_matches(regex, contigs, overlap=True)
    for m in matches:
        print(m)
        for seq in utils.disambiguate_iupac(m[4]):
            print(seq, len(seq), hsu.calcHitScore(target, seq))
    
    # Test Hsu score
    print("=== Hsu 2013 ===")
    a = 'CGATGGCTWGGATCGATTGAC'
    b = 'AAGTGCTCTTAAGAGAAATTC'
    c = 'ATGSCTCGGATCGATTGAC'
    print(hsu.calcHitScore(a, a), hsu.hsu_score(a, a), hsu.hsu_score(a, a, True))
    print(hsu.calcHitScore(a, b), hsu.hsu_score(a, b), hsu.hsu_score(a, b, True))
    print(hsu.calcHitScore(a, c), hsu.hsu_score(a, c), hsu.hsu_score(a, c, True))
    
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

def main():
    """Function to run complete AddTag analysis"""
    
    # Obtain command line arguments and parse them
    args = parse_arguments()
    
    if args.test:
        # Perform test code
        test(args)
    else:
        # Load the FASTA file specified on the command line
        contigs = utils.load_fasta_file(args.fasta)
        
        # Open and parse the GFF file specified on the command line
        features = utils.load_gff_file(args.gff, args.features, args.tag)
        
        # Merge features?
        features = merge_features(features)
        
        # Set gRNA target_length to be 20 nt, despite what the command line specifies
        target_length = 20
        
        # Index the reference FASTA
        if (args.aligner == 'bowtie2'):
            index_file = bowtie2.index_reference(args.fasta, tempdir=args.temp, threads=args.processors)
        
        # Ideally, this code would procedurally write to the query file
        # as new sequences were being added, thus limiting the amount of memory
        # used to store sequences.
        
        # Find unique gRNA sites within each feature
        targets = []
        # Use a sliding window to make a list of queries
        for feature in features:
            contig, start, end, strand = features[feature]
            # Make sure the contig the feature is on is present in the FASTA
            if contig in contigs:
                print(feature, features[feature], file=sys.stderr)
                # Find a site within this feature that will serve as a unique gRNA
                targets.extend(map(lambda x: (feature, contig,)+x, utils.sliding_window(contigs[contig], window=target_length, start=start, stop=end)))
                
                #for target in utils.sliding_window(contigs[contig][start:end], target_length):
                #    Use regex (slow) to find all matches in the genome
                #    regex = utils.build_regex(target, max_substitutions=2)
                #    matches = utils.find_target_matches(regex, contigs, overlap=True)
                #    for m in matches:
                #        print(m)
        
        name = 'temp_alignment'
        # os.path.join(args.temp, name+'.fasta')
        query_file = bowtie2.generate_query(name+'.fasta', targets)
        
        # Use bowtie2 to find all matches in the genome
        sam_file = bowtie2.align(query_file, index_file)
        
