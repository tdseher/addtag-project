#!/usr/bin/env python3

"""Program for identifying unique endogenous gRNA sites
and creating unique synthetic gRNA sites."""

# source/main.py

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

# Define meta variables
__author__ = "Thaddeus D. Seher (@tdseher), & Aaron Hernday"
__date__ = "2017-02-03"
__version__ = "1"
__program__ = os.path.basename(sys.argv[0])
__description__ = """\
Program for identifying unique endogenous gRNA sites
and creating unique synthetic gRNA sites. \
Copyright (c) {__date__} {__author__}. All rights reserved. \
""".format(**locals())
__epilog__ = """\
example:
 $ python3 {__program__} genome.fasta genome.gff > results.txt
""".format(**locals())

def parse_arguments():
    # Create the argument parser
    parser = argparse.ArgumentParser(
        description=__description__,
        epilog=__epilog__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Add mandatory arguments
    parser.add_argument("fasta", type=str, help="FASTA file containing all contigs to derive gRNA sites for. All FASTA sequences should have unique primary headers (everything between the '>' and the first ' ' should be unique).")
    parser.add_argument("gff", type=str, help="GFF file specifying chromosomal features")
    
    # Add optional arguments
    #parser.add_argument("--homologs", metavar="FILE", type=str, default=None, help="Path to text file homologous contigs on the same line, separated by TAB characters")
    parser.add_argument("--min_contig_edge_distance", metavar="N", type=int, default=500, help="Minimum distance from contig edge a site can be found")
    #parser.add_argument("--features", metavar="FEATURE,FEATURE", type=str, default="ORF", help="Features to design gRNA sites against")
    parser.add_argument("--features", metavar="FEATURE", type=str, nargs="+", default=["ORF"], help="Features to design gRNA sites against. Must exist in GFF file.")
    parser.add_argument("--target_lengths", nargs=2, metavar=('min', 'max'), type=int, default=[19, 21], help="The length range of the 'target'/'spacer'/gRNA site")
    parser.add_argument("--donor_lengths", nargs=2, metavar=('min', 'max'), type=int, default=[80, 100], help="The length range of the final computed donor DNA for each site")
    parser.add_argument("--min_feature_edge_distance", metavar="N", type=int, default=24, help="The minimum distance a gRNA site can be from the edge of the feature. If negative, the maximum distance a gRNA site can be outside the feature.")
    parser.add_argument("--min_donor_insertions", metavar="N", type=int, default=2, help="The uniqueness of final donor DNA compared to the rest of the genome")
    parser.add_argument("--min_donor_deletions", metavar="N", type=int, default=2, help="The uniqueness of final donor DNA compared to the rest of the genome")
    parser.add_argument("--min_donor_mismatches", metavar="N", type=int, default=2, help="The uniqueness of final donor DNA compared to the rest of the genome")
    parser.add_argument("--min_donor_differences", metavar="N", type=int, default=3, help="The uniqueness of final donor DNA compared to the rest of the genome")
    parser.add_argument("--min_donor_distance", metavar="N", type=int, default=36, help="The minimum distance in bp a difference can exist from the edge of donor DNA") # homology with genome
    parser.add_argument("--strands", type=str, choices=["+", "-", "both"], default="both", help="Strands to search for gRNAs")
    parser.add_argument("--overlap", action="store_true", help="Include exhaustive search for overlapping sites. May increase computation time.")
    
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
    # Test Hsu score
    a = 'CGATGGCTWGGATCGATTGAC'
    b = 'AAGTGCTCTTAAGAGAAATTC'
    c = 'ATGSCTCGGATCGATTGAC'
    print(hsu.calcHitScore(a, a), hsu.hsu_score(a, a), hsu.hsu_score(a, a, True))
    print(hsu.calcHitScore(a, b), hsu.hsu_score(a, b), hsu.hsu_score(a, b, True))
    print(hsu.calcHitScore(a, c), hsu.hsu_score(a, c), hsu.hsu_score(a, c, True))
    
    #contigs = utils.read_fasta_file(args.fasta)
    #contigs = utils.read_fasta_file(r'C:\Users\thaddeus\Labs\Nobile lab\CRISPR-Cas9 AddTag Project\data\C_albicans_SC5314_A22_current_chromosomes.fasta')
    contigs = utils.read_fasta_file(r'C:\Users\thaddeus\Labs\Nobile lab\CRISPR-Cas9 AddTag Project\data\test.fasta')
    target = 'TCCGGTACAKTGAKTTGTAC' #AAAGTCAGAGTAGTTGTAAACRAGAAGAGAGTTTTAAACT
    regex = utils.build_regex(target, max_errors=2)
    matches = utils.find_target_matches(regex, contigs, overlap=True)
    for m in matches:
        print(m)
        for seq in utils.disambiguate_iupac(m[4]):
            print(seq, len(seq), hsu.calcHitScore(target, seq))
    
    # Test Doench score
    a = 'CGATGGCTTGGATCGATTGA' + 'AGG'
    b = 'CGTTGGCTTGGATCGATTGA' + 'AGG'
    c = 'CGATGGCTTCGATCGATTGA' + 'AGG'
    d = 'CGATGGCTTCGAGCGATTGA' + 'AGG'
    print(doench.on_target_score_2016_loader(a, a))
    print(doench.on_target_score_2016_loader(a, b))
    print(doench.on_target_score_2016_loader(a, c))
    print(doench.on_target_score_2016_loader(a, d))
    
    

def main():
    """Function to run complete AddTag analysis"""
    # Get timestamp
    start = time.time()
    
    # Obtain command line arguments and parse them
    args = parse_arguments()
    print(args, file=sys.stderr)
    
    # Convert input files to memory structures
    #contigs = utils.read_fasta_file(args)
    
    # Perform test code
    test(args)
    
    # Print time taken for program to complete
    print('Runtime: {}s'.format(time.time()-start), file=sys.stderr)
