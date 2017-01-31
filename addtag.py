#!/usr/bin/env python3

"""Program for identifying unique endogenous gRNA sites
and creating unique synthetic gRNA sites."""

# Import standard packages
import sys
import argparse
import textwrap
import time

# Import non-standard packages
import regex

# Define meta variables
__author__ = "Aaron Hernday & Thaddeus D. Seher"
__twitter__ = "@tdseher"
__date__ = "2017-01-28"
__version__ = "1"

__description__ = """\
Program for identifying unique endogenous gRNA sites
and creating unique synthetic gRNA sites. \
Copyright (c) {__date__} {__author__} ({__twitter__}). All rights reserved. \
""".format(**locals())



def parse_arguments():
    # Create the parser
    parser = argparse.ArgumentParser(
        description="""Counts how many times each read maps with certain \
                    orientations in a set of SAM files. \
                    Copyright 2015 Thaddeus Seher.""",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Add mandatory arguments
    parser.add_argument("fasta", type=str, help="FASTA file containing all contigs to derive gRNA sites for. All FASTA sequences should have unique primary headers (everything between the '>' and the first ' ' should be unique).")
    parser.add_argument("gff", type=str, help="GFF file specifying chromosomal features")
    
    # Add optional arguments
    parser.add_argument("--min_contig_edge_distance", metavar="N", type=int, default=500, help="Minimum distance from contig edge a site can be found")
    #parser.add_argument("--features", metavar="FEATURE,FEATURE", type=str, default="ORF", help="Features to design gRNA sites against")
    parser.add_argument("--features", metavar="FEATURE", type=str, nargs="+", default=["ORF"], help="Features to design gRNA sites against")
    parser.add_argument("--min_donor_length", metavar="N", type=int, default=80, help="The minimum length of the final computed donor DNA for each site")
    parser.add_argument("--max_donor_length", metavar="N", type=int, default=90, help="The maximum length of the final computed donor DNA for each site")
    parser.add_argument("--min_feature_edge_distance", metavar="N", type=int, default=24, help="The minimum distance a gRNA site can be from the edge of the feature. If negative, the maximum distance a gRNA site can be outside the feature.")
    parser.add_argument("--min_donor_insertions", metavar="N", type=int, default=2, help="The uniqueness of final donor DNA compared to the rest of the genome")
    parser.add_argument("--min_donor_deletions", metavar="N", type=int, default=2, help="The uniqueness of final donor DNA compared to the rest of the genome")
    parser.add_argument("--min_donor_mismatches", metavar="N", type=int, default=2, help="The uniqueness of final donor DNA compared to the rest of the genome")
    parser.add_argument("--min_donor_differences", metavar="N", type=int, default=3, help="The uniqueness of final donor DNA compared to the rest of the genome")
    parser.add_argument("--min_donor_distance", metavar="N", type=int, default=36, help="The minimum distance in bp a difference can exist from the edge of donor DNA")
    parser.add_argument("--min_target_length", metavar="N", type=int, default=23, help="The minimum length of the 'target'/'spacer'/gRNA site")
    parser.add_argument("--max_target_length", metavar="N", type=int, default=28, help="The maximum length of the 'target'/'spacer'/gRNA site")
    
    return parser.parse_args()

def read_fasta_file(args):
    """Load contig sequences from file into dict()
    Primary sequence headers must be unique
    """
    contigs = {}
    with open(args.fasta, 'r') as flo:
        name = None
        for line in flo:
            line = line.rstrip()
            if line.startswith('>'):
                name = regex.split(r'\s+', line[1:], 1)[0]
                contigs[name] = ''
            else:
                # Handle malformatted FASTA
                if ((name == None) or (name == "")):
                    raise ValueError('FASTA file malformatted')
                else:
                    contigs[name] += line
    
    return contigs

def lcs(string1, string2):
    """Find the longest common substring between two strings"""
    import difflib
    
    matcher = difflib.SequenceMatcher(None, string1, string2, True)
    match = matcher.find_longest_match(0, len(string1), 0, len(string2))
    # Match(a=0, b=15, size=9)
    return match
    # print(string1[match.a: match.a + match.size])
    # print(string2[match.b: match.b + match.size])

def generate_excise_target(args, feature):
    """Finds the gRNA sequence to cut within the specified places
    on the feature.
    """
    pass

def generate_revert_target(args, feature):
    """
    Creates the gRNA sequence that targets the excise donor DNA oligo
    """
    pass

def generate_excise_donor(args, feature, revert_target):
    """Creates the DNA oligo with the structure:
    homology--unique gRNA--homology
    that excises the target feature
    """
    pass

def generate_revert_donor(args, feature):
    """
    Use template DNA sequence to create oligo that will be used for fixing
    the DSB
    """
    # Optionally provide a list of restriction enzyme targeting sites
    # on either side for easier cloning
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

def main():
    # Get timestamp
    start = time.time()
    
    # Obtain command line arguments
    args = parse_arguments()
    print(args, file=sys.stderr)
    
    # Convert input files to memory structures
    contigs = read_fasta_file(args)
    
    
    # Print time taken for program to complete
    print(time.time()-start, file=sys.stderr)

if (__name__ == "__main__"):
    main()
