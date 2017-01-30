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
    parser.add_argument("--min_target_length", metavar="N", type=int, default=23, help="The minimum length of the 'target'/'spacer'/gRNA site")
    parser.add_argument("--max_target_length", metavar="N", type=int, default=28, help="The maximum length of the 'target'/'spacer'/gRNA site")
    
    return parser.parse_args()

def read_fasta_file(args):
    contigs = {}
    with open(args.fasta, 'r') as flo:
        name = None
        for line in flo:
            line = line.rstrip()
            if line.startswith('>'):
                name = regex.split(r'\s+', line[1:], 1)[0]
                contigs[name] = ''
            else:
                contigs[name] += line
    
    return contigs

def process(args):
    # The Cas9 cuts 3-4bp upstream of the PAM sequence.
    
    # Cut site can be anywhere within target feature
    # genome          ACGGATTAGAGAGAGGCCTCCTCCCGGAT-GAATGGAAGACTAAACGGTAGATATAAG
    # feature         -------|ORF->                                <-ORF|-------
    # PAM                                              PAM         PAM
    # cut target                            CCCGGAT^GAATGG
    # cut genome      ACGGATTAGAGAGAGGCCTCCTCCCGGAT^GAATGGAAGACTAAACGGTAGATATAAG
    # donor           ACGGATT-----------------------------------AAACGGTAGATATAAG
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
