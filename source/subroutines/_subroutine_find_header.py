#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/subroutines/_subroutine_find_header.py

# Import non-standard packages
import regex

# Import included AddTag-specific modules
from . import subroutine

class ExtractParser(subroutine.Subroutine):
    def __init__(self, subparsers):
        self.subparsers = subparsers
        
        self.name = 'find_header'
        self.description = (
            "description:" "\n"
            "  Extracts selected sequences from input FASTA by matching their primary" "\n"
            "  sequence header (Everything between the '>' and the first whitespace), and" "\n"
            "  outputs in FASTA format." "\n"
        )
        self.help = "Search FASTA headers for specific text."
        self.epilog = (
            "example:" "\n"
            "  Try running AddTag with the following arguments:" "\n"
            "   $ python3 {__program__} {__subroutine__} --fasta excision-spacers.fasta" "\n"
            "   --query exTarget-33 exTarget-21 > extract.fasta" "\n"
        ).format(**dict(list(subroutine.__dict__.items()) + list({"__subroutine__": self.name}.items()))) # key:value pairs in the latter dict will replace any instances in the prior dict
        
        self.define_parser()
        self.define_arguments()
    
    def define_arguments(self):
        # Add mandatory arguments
        required_group = self.parser.add_argument_group('required arguments')
        required_group.add_argument("--fasta", required=True, nargs="+", metavar="*.fasta", type=str,
            help="FASTA files with contigs. All FASTA sequences should have unique \
            primary headers (everything between the '>' symbol and the first \
            whitespace should be unique).")
        
        required_group.add_argument("--query", required=True, nargs="+", metavar="TEXT",
            type=str, help="One or more words to search for within sequence headers.")
        
        # Add optional arguments
        self.parser.add_argument("--allow_errors", action="store_true",
            default=False, help="Include matches with minor differences from the query.")
    
    def compute(self, args):
        """
        Search input FASTA headers for arbitrary text
        Print out headers+sequences that match
        """
        
        # Build the regex pattern
        q_patterns = []
        for q in args.query:
            if args.allow_errors:
                errors = len(q)//4
                q_pattern = '(?:'+q+'){e<='+str(errors)+'}'
            else:
                q_pattern = q
            q_patterns.append(q_pattern)
        
        pattern = '|'.join(q_patterns)
        
        print_on = False
        for fn in args.fasta:
            with open(fn, 'r') as flo:
                for line in flo:
                    line = line.rstrip()
                    if (len(line) > 0):
                        if line.startswith('>'):
                            m = regex.search(pattern, line, flags=regex.IGNORECASE | regex.ENHANCEMATCH) # regex.BESTMATCH
                            if m:
                                print_on = True
                            else:
                                print_on = False
                        if print_on:
                            print(line)
        
        # End 'compute()'
        