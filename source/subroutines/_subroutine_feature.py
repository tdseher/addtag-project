#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/subroutines/_subroutine_feature.py

# Import standard packages
import sys

# Import non-standard packages
import regex

# Import included AddTag-specific modules
from . import subroutine
from .. import feature

class FeatureParser(subroutine.Subroutine):
    def __init__(self, subparsers):
        self.subparsers = subparsers
        
        self.name = 'feature'
        self.description = (
            "description:" "\n"
            "  Search GFF features for specific text." "\n"
        )
        self.help = "Search GFF features for specific text."
        self.epilog = (
            "example:" "\n"
            "  In general, you can search all feature attributes for text as follows:" "\n"
            "   $ python3 {__program__} feature --gff genome.gff --query HSP90 > features.gff" "\n"
        ).format(**subroutine.__dict__)
        #  If you want to limit your search to specific tags, then you could use
        #  these parameters:
        #   $ python3 {__program__} feature --gff genome.gff --query Gene=GAL4 > features.gff
        
        self.define_parser()
        self.define_arguments()
    
    def define_arguments(self):
        # Add mandatory arguments
        required_group = self.parser.add_argument_group('required arguments')
        required_group.add_argument("--gff", required=True, metavar="*.gff", type=str,
            help="GFF file specifying chromosomal features that will be searched.")
        
        required_group.add_argument("--query", required=True, nargs="+", metavar="TEXT",
            type=str, help="One or more words to search for within the GFF file.")
        
        # Add optional arguments
        self.parser.add_argument("--linked_tags", metavar="TAG", nargs="*",
            type=str, default=['ID', 'Name', 'Alias', 'Parent', 'Gene'],
            help="If a feature is found that matches the query, include other \
            features that have similar values for these tags.")
        
        self.parser.add_argument("--excluded_tags", metavar="TAG", nargs="*",
            type=str, default=['Note', 'orf_classification', 'parent_feature_type'],
            help="Do not search for within these tags.")
        
        self.parser.add_argument("--allow_errors", action="store_true",
            default=False, help="Include matches with minor differences from the query.")
        
        self.parser.add_argument("--header", action="store_true",
            default=False, help="Begin output with a commented line containing field names.")
    
    def compute(self, args):
        """
        Search GFF attribute field.
        Allow for fuzzy matches depending on the length of the input.
        If a match is found, then try to match the 'Gene', 'ID', 'Parent',
        etc tags to other features.
        Prints lines from input GFF that match.
        """
        
        matched_lines = set()
        includes = []
        
        # Search initial queries
        with open(args.gff, 'r') as flo:
            for i, line in enumerate(flo):
                line = line.rstrip()
                obj = feature.Feature.parse_gff_line(line)
                if obj:
                    for q in args.query:
                        if args.allow_errors:
                            errors = len(q)//4
                            q_pattern = '(?:'+q+'){e<='+str(errors)+'}'
                        else:
                            q_pattern = q
                        #for k in link_tags:
                        #    if k in obj.attributes:
                        for k in obj.attributes:
                            if k not in args.excluded_tags:
                                m = regex.search(q_pattern, obj.attributes[k], flags=regex.IGNORECASE)
                                if m:
                                    #matched_lines.add((line, 0))
                                    matched_lines.add(line)
                                    includes.append({ tag_key: obj.attributes[tag_key] for tag_key in args.linked_tags if tag_key in obj.attributes })
                                    #print('k={}, m={}'.format(k, m), file=sys.stderr)
                                    #print('includes[-1]={}'.format(includes[-1]), file=sys.stderr)
                                    break
        #for line in sorted(matched_lines):
        #    print(line)
        #for inc in includes:
        #    print(inc)
        # search linked attributes
        with open(args.gff, 'r') as flo:
            for i, line in enumerate(flo):
                line = line.rstrip()
                obj = feature.Feature.parse_gff_line(line)
                if obj:
                    for inc in includes:
                        for k, v in inc.items():
                            #errors = len(q)//4
                            #inc_pattern = '(?:'+v+'){e<='+str(errors)'}'
                            inc_patterns = [v]
                            if (k == 'Alias'):
                                inc_patterns = v.split(',')
                            
                            for inc_pattern in inc_patterns:
                                #if k in obj.attributes:
                                for obj_k, obj_v in obj.attributes.items():
                                    if ((obj_k not in args.excluded_tags) and (k not in args.excluded_tags)):
                                        m = regex.search(inc_pattern, obj_v, flags=regex.IGNORECASE)
                                        if m:
                                            #matched_lines.add((line, 1))
                                            matched_lines.add(line)
                                            #print('k={}, m={}'.format(k, m), file=sys.stderr)
                                            #print('inc_pattern={}'.format(inc_pattern), file=sys.stderr)
                                            break
        
        if args.header:
            print('# '+'\t'.join(['seqid', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']))
        
        for line in sorted(matched_lines):
            print(line)
        
        # End 'compute()'
        