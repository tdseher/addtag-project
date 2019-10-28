#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/subroutines/_subroutine_list_aligners.py

# Import included AddTag-specific modules
from . import subroutine
from .. import aligners

class ListAlignersParser(subroutine.Subroutine):
    def __init__(self, subparsers):
        self.subparsers = subparsers
        
        self.name = 'list_aligners'
        self.description = (
            "description:" "\n"
            "  Show list of all supported alignment programs." "\n"
        )
        self.help = "Show list of all supported alignment programs."
        self.epilog = (
            "example:" "\n"
            "  You can show details on all implemented aligners by running:" "\n"
            "   $ python3 {__program__} aligners" "\n"
        ).format(**subroutine.__dict__)
        
        self.define_parser()
        self.define_arguments()
    
    def define_arguments(self):
        pass
    
    def compute(self, args):
        """Print information about the supported aligners"""
        # Aligners intended to implement = ['addtag', 'blast+', 'blat', 'bowtie', 'bowtie2', 'bwa', 'cas-offinder']
        # Other aligners to consider: 'rmap', 'maq', 'shrimp2', 'soap2', 'star', 'rhat', 'mrsfast', 'stampy'
        for x in aligners.aligners:
            print('==========', x.name, '==========')
            print('    Authors:', x.authors)
            print('      Title:', x.title)
            print('    Journal:', x.journal)
            print('    Issuing:', x.issuing)
            print('       Year:', x.year)
            print('        doi:', x.doi)
            #print('   Citation:', x.citation)
            print('      Input:', x.input)
            print('     Output:', x.output)
            print('  Truncated:', x.truncated)
            print('')
        
        # End 'compute()'
        