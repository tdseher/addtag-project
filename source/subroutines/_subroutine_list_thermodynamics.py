#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/subroutines/_subroutine_list_thermodynamics.py

# Import included AddTag-specific modules
from . import subroutine
from .. import thermodynamics

class ListThermodynamicsParser(subroutine.Subroutine):
    def __init__(self, subparsers):
        self.subparsers = subparsers
        
        self.name = 'list_thermodynamics'
        self.description = (
            "description:" "\n"
            "  Show list of all supported oligonucleotide thermodynamics property programs." "\n"
        )
        self.help = "Show list of all supported oligonucleotide thermodynamics property programs."
        self.epilog = (
            "example:" "\n"
            "  You can show details on all implemented aligners by running:" "\n"
            "   $ python3 {__program__} {__subroutine__}" "\n"
        ).format(**dict(list(subroutine.__dict__.items()) + list({"__subroutine__": self.name}.items()))) # key:value pairs in the latter dict will replace any instances in the prior dict
        
        self.define_parser()
        self.define_arguments()
    
    def define_arguments(self):
        pass
    
    def compute(self, args):
        """Print information about the supported thermodynamics programs"""
        for x in thermodynamics.oligos:
            print('========== {} =========='.format(x.name))
            print('  Available: {}'.format(x.available))
            print('    Authors: {}'.format(x.authors))
            print('      Title: {}'.format(x.title))
            print('    Journal: {}'.format(x.journal))
            print('    Issuing: {}'.format(x.issuing))
            print('       Year: {}'.format(x.year))
            print('        doi: {}'.format(x.doi))
            print('')
        
        # End 'compute()'
        