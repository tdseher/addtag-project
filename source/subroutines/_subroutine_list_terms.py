#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/subroutines/_subroutine_list_terms.py

# Import included AddTag-specific modules
from . import subroutine
from .. import utils

class ListTermsParser(subroutine.Subroutine):
    def __init__(self, subparsers):
        self.subparsers = subparsers
        
        self.name = 'list_terms'
        self.description = (
            "description:" "\n"
            "  Show glossary of common CRISPR/Cas terms, then exit." "\n"
        )
        self.help = "Show glossary of common CRISPR/Cas terms."
        self.epilog = (
            "example:" "\n"
            "  You can show a list of CRISPR/Cas terms by running:" "\n"
            "   $ python3 {__program__} {__subroutine__}" "\n"
        ).format(**dict(list(subroutine.__dict__.items()) + list({"__subroutine__": self.name}.items()))) # key:value pairs in the latter dict will replace any instances in the prior dict
        
        self.define_parser()
        self.define_arguments()
    
    def define_arguments(self):
        pass
    
    def compute(self, args):
        """Print the Glossary"""
        utils.print_local_file('glossary.txt')
        
        # End 'compute()'
        