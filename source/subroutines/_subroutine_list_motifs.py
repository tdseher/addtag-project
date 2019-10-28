#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/subroutines/_subroutine_list_motifs.py

# Import included AddTag-specific modules
from . import subroutine
from .. import utils

class ListMotifsParser(subroutine.Subroutine):
    def __init__(self, subparsers):
        self.subparsers = subparsers
        
        self.name = 'list_motifs'
        self.description = (
            "description:" "\n"
            "  Show list of common CRISPR/Cas SPACER>PAM arrangements, then exit." "\n"
        )
        self.help = "List common CRISPR/Cas SPACER>PAM arrangements."
        self.epilog = (
            "example:" "\n"
            "  You can list CRISPR/Cas motifs by running:" "\n"
            "   $ python3 {__program__} {__subroutine__}" "\n"
        ).format(**dict(list(subroutine.__dict__.items()) + list({"__subroutine__": self.name}.items()))) # key:value pairs in the latter dict will replace any instances in the prior dict
        
        self.define_parser()
        self.define_arguments()
    
    def define_arguments(self):
        pass
    
    def compute(self, args):
        """Print the list of common CRISPR/Cas motifs"""
        utils.print_local_file('motifs.txt')
        
        # End 'compute()'
        