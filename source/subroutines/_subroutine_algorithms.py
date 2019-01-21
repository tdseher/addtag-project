#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/subroutines/_subroutine_algorithms.py

# Import included AddTag-specific modules
from . import subroutine
from .. import algorithms

class AlgorithmsParser(subroutine.Subroutine):
    def __init__(self, subparsers):
        self.subparsers = subparsers
        
        self.name = 'algorithms'
        self.description = (
            "description:" "\n"
            "  Show list of all implemented gRNA evaluation algorithms." "\n"
        )
        self.help = "Show list of all implemented gRNA evaluation algorithms."
        self.epilog = (
            "example:" "\n"
            "  You can show details on all implemented algorithms by running:" "\n"
            "   $ python3 {__program__} algorithms" "\n"
        ).format(**subroutine.__dict__)
        
        self.define_parser()
        self.define_arguments()
    
    def define_arguments(self):
        pass
    
    def compute(self, args):
        """Print information on the implemented algorithms"""
        for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms:
            print('==========', C.name, '==========')
            print('     Author:', C.author)
            print('       Year:', C.year)
            print('   Citation:', C.citation)
            print(' Off-target:', C.off_target)
            print('  On-target:', C.on_target)
            print('  Prefilter:', C.prefilter)
            print(' Postfilter:', C.postfilter)
            print('   Min, max:', C.minimum, C.maximum)
            print('    Default:', C.default)
            print('')
        #prefilter_choices = [C for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.prefilter]
        #off_target_choices = [C for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.off_target]
        #on_target_choices = [C for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.on_target]
        #postfilter_choices = [C for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.postfilter]
        
        # End 'compute()'
        