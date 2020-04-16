#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/subroutines/_subroutine_list_algorithms.py

# Import included AddTag-specific modules
from . import subroutine
from .. import algorithms

class ListAlgorithmsParser(subroutine.Subroutine):
    def __init__(self, subparsers):
        self.subparsers = subparsers
        
        self.name = 'list_algorithms'
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
            print('========== {} =========='.format(C.name))
            print('  Available: {}'.format(C.available))
            print('    Authors: {}'.format(C.authors))
            print('      Title: {}'.format(C.title))
            print('    Journal: {}'.format(C.journal))
            print('    Issuing: {}'.format(C.issuing))
            print('       Year: {}'.format(C.year))
            print('        doi: {}'.format(C.doi))
            print(' Off-target: {}'.format(C.off_target))
            print('  On-target: {}'.format(C.on_target))
            print('  Prefilter: {}'.format(C.prefilter))
            print(' Postfilter: {}'.format(C.postfilter))
            print('   Min, Max: {}, {}'.format(C.minimum, C.maximum))
            print('    Default: {}'.format(C.default))
            print('     Weight: {}'.format(C.weight_str))
            print('       RGNs: {}'.format(C.rgn_list))
            print('Description: {}'.format(C.description))
            print('')
        #prefilter_choices = [C for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.prefilter]
        #off_target_choices = [C for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.off_target]
        #on_target_choices = [C for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.on_target]
        #postfilter_choices = [C for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.postfilter]
        
        # End 'compute()'
        