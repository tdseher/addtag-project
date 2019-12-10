#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/subroutines/subroutine.py

# Import standard packages
import sys
import os
import argparse
import textwrap

# Import included AddTag-specific modules
from .. import utils

# Define meta variables
__author__ = "Thaddeus D. Seher (@tdseher) & Aaron Hernday"
__summary__ = "CRISPR/Cas-directed HDR genome editing suite: finds+scores gRNA targets, generates donor DNAs, & produces optimal cPCR primer designs."
# '{__wrap_summary__}'.format(__wrap_summary__='\n  '.join(textwrap.wrap(__summary__)))
__date__ = utils.load_git_date()
__fullversion__ = utils.load_git_version()
__version__ = __fullversion__[:7]
__revision__ = utils.load_git_revision()
__program__ = os.path.basename(sys.argv[0])
__repository__ = 'https://github.com/tdseher/addtag-project.git'
#__citation__ = "{__author__}. AddTag. In preparation ({__date__})".format(**locals())
__citation__ = "Thaddeus D. Seher and Aaron D. Hernday. Addtag: Program for identifying exclusive endogenous gRNA sites and creating unique synthetic gRNA sites. University of California, Merced. Retrieved from <https://github.com/tdseher/addtag-project> (2019)."
__wrap_citation__ = '\n    '.join(textwrap.wrap(__citation__))
__copyright__ = "Copyright (c) 2016 {__author__}".format(**locals())
__license__ = '''This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Affero General Public License for more details.
  
  You should have received a copy of the GNU Affero General Public License
  along with this program. If not, see <https://www.gnu.org/licenses/>.
  
  Some features rely on code written by other people, provided under 
  different licenses. Please review them individually for more information.
'''.format(**locals())

class Subroutine():
    '''
    Class that serves as a template for a Subroutine that will be added to the program usage.
    '''
    def __init__(self, subparsers):
        self.subparsers = subparsers
        
        self.name = None
        self.description = None
        self.help = None
        self.epilog = None
        
        self.define_parser()
    
    def define_parser(self):
        self.parser = self.subparsers.add_parser(
            self.name,
            description=self.description,
            epilog=self.epilog,
            formatter_class=CustomHelpFormatter,
            help=self.help
        )
        self.parser.set_defaults(func=self.compute)
        
        # Change the help text of the "-h" flag
        self.parser._actions[0].help='Show this help message and exit.'
        
        # Special version action optional argument
        self.parser.add_argument("-v", "--version", action='version',
            help="Show program's version number and exit.",
            version='{__program__} {__version__} (revision {__revision__})'.format(**globals()))
    
    def define_arguments(self):
        '''
        Overload this method
        '''
        pass
    
    def compute(self, args):
        '''
        Overload this method
        '''
        pass

class CustomHelpFormatter(argparse.HelpFormatter):
    '''
    Help message formatter which retains any formatting in descriptions
    and adds default values to argument help.
    
    Only the name of this class is considered a public API. All the methods
    provided by the class are considered an implementation detail.
    '''
    # This class combines:
    #   argparse.ArgumentDefaultsHelpFormatter
    #   argparse.RawDescriptionHelpFormatter
    
    def _fill_text(self, text, width, indent):
        return ''.join([indent + line for line in text.splitlines(True)])
    
    def _get_help_string(self, action):
        help = action.help
        if '%(default)' not in action.help:
            if action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += ' (default: %(default)s)'
        return help

class ValidateCycleStop(argparse.Action):
    def __call__(self, parser, args, value, option_string=None):
        # print '{n} {v} {o}'.format(n=args, v=values, o=option_string)
        #f = '--'+self.dest.replace('_', '-')
        f = '--'+self.dest
        
        if (value == None):
            value = args.cycle_start
            
        elif (value < args.cycle_start):
            parser.error("argument {f}: must be greater than or equal to '--cycle_start {p}'".format(f=f, p=args.cycle_start))
        
        setattr(args, self.dest, value)

class ValidateFlanktags(argparse.Action):        
    def __call__(self, parser, args, values, option_string=None):
        # print '{n} {v} {o}'.format(n=args, v=values, o=option_string)
        f = '--'+self.dest.replace('_', '-')
        
        nmin = 1
        nmax = 2
        if not (nmin <= len(values) <= nmax):
            #raise argparse.ArgumentTypeError('argument --{f}: expected between {nmin} and {nmax} arguments'.format(f=self.dest, nmin=nmin, nmax=nmax))
            parser.error('argument {f}: expected between {nmin} and {nmax} arguments'.format(f=f, nmin=nmin, nmax=nmax))
        
        # If only 1 argument, then it is a filename
        if (len(values) == 1):
            # Look to see if the path exists and is a file (not a directory)
            if not os.path.isfile(values[0]):
                parser.error("argument {f}: no such file: '{p}'".format(f=f, p=values[0]))
        elif (len(values) == 2):
            # If 2 arguments, then the first is either "uniform" or "specific"
            if (values[0] not in ['uniform', 'specific']):
                parser.error("argument {f}: invalid design TYPE: '{t}' (choose from 'uniform', 'specific')".format(f=f, t=values[0]))
            
            if (vlaues[1] not in ['single', 'paired']):
                parser.error("argument {f}: invalid design TYPE: '{t}' (choose from 'single', 'paired')".format(f=f, t=values[0]))
        
        setattr(args, self.dest, values)

class ValidateKodDNA(argparse.Action):
    def __call__(self, parser, args, value, option_string=None):
        # print '{n} {v} {o}'.format(n=args, v=values, o=option_string)
        f = '--'+self.dest.replace('_', '-')
        
        if (value not in ['mintag', 'addtag', 'unitag', 'bartag']):
            # Look to see if the path exists and is a file (not a directory)
            if not os.path.isfile(value):
                parser.error("argument {f}: no such file: '{p}'".format(f=f, p=value))
        
        setattr(args, self.dest, value)

class ValidateInternalPrimersRequired(argparse.Action):
    def __call__(self, parser, args, values, option_string=None):
        # print '{n} {v} {o}'.format(n=args, v=values, o=option_string)
        f = '--'+self.dest
        
        if (values != None):
            if (2*len(args.dDNAs) != len(values)):
                parser.error('argument {f}: expected {n} arguments'.format(f=f, n=2*len(args.dDNAs)))
        
        setattr(args, self.dest, values)

class ValidatePrimersRequired(argparse.Action):
    def __call__(self, parser, args, values, option_string=None):
        # print '{n} {v} {o}'.format(n=args, v=values, o=option_string)
        f = '--'+self.dest
        if (values != None):
            if (len(args.dDNAs)+1 != len(values)):
                parser.error('argument {f}: expected {n} arguments'.format(f=f, n=len(args.dDNAs)+1))
        
        setattr(args, self.dest, values)

# class ValidateDesign(argparse.Action):        
#     def __call__(self, parser, args, values, option_string=None):
#         # print '{n} {v} {o}'.format(n=args, v=values, o=option_string)
#         nmin = 1
#         nmax = 2
#         if not (nmin <= len(values) <= nmax):
#             raise argparse.ArgumentTypeError('argument --{f}: expected between {nmin} and {nmax} arguments'.format(f=self.dest, nmin=nmin, nmax=nmax))
#         
#         s = set(values) # in case there are duplicates, reduce them to a single occurrence with set()
#         valid_values = {'cut', 'addtag', 'sigtag'} # set
#         
#         d = s.difference(valid_values)
#         if (len(d) > 0):
#             raise ValueError('invalid Design TYPE: %s (choose from {cut, addtag, sigtag})' % d.pop())
#         
#         setattr(args, self.dest, list(s))
# 
# class ValidateShowMotifs(argparse.Action):
#     def __call__(self, parser, args, values, option_string=None):
#         # print '{n} {v} {o}'.format(n=args, v=values, o=option_string)
#         if (values == []):
#             utils.print_local_file('motifs.txt')
#             sys.exit()
#             #raise ValueError('invalid cores N: %s (choose N > 0)' % value)
#             #raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
#         
#         setattr(args, self.dest, values)
# 
# class ValidateShowGlossary(argparse.Action):
#     def __call__(self, parser, args, values, option_string=None):
#         # print '{n} {v} {o}'.format(n=args, v=values, o=option_string)
#         if (values == []):
#             utils.print_local_file('glossary.txt')
#             sys.exit()
#             #raise ValueError('invalid cores N: %s (choose N > 0)' % value)
#             #raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
#         
#         setattr(args, self.dest, values)