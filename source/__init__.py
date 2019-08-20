#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/__init__.py

# Import standard packages
import sys
import os
import argparse
import time
import datetime
import logging
import random
import math
import copy
import signal
from collections import namedtuple

# Import non-standard packages
import regex

# Import included AddTag-specific modules
from . import utils
from . import nucleotides
#from . import scores
from . import algorithms
from . import aligners
from . import thermodynamics
from . import bartag
#from . import evalues

from . import subroutines
#from . import subroutine
#from . import _subroutine_glossary
#from . import _subroutine_motifs
#from . import _subroutine_algorithms
#from . import _subroutine_aligners
#from . import _subroutine_thermodynamics
#from . import _subroutine_search
#from . import _subroutine_feature
#from . import _subroutine_extract
#from . import _subroutine_evaluate
#from . import _subroutine_generate
#from . import _subroutine_confirm

from .donors import Donor, ExcisionDonor, ReversionDonor
from .targets import Target, ExcisionTarget, ReversionTarget
from .motifs import OnTargetMotif, OffTargetMotif
from .feature import Feature

# Define meta variables
__author__ = "Thaddeus D. Seher (@tdseher) & Aaron Hernday"
__date__ = utils.load_git_date()
__fullversion__ = utils.load_git_version()
__version__ = __fullversion__[:7]
__revision__ = utils.load_git_revision()
__program__ = os.path.basename(sys.argv[0])
__citation__ = "{__author__}. AddTag. Unpublished ({__date__})".format(**locals())
__description__ = """\
description:
  Program for identifying exclusive endogenous gRNA sites and creating unique
  synthetic gRNA sites.
  
  AddTag can produce single or dual gRNA targeting 5' exon or essential protein
  domains. Excision (knock-out) and reversion (knock-in).
  
  It is important to note that although the on- and off-target scores are
  provided for other nucleases (SaCas9, Cpf1, etc), the algorithms used
  were trained on SpCas9. Therefore, the scores provided likely do not
  reflect the behavior and specificity of nucleases other than SpCas9.
  AddTag can still be used to design guides for these nucleases based on
  complementarity and enzyme-dependent PAM sites, but the efficacy and
  specificity of these guides are unpredictable in this version.
  
  Diagram of DNA-RNA hybridization and nuclease catalyzed by Cas9:
    (sgRNA = crRNA + linker + tracrRNA) = (gRNA = spacer + scaffold)
                                       gRNA ┌────────────────────────┐
                                      sgRNA ┌───tracrRNA────┐┌linker┐│
                            cut┐   ┌┬┬┬┐  3'╤╤╤╤╗            ╔╤╤╗   ││
           ┌──╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥ ╥╥╥┘   │        ╚╦╦╦╦╦╦╦╦╦╦╦╦╝  ╢   ││
           │5'╩╩╩╩╩╩╩╩╩╩╩╩╩╩╩╩╩═╩╩╩╧╧╧╧╪╧╧╧╧╧╧╧╧╧╩╩╩╩╩╩╩╩╩╩╩╩╧╧╧╝   ││
           │  └────────────────────────│────────crRNA───────┘└──────┘│
           │  └──────spacer───────┘└───│─────────────scaffold────────┘
           │                     Cas┬─┐│
    3'╥╥╥╥╥┘                cut┐    xxx└╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥5'
    5'╨╨╨╨╨───┴┴┴┴┴┴┴┴┴┴┴┴┴┴┴┴┴ ┴┴┴─╨╨╨─╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨3' genome
              └──────target───────┘ └─┴PAM
  
  Diagram of genome labels for excision:
                us_feature_trim┐           ┌ds_feature_trim
   ─genome┐┌──ex_us_homology─┐┌┴┐┌feature┐┌┴┐┌─ex_ds_homology──┐┌genome─
   ╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥
   ╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨
                           ex_target┴─┘
  
  Diagram of excision dDNA labels:
              ┌──ex_us_homology─┐┌insert┐┌─ex_ds_homology──┐
              ╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥
              ╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨
              └──────────────ex_donor_length───────────────┘
  
  Diagram of genome labels for reversion:
  ─genome┐┌────re_us_homology───┐┌insert┐┌───re_ds_homology────┐┌genome─
  ╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥
  ╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨
                              └──re_target─┘
  
  Diagram of reversion dDNA labels:
       ┌────re_us_homology───┐   ┌feature┐   ┌───re_ds_homology────┐
       ╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥
       ╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨
       └──────────────────────re_donor_length──────────────────────┘

version:
  short    {__version__}
  full     {__fullversion__}
  revision {__revision__}
  date     {__date__}

citation:
  {__citation__}

copyright:
  AddTag Copyright (c) 2017 {__author__}.
  All rights reserved.

license:
  You may not hold the authors liable for damages or data loss regarding the
  use or inability to use the software. AddTag is provided without any
  warranty. You may use the software for academic or private purposes without
  purchasing a license. However, using the software for commercial purposes
  requires purchase of a license. You may modify and distribute the software
  as long as you make the modifications open source, and grant patent license
  to the authors for inclusion in AddTag.
  
  Some AddTag features rely on code written by other people, provided under
  different licenses. Please review them individually for more information.

protein:
  The Cas9 or Cpf1 protein you use should be engineered specifically for your
  organism. It should be codon-optomized, and if using eukarya, contain an
  appropriate nuclear localization sequence. Cas9 from different species bind
  to different PAM sequences, which is useful when no suitable PAM sequence is
  present within your gene of interest. Additionally, the different Cas9 gene
  sequences can have huge length differences. Remember that each Cas9 is only
  compatible with the tracrRNA and crRNA (or synthetic gRNA) derived from the
  same species.

outputs:
  STDOUT                            Tab-delimited analysis results
  STDERR                            Errors
  folder/                           Folder holding generated files
  folder/log.txt                    Log of important program operations
  folder/aligner-index/*            Index of input FASTA created by aligner
  folder/excision-query.fasta       FASTA file of candidate knock-out spacers
  folder/excision-query.OUT         SAM/BLASTn alignment of candidate spacers
  folder/excision-query.err         STDOUT/STDERR from alignment
  folder/excision-spacers.fasta     Spacer sequences with scores
  folder/excision-constructs.fasta  Construct sequences to be synthesized
  folder/excision-dDNAs.fasta       dDNA sequences to be synthesized for
                                    knock-out
  folder/reversion-query.fasta      FASTA file of candidate knock-in spacers
  folder/reversion-query.OUT        SAM/BLASTn alignment of candidate spacers
  folder/reversion-query.err        STDOUT/STDERR from alignment
  folder/reversion-spacers.fasta    Spacer sequences with scores
  folder/reversion-constructs.fasta Construct sequences to be synthesized
  folder/reversion-dDNAs.fasta      dDNA sequences to be synthesized for
                                    knock-in
  folder/reversion-primers.fasta    Primers for amplifying knock-in dDNAs
  folder/protection-dDNAs.fasta     wt dDNAs of off-target sites to prevent
                                    mutations at off-target Cas9 binding sites
                                    (contains edit distance in header s/i/d)
  folder/protection-primers.fasta   Primers for amplifying protection DNAs
  folder/primers.fasta              Primers to check for KO/KI (contains
                                    expected amplicon sizes in header)
""".format(**locals())
__epilog__ = """\
example:
  At a minimum, AddTag needs to be executed with the '--fasta', '--gff', and
  '--folder' options:
   $ python3 {__program__} --fasta genome.fasta --gff genome.gff
     --folder output
  
  You will often want to save the analysis results. You can do so by redirecting
  STDOUT to a file. Here is another example usage with a few more parameters:
   $ python3 {__program__} --fasta chromosomes.fasta --gff features.gff
     --feature_homologs homologs.txt --excise_insert_lengths 0 3
     --excise_downstream_homology 47 50 --folder analysis > analysis.out
     2> analysis.err
""".format(**locals())



class TempPrimer(object):
    sequences = {} # key = nucleotide sequence, value = Primer object

class Main(object):
    def __init__(self):
        """Create argument parser, then run the selected sub-program"""
        
        # Check to make sure that STDOUT redirects in Windows will work.
        # We restart the command if necessary
        self.win_restart()
        
        # Obtain command line arguments and parse them
        args = self.parse_arguments()
        
        self.process_general_arguments(args)
        
        # Get timestamp for analysis beginning
        start_time = time.time()
        
        # Echo the command line parameters to STDOUT
        #print(args, flush=True)
        
        # Print command line parameters (arguments) to log file
        if hasattr(args, 'folder'):
            # Create the logger
            global logger
            logger = logging.getLogger(__name__)
            Main.logger = logger.getChild('Main')
            self.logger.info(args)
        
        # call the function for which action was used
        args.func(args)
        
        if hasattr(args, 'folder'):
            # Print time taken for program to complete
            end_time = time.time()
            elapsed = end_time-start_time
            self.logger.info('{} finished'.format(__program__))
            self.logger.info('Start time: {}s'.format(start_time))
            self.logger.info('End time: {}s'.format(end_time))
            self.logger.info('Runtime: {}s'.format(elapsed))
            self.logger.info('Runtime: {}'.format(str(datetime.timedelta(seconds=elapsed))))
    
    def win_restart(self):
        '''
        Check to make sure that STDOUT redirects in Windows will work.
        The environmental variable PYTHONIOENCODING needs to be set to 'utf-8'.
        If PYTHONIOENCODING has not been set or has a different value, then
        we restart the program, replacing the PID.
        '''
        # In Windows, exec() does not really replace the current process.
        # It creates a new process (with a new pid), and exits the current one.
        # When a call to an _exec function is successful, the new process is
        # placed in the memory previously occupied by the calling process.
        # Sufficient memory must be available for loading and executing the new
        # process. 
        
        #print('Executable: {}'.format(sys.executable))
        #print('ARGV: {}'.format(sys.argv))
        #print('Parent process id: {}'.format(os.getppid()))
        #print('Current process id: {}'.format(os.getpid()))
        #try:
        #    print('Env var: {}'.format(os.environ['PYTHONIOENCODING']))
        #except KeyError:
        #    print('Env var: None')
        #print('')
        
        if sys.platform.startswith('win'):
            if (('PYTHONIOENCODING' not in os.environ) or (os.environ['PYTHONIOENCODING'] != 'utf-8')):
                #print('Restarting')
                os.environ['PYTHONIOENCODING'] = 'utf-8'
                os.execvp(sys.executable, [sys.executable] + sys.argv)
    
    def calculate_amplicons(self, args, primer_sets, contigs):
        """
        Calculates number of amplicons for each primer set
        """
        greek_labels = {
            1: 'haploid', # mono
            2: 'diploid',
            3: 'triploid',
            4: 'tetraploid',
            5: 'pentaploid',
            6: 'hexaploid',
            7: 'heptaploid', # septaploid
            8: 'octaploid',
            9: 'ennaploid',
            10: 'decaploid',
        }
        
        #output_list = []
        for ps in primer_sets:
            # [[('Ca22chr3A_C_albicans_SC5314-r0[exDonor-42]', 943006, 943997, 991), ('Ca22chr3B_C_albicans_SC5314-r0[exDonor-42]', 942984, 943975, 991)], None, [], None, []]
            # uf_dr_pair, uf_fr_pair, uf_ir_pair, ff_dr_pair, if_dr_pair = ps
            amp_sets = []
            for pp in ps:
                if (pp != None):
                    amp_sets.append(args.selected_oligo.simulate_amplification(pp, contigs))
                else:
                    amp_sets.append(None)
            
            amp_lengths = [len(a) for pp, a in zip(ps, amp_sets) if pp]
            if (len(set(amp_lengths)) == 1):
                for n, label in greek_labels.items():
                    #if all([ len(a)==n for pp, a in zip(ps, amp_sets) if pp ]):
                    if (amp_lengths[0] == n):
                        ploidy = label + ' (' + str(n) + 'n)'
                        break
                else:
                    ploidy = 'polyploid' + ' (' + str(amp_lengths[0]) + 'n)'
            else:
                ploidy = 'ambiguous'
                
            # # Requires there to be only 1 amplification
            # if all([ len(a)==1 for pp, a in zip(ps, amp_sets) if pp ]):
            #     ploidy = 'haploid'
            # 
            # # Requires there to be only 2 amplifications
            # elif all([ len(a)==2 for pp, a in zip(ps, amp_sets) if pp ]):
            #     ploidy = 'diploid'
            #     
            # # Allows for any number of amplifications
            # else:
            #     ploidy = 'polyploid'
            
            print((ploidy, amp_sets), flush=True)
            #output_list.append((ploidy, amp_sets))
        
        #return output_list
    
    def _parser_general(self):
        '''general parser'''
        # Create the parent argument parser
        parser = argparse.ArgumentParser(
            description=__description__,
            epilog=__epilog__,
            formatter_class=subroutines.CustomHelpFormatter
        )
        
        # Create the subparsers
        subparsers = parser.add_subparsers(metavar='action', help='choose an action to perform (required)', dest='action')
        
        # Change the help text of the "-h" flag
        parser._actions[0].help='Show this help message and exit.'
        
        # Special version action optional argument
        parser.add_argument("-v", "--version", action='version',
            help="Show program's version number and exit.",
            version='{__program__} {__version__} (revision {__revision__})'.format(**globals()))
        
        return parser, subparsers
    
    def parse_arguments(self):
        '''
        Create parser object, populate it with options, then parse command line
        input.
        '''
        
        parser, subparsers = self._parser_general()
        #parser_glossary = self._parser_glossary(subparsers)
        #parser_motifs = self._parser_motifs(subparsers)
        #parser_algorithms = self._parser_algorithms(subparsers)
        #parser_aligners = self._parser_aligners(subparsers)
        #parser_oligos = self._parser_oligos(subparsers)
        #parser_search = self._parser_search(subparsers)
        #parser_feature = self._parser_feature(subparsers)
        #parser_extract = self._parser_extract(subparsers)
        #parser_evaluate = self._parser_evaluate(subparsers)
        #parser_generate = self._parser_generate(subparsers)
        #parser_confirm = self._parser_confirm(subparsers)
        
        # 'subroutine.py' should automatically search all '_subroutine_*.py' files for 'Subroutine' subclasses, then add them to 'subparsers'
        #parser_glossary = _subroutine_glossary.GlossaryParser(subparsers)
        #parser_motifs = _subroutine_motifs.MotifsParser(subparsers)
        #parser_algorithms = _subroutine_algorithms.AlgorithmsParser(subparsers)
        #parser_aligners = _subroutine_aligners.AlignersParser(subparsers)
        #parser_oligos = _subroutine_thermodynamics.ThermodynamicsParser(subparsers)
        #parser_search = _subroutine_search.SearchParser(subparsers)
        #parser_feature = _subroutine_feature.FeatureParser(subparsers)
        #parser_extract = _subroutine_extract.ExtractParser(subparsers)
        #parser_evaluate = _subroutine_evaluate.EvaluateParser(subparsers)
        #parser_generate = _subroutine_generate.GenerateParser(subparsers)
        #parser_confirm = _subroutine_confirm.ConfirmParser(subparsers)
        
        subroutines.make_subroutines(subparsers)
        
        
        
        
        # Add special arguments for additional help messages: These have been deprecated
        #parser.add_argument("-s", "--show", nargs=0, action=ValidateShowMotifs, default=argparse.SUPPRESS,
        #    help="Show list of common RGN motifs, then exit.")
        #parser.add_argument("-g", "--glossary", nargs=0, action=ValidateShowGlossary, default=argparse.SUPPRESS,
        #    help="Show glossary of common CRISPR/RGN terms, then exit.")
        
        # Add optional arguments: Deprecated
        #parser.add_argument("--design", nargs='+', type=str, action=ValidateDesign, default=['cut','addtag'],
        #    help="Choose 1 or 2 from {cut, addtag, sigtag}")
        
        
        #parser.add_argument("--python2_path", type=str, default="python",
        #    help="Path to the Python 2.7+ program")
        #parser.add_argument("--bowtie_path", type=str, default="bowtie",
        #    help="Path to the 'bowtie' executable")
        #parser.add_argument("--bowtie-build_path", type=str, default="bowtie-build",
        #    help="Path to the 'bowtie-build' executable")
        #parser.add_argument("--bowtie2_path", type=str, default="bowtie2",
        #    help="Path to the 'bowtie2' executable")
        #parser.add_argument("--bowtie2-build_path", type=str, default="bowtie2-build",
        #    help="Path to the 'bowtie2-build' executable")
        #parser.add_argument("--bwa_path", type=str, default="bwa",
        #    help="Path to the 'bwa' executable")
        #parser.add_argument("--blastn_path", type=str, default="blastn",
        #    help="Path to the 'blastn' executable")
        #parser.add_argument("--blat_path", type=str, default="blat",
        #    help="Path to the 'blat' executable")
        
        # Parse the arguments
        args = parser.parse_args()
        
        # If no arguments provided, then print the mini-help usage
        if (len(sys.argv) < 2):
            parser.print_usage()
            sys.exit(1)
        
        # Run Action on arguments with default values
        # Does not check for action 'conflicts', and does not update the
        # 'seen' actions list
        for sub in subroutines.subroutines:
            if (sub.name == args.action):
                subparser = sub.parser
                for action in subparser._actions:
                    if hasattr(args, action.dest):
                        if (hasattr(action, 'process_action_on_default') and (getattr(action, 'process_action_on_default') == True)):
                            v = getattr(args, action.dest)
                            if (v == action.default):
                                #sys.stderr.write("YES  " + str(action) + "\n")
                                action(subparser, args, v)
                            #else:
                            #    sys.stderr.write("NO   " + str(action) + "\n")
        
        
        # Return the parsed arguments
        return args
    
    def process_general_arguments(self, args):
        """Perform ubiquitous argument parsing things"""
        if hasattr(args, 'folder'):
            # Create the project directory if it doesn't exist
            os.makedirs(args.folder, exist_ok=True)
        
            # Create the logger, and have it write to 'folder/log.txt'
            logging.basicConfig(filename=os.path.join(args.folder, 'log.txt'), level=logging.INFO, format='[%(asctime)s] [%(name)s.%(funcName)s (%(lineno)d)] %(message)s') # format='%(levelname)s %(asctime)s: %(message)s' #### %(pathname)s | %(filename)s | %(module)s
        
        if hasattr(args, 'aligner'):
            # Add 'args.selected_aligner' to hold the actual aligner object
            for a in aligners.aligners:
                if (a.name == args.aligner):
                    args.selected_aligner = a
                    break
        
        if hasattr(args, 'oligo'):
            # Add 'args.selected_oligo' to hold the actual oligo object
            for o in thermodynamics.oligos:
                if (o.name == args.oligo):
                    args.selected_oligo = o
                    break
        
        ############################# Motif stuff ##############################
        # Old colde for compiling regex for motifs
        ## Note that any logging these functions perform will not be written to disk,
        ## as no 'handler' has been specified.
        #args.parsed_motifs = []
        #args.compiled_motifs = []
        #for motif in args.motifs:
        #    spacers, pams, side = self.parse_motif(motif) # Parse the motif
        #    args.parsed_motifs.append((spacers, pams, side)) # Add to args
        #    args.compiled_motifs.append(nucleotides.compile_motif_regex(spacers, pams, side, anchored=False)) # Add to args
        
        # populate --off_target_motifs with --motifs if None
        #if (args.off_target_motifs == None):
        #    args.off_target_motifs = args.motifs
        
        if hasattr(args, 'motifs'):
            # Parse the on-target motifs
            for motif in args.motifs:
                OnTargetMotif(motif)
        
        if hasattr(args, 'off_target_motifs'):
            # Parse the off-target motifs
            for motif in set(args.off_target_motifs).difference(args.motifs): # off-target motifs don't overlap with on-target ones
                OffTargetMotif(motif)
        
        # Motif.motifs           list of ALL motifs, both on-target and off-target?  <-- not implemented yet
        # OnTargetMotif.motifs   list of on-target motifs
        # OffTargetMotif.motifs  list of off-target motifs
        ########################################################################
        
        #if (hasattr(args, 'internal_primers_required') and hasattr(args, 'dDNAs')):
        #    if args.internal_primers_required:
        #        if (2*len(args.dDNAs) != len(args.internal_primers_required)):
        #            parser.error('argument {f}: expected {n} arguments'.format(f='--internal_primers_required', n=2*len(args.dDNAs)))
    
    def filter_features(self, features, selection):
        """
        Reduce the total number of features to just the ones indicated in the selection
        """
        # Current implementation limitations
        # GFF file only has 'Gene=' tag on ONE of the homologs, and not the other
        # User will have to specify the feature twice if he wants to target both homologs
        
        # Require at least one in order to filter
        # Feature has this format: key=gene/tag, value = (contig, start(bp), end(bp), strand)
        if selection:
            new_features = {}
            for s in selection:
                f = features.get(s)
                if f:
                    new_features[s] = f
                else:
                    raise Exception("Feature specified as '--selection "+s+"' does not exist in input '--gff'.")
                    #sys.exit(1)
            if (len(new_features) == 0):
                raise Exception("Input '--gff' and '--selection' combination returns no valid features.")
                #sys.exit(1)
            return new_features
        else:
            if (len(features) == 0):
                raise("Input '--gff' has no valid features.")
                #sys.exit(1)
            return features

