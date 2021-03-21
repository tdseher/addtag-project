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

# Import included AddTag-specific modules
from . import utils
from . import nucleotides
from . import algorithms
from . import aligners
from . import thermodynamics
from . import bartag

from . import subroutines
from .subroutines import subroutine

from .donors import Donor, ExcisionDonor, ReversionDonor
from .targets import Target, ExcisionTarget, ReversionTarget
from .motifs import OnTargetMotif, OffTargetMotif
from .feature import Feature

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
  If you use AddTag for your research, please cite us. Because the manuscript
  is currently in preparation, you will need to cite the code repository
  instead.
  
    {__wrap_citation__}

copyright:
  AddTag {__copyright__}.

license:
  {__license__}

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
  folder/excision-targets.fasta     Target sequences with scores
  folder/excision-constructs.fasta  Construct sequences to be synthesized
  folder/excision-dDNAs.fasta       dDNA sequences to be synthesized for
                                    knock-out
  folder/reversion-query.fasta      FASTA file of candidate knock-in spacers
  folder/reversion-query.OUT        SAM/BLASTn alignment of candidate spacers
  folder/reversion-query.err        STDOUT/STDERR from alignment
  folder/reversion-targets.fasta    Target sequences with scores
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
""".format(**subroutine.__dict__) #.format(**locals())
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
""".format(**subroutine.__dict__) #.format(**locals())

class CustomFormatter(logging.Formatter):
    # logging.Formatter methods:
    #   self.formatTime(record, datefmt=None)
    #   self.formatException(ei)
    #   self.usesTime()
    #   self.formatMessage(record)
    #   self.formatStack(stack_info)
    #   self.format(record)
    
    def format(self, record):
        # logging.LogRecord attributes/methods:
        #   self.args
        #   self.name
        #   self.levelname
        #   self.levelno
        #   self.pathname
        #   self.filename
        #   self.module
        #   self.lineno
        #   self.created           # this is a time.time() value
        #   self.relativeCreated   # number milliseconds since created
        #   self.thread
        #   self.threadName
        #   self.process
        #   self.processName
        #   self.msg
        #   self.exec_info
        #   self.stack_info
        #   self.funcName
        #   self.getMessage()
        
        record.extname = '{name}.{funcName}()'.format(**record.__dict__)
        record.message = record.getMessage()
        #if self.usesTime():
        record.asctime = self.formatTime(record, self.datefmt)
        #s = self.formatMessage(record)
        
        return '[{asctime}] [{levelname:<8}] [{lineno:>4}] [{extname:<110}] {message}'.format(**record.__dict__)
        #else:
        #    return '[{levelname:<8}] [{extname:<110}] [{lineno:>4}] {message}'.format(**record.__dict__)

class Main(object):
    def __init__(self):
        """Create argument parser, then run the selected sub-program"""
        
        # Check to make sure that STDOUT redirects in Windows will work.
        # We restart the command if necessary
        self.encoding_restart()
        
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
            self.logger.info('{__program__} {__version__} (revision {__revision__})'.format(**subroutine.__dict__))
            self.logger.info(args)
        
        # call the function for which action was used
        args.func(args)
        
        if hasattr(args, 'folder'):
            # Print time taken for program to complete
            end_time = time.time()
            elapsed = end_time-start_time
            self.logger.info('{__program__} finished'.format(**subroutine.__dict__))
            self.logger.info('Start time: {}s'.format(start_time))
            self.logger.info('End time: {}s'.format(end_time))
            self.logger.info('Runtime: {}s'.format(elapsed))
            self.logger.info('Runtime: {}'.format(str(datetime.timedelta(seconds=elapsed))))
    
    def encoding_restart(self):
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
        # Sufficient memory must be available for loading and executing the new process. 
        # TODO: make '--no_encoding_restart' an option via the command-line, because it can mess things up
        # This means that if launch AddTag from another process
        # that needs to monitor its progress, then it will only function correctly if
        # you set the PYTHONIOENCODING variable before calling AddTag.
        # For example, running 'python3 -m cProfile addtag-project\addtag' will not work because of this.
        
        #print('Executable: {}'.format(sys.executable))
        #print('ARGV: {}'.format(sys.argv))
        #print('Parent process id: {}'.format(os.getppid()))
        #print('Current process id: {}'.format(os.getpid()))
        #try:
        #    print('Env var: {}'.format(os.environ['PYTHONIOENCODING']))
        #except KeyError:
        #    print('Env var: None')
        #print('')
        
        #if sys.platform.startswith('win'):
        if (sys.stdout.encoding != 'utf-8'):
            if ((sys.version_info >= (3, 7)) and hasattr(sys.stdout, 'reconfigure')):
                sys.stdout.reconfigure(encoding='utf-8')
            elif (('PYTHONIOENCODING' not in os.environ) or (os.environ['PYTHONIOENCODING'] != 'utf-8')):
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
        subparsers = parser.add_subparsers(metavar='action', help='Choose an action to perform (required)', dest='action')
        
        # Change the help text of the "-h" flag
        parser._actions[0].help='Show this help message and exit.'
        
        # Special version action optional argument
        parser.add_argument("-v", "--version", action='version',
            help="Show program's version number and exit.",
            version='{__program__} {__version__} (revision {__revision__})'.format(**subroutine.__dict__)) # .format(**globals()))
        
        # TODO: add '-i', '--info' option that shows the PATH of the python executable running AddTag
        #       as well as the PATH of the 'addtag' executable that is being run
        
        return parser, subparsers
    
    def parse_arguments(self):
        '''
        Create parser object, populate it with options, then parse command line
        input.
        '''
        # Create the main parser
        parser, subparsers = self._parser_general()

        # Make all subroutines
        subroutines.make_subroutines(subparsers)
        
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
            
            # TODO: If the folder DOES exist, then check if a 'log.txt' file exists already.
            #       If there is already a 'log.txt' file, then rename it to 'log-1.txt'
            #       Do this for every previous log file.
            
            # Create the logger, and have it write to 'folder/log.txt'
            #logging.basicConfig(filename=os.path.join(args.folder, 'log.txt'), level=logging.INFO, format='[%(asctime)s] [%(name)s.%(funcName)s] [%(lineno)d] %(message)s') # format='%(levelname)s %(asctime)s: %(message)s' #### %(pathname)s | %(filename)s | %(module)s
            #logging.basicConfig(filename=os.path.join(args.folder, 'log.txt'), level=logging.INFO, style='{', format='[{asctime}] [{name:<70}] [{funcName:<40}] [{lineno:>4}] {message}')
            logging.basicConfig(filename=os.path.join(args.folder, 'log.txt'), level=logging.INFO, style='{')
            formatter = CustomFormatter(style='{')
            root_logger = logging.getLogger(name=None)
            for h in root_logger.handlers:
                h.setFormatter(formatter)
        
        if hasattr(args, 'aligner'):
            # Add 'args.selected_aligner' to hold the actual aligner object
            for a in aligners.pw_aligners:
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
        
        # Actions don't activate when the argument retains its default value.
#        if hasattr(args, 'cycle_stop'):
#            if (args.cycle_stop == None):
#                args.cycle_stop = args.cycle_start
        
        if hasattr(args, 'weights'):
            for w in args.weights:
                w_name, w_pars = w.split(':')
                for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms:
                    if (C.name == w_name):
                        C.weight_str = w
                        C.weight_parameters = C.parse_weight(C.weight_str)
                        break
                else:
                    raise Exception('No Algorithm named {} exists.'.format(w_name))
                
    
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

