#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/__init__.py

# Import standard packages
import sys
import os
import argparse
import time
import logging
import random
import pickle
import math
import copy
from collections import namedtuple

# Import non-standard packages
import regex

# Import included AddTag-specific modules
from . import utils
from . import nucleotides
from . import scores
from . import algorithms
from . import aligners
from . import oligos
from . import bartag
#from . import evalues

# Create the logger
#logger = logging.getLogger(__name__)

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

class Alignment(object):
    """Class representing an alignment"""
    __slots__ = [
        'sequence',
        'target',
        'pam',
        'motif',
        'contig',
        'start',
        'end',
        'orientation',
        'upstream',
        'downstream',
        'postfilter',
        'score',
        'action',
    ]
    
    def __init__(self, sequence, target, pam, motif, contig, start, end, orientation, upstream='', downstream=''):
        # Most attributes derived from SAM output
        self.sequence = sequence
        self.target = target
        self.pam = pam
        self.motif = motif
        self.contig = contig
        self.start = start
        self.end = end
        self.orientation = orientation
        
        self.upstream = upstream
        self.downstream = downstream
        
        # Variable to hold whether or not this alignment should be included
        # in the off-target scoring
        self.postfilter = None
        
        # Variable to hold the scores
        self.score = {}
        self.action = 'None'
        
        # Haven't yet added:
        #  lidentities, ridentities, r_score, bae, chari, oof, proxgc, want, xu
        #self.ridentities = nucleotides.ridentities(self.contig_target, aligned_target)
        #self.r_scores = {}
        #for i in [4, 8, 12, 16]:
        #    self.r_scores[i] = scores.r_score(self.contig_target, aligned_target, i)
    
    def calculate_scores(self, parent):
        # parent = (sequence, target, pam, upstream, downstream)
        this = (self.sequence, self.target, self.pam, self.upstream, self.downstream)
        postfilter = []
        for C in algorithms.single_algorithms:
            c_score = C.calculate(this, parsed_motif=self.motif.parsed_list)
            self.score[C.name] = c_score
            if C.postfilter:
                if (C.minimum <= c_score <= C.maximum):
                    postfilter.append(True)
                else:
                    postfilter.append(False)
        for C in algorithms.paired_algorithms:
            c_score = C.calculate(parent, this, parsed_motif=self.motif.parsed_list)
            self.score[C.name] = c_score
            if C.postfilter:
                if (C.minimum <= c_score <= C.maximum):
                    postfilter.append(True)
                else:
                    postfilter.append(False)
        
        # If alignment meets all postfilter criteria, then set it as True
        self.postfilter = all(postfilter) # all([]) returns True
    
    def get_upstream_sequence(self, length, contigs):
        """Returns upstream sequence"""
        return ''
    
    def get_downstream_sequence(self, length, contigs):
        """Returns downstream sequence"""
        return ''
    
    def __repr__(self):
        return self.__class__.__name__ + '(' + ' '.join([
            self.target + '|' + self.pam,
            self.contig + ':' + self.orientation + ':' + str(self.start) + '..' + str(self.end),
            'action=' + str(self.action),
            'motif=' + self.motif.motif_string,
            'postfilter=' + str(self.postfilter)] +
            [x + '=' + str(round(self.score[x], 2)) for x in self.score]
            ) + ')'

class Donor(object):
    prefix = 'Donor'
    sequences = {}
    indices = {}
    
    @classmethod
    def get_contig_dict(cls):
        contigs = {}
        for sequence, obj in cls.sequences.items():
            contigs[obj.name] = sequence
        return contigs
    
    @classmethod
    def generate_fasta(cls, filename, sep=':'):
        """
        Creates a FASTA file
        >id feature:contig:orientation:start1..end1:mAT:start2..end2 ...
        """
        with open(filename, 'w') as flo:
            for sequence, obj in sorted(cls.sequences.items(), key=lambda x: int(x[1].name.split('-')[1])):
            #for sequence, obj in cls.sequences.items():
                #don_entry = tuple(["dDNA-"+str(i), feature, contig, '+', start1, start2, end1, end2, dDNA])
                print(' '.join(['>'+obj.name, 'spacers='+str(len(obj.spacers))] + sorted([obj.format_location(x, sep) for x in obj.locations])), file=flo)
                print(sequence, file=flo)
        logging.info(cls.__name__ + ' dDNA FASTA generated: {!r}'.format(filename))
        return filename
    
    def format_location(self, location, sep=':'):
        feature, contig, orientation, *segments = location
        output_list = [feature, contig, orientation]
        for x in segments:
            if isinstance(x, str):
                x = (x,)
            output_list.append('..'.join(map(str, x)))
        return sep.join(output_list)
    
    def get_features(self):
        """Return (sorted) list of all feature names this Donor maps to"""
        return sorted(set(x[0] for x in self.locations))
    
    def get_location_features(self):
        """Return (sorted) list of all feature names this Donor maps to"""
        return sorted(set(x[0] for x in self.locations))
    
    def __init__(self, feature, contig, orientation, sequence, *segments, spacer=None):
    #def __init__(self, feature, contig, orientation, segment1, insert, segment2, sequence, spacer=None):
        # List of all genomic locations and features this dDNA corresponds to
        location = (feature, contig, orientation, *segments)
        self.locations = set()
        
        #self.contig = contig, # The name of the contig (not the sequence)
        #self.orientation = orientation
        #self.segment1 = segment1 # (start1, end1)
        #self.insert = insert
        #self.segment2 = segment2 # (start2, end2)
        self.sequence = sequence
        
        # List to hold all spacer objects that map to this dDNA
        self.spacers = set()
        
        # Get the index number
        self.index = len(self.sequences)
        self.name = self.prefix + '-' + str(self.index)
        
        # Add this to the non-redundant list of excision dDNAs if it doesn't already exist
        # otherwise, add this locus to the pre-existing one
        the_dDNA = self.sequences.setdefault(self.sequence, self)
        if spacer:
            the_dDNA.spacers.add(spacer)
        the_dDNA.locations.add(location)
        self.indices.setdefault(the_dDNA.name, the_dDNA)
        
    # >dDNA-0 C2_10210C_B:Ca22chr2B_C_albicans_SC5314:+:2096471..2096521:CCA:2097397:2097444
#    def build_sequence(self, contig_sequence, orientation, segment1, segment2, insert):
#        """Takes substrings of contig to generate the sequence"""
#        sequence = contig_sequence[segment1[0]:segment1[1]] + insert + contig_sequence[segment2[0]:segment2[1]]
#        if (orientation == '-'):
#            print(self.name, "has '-' orientation", file=sys.stderr)
#        else:
#            pass
#        return sequence
    
    #def __repr__(self):
    #    return self.__class__.__name__ + '(' + ' '.join([
    #        self.name,
    #        self.contig + ':' + self.orientation + ':' +
    #        str(self.segment1[0]) + '..' + str(self.segment1[1]) + ':' +
    #        self.insert + ':' + str(self.segment2[0]) + '..' + str(self.segment2[1]),
    #        'features=' + ','.join(self.features) or 'None',
    #        'spacers=' + str(len(self.spacers))
    #        ]) + ')'

    def __repr__(self):
        """Return a string containing a printable representation of the Target object."""
        return self.__class__.__name__ + '(' + ' '.join(
            [self.name, 'spacers='+str(len(self.spacers))] +
            [self.format_location(x) for x in sorted(self.locations)]
            ) + ')'

class ExcisionDonor(Donor):
    prefix = 'exDonor'
    sequences = {}
    indices = {}
    
    # @classmethod
    # def get_targets(cls, args, sequence):
    #     targets = set()
    #     #for seq_i, sequence in enumerate(dDNAs):
    #     for orientation in ['+', '-']:
    #         if (orientation == '-'):
    #             sequence = nucleotides.rc(sequence)
    #         
    #         #for i in range(len(args.parsed_motifs)):
    #         for mymotif in OnTargetMotif.motifs:
    #             #spacers, pams, side = args.parsed_motifs[i]
    #             spacers, pams, side = mymotif.parsed_list
    #             #compiled_regex = args.compiled_motifs[i]
    #             #matches = nucleotides.motif_search(sequence, spacers, pams, side)
    #             matches = nucleotides.motif_search2(sequence, side, mymotif.compiled_regex)
    #             for seq, start, end, spacer, pam in matches:
    #                 if (orientation == '-'):
    #                     start, end = len(sequence) - end, len(sequence) - start
    #                 t_upstream = sequence[start-10:start]
    #                 t_downstream = sequence[end:end+10]
    #                 filtered_targets = target_filter(seq, spacer, pam, t_upstream, t_downstream, args)
    #                 for filt_seq, filt_spacer, filt_pam in filtered_targets:
    #                     targets.add((orientation, start, end, t_upstream, t_downstream, filt_seq, side, filt_spacer, filt_pam, mymotif.motif_string, tuple([tuple(x) if isinstance(x, list) else x for x in mymotif.parsed_list])))
    #     return sorted(targets) # becomes a list
    
    @classmethod
    def mintag_exhaustive_site_search(cls, args, f, contig_sequence, orientation):
        """
        Method for generating all potential mAT given the us/ds trims and
        insert size. Also called "Brute force" method
        
        This code is run once for each feature
        """
        for us_trim in range(args.excise_upstream_feature_trim[0], args.excise_upstream_feature_trim[1]+1):
            for ds_trim in range(args.excise_downstream_feature_trim[0], args.excise_downstream_feature_trim[1]+1):
                for us_hom in range(args.excise_upstream_homology[0], args.excise_upstream_homology[1]+1):
                    for ds_hom in range(args.excise_downstream_homology[0], args.excise_downstream_homology[1]+1):
                        for insert_length in range(args.excise_insert_lengths[0], args.excise_insert_lengths[1]+1):
                            if (args.excise_donor_lengths[0] <= us_hom+insert_length+ds_hom <= args.excise_donor_lengths[1]):
                                start1, end1 = f.start-us_hom-us_trim, f.start-us_trim
                                start2, end2 = f.end+ds_trim, f.end+ds_hom+ds_trim
                                upstream = contig_sequence[start1:end1]
                                downstream = contig_sequence[start2:end2]
                                #upstream = contigs[contig][start - args.excise_donor_homology[1]:start]
                                #downstream = contigs[contig][end:end + args.excise_donor_homology[1]]
                                
                                # when insert_length = 0, then the kmers are [''] (single element, empty string)
                                for mAT in nucleotides.kmers(insert_length):
                                    # Add this candidate dDNA to the list of all candidate dDNAs
                                    dDNA = upstream + mAT + downstream
                                    
                                    #if (args.excise_donor_lengths[0] <= len(dDNA) <= args.excise_donor_lengths[1]):
                                    #dDNAs.append(dDNA)
                                    new_targets = Target.get_targets(args, dDNA) # [(orientation, start, end, filt_seq, side, filt_spacer, filt_pam), ...]
                                    
                                    # The ReversionTarget must overlap the junction. These are valid cases:
                                    # dDNA    uuuuuuuuuuuuiiiidddddddddddd
                                    # valid      ----------
                                    # valid             ---------
                                    # valid                  ------
                                    # invalid   ----------
                                    # invalid                 -------
                                    # dDNA    uuuuuuuuuuuudddddddddddd
                                    # valid        --------
                                    # valid              --------
                                    
                                    for t in new_targets:
                                        # Target is a tuple: (orientation, start, end, upstream, downstream, sequence, side, spacer, pam, motif)
                                        # If it is not one of the two failure states, then add it (that is, it MUST overlap the junction)
                                        if not ((t[2] < len(upstream)) or (t[1] >= len(upstream) + len(mAT))):
                                            cls(f.name, f.contig, orientation, dDNA, (start1, end1), mAT, (start2, end2), spacer=t)
    
    @classmethod
    def bartag(cls, args, f, contig_sequence, orientation, bartags):
        """
        Method for generating dDNA with bartags.
        Does not check for ki-gRNA target sequences.
        
        This code is run once for each feature
        Does not respect:
          args.excise_upstream_homology
          args.excise_downstream_homology
          args.excise_upstream_feature_trim
          args.excise_downstream_feature_trim
          args.excise_insert_lengths
          args.excise_donor_lengths
        """
        
        # Use this instead of 'args.excise_donor_lengths'
        slen = 100
        
        # when insert_length = 0, then the kmers are [''] (single element, empty string)
        for bartag in bartags:
            # Add this bartag dDNA substring to the center of the up/down-stream DNA
            # dDNA = upstream[-slen//2-(-len(bartag)//2):] + bartag + downstream[:slen//2-len(bartag)//2]
            # start1=15, start2=5
            # dDNA = sequence[15+(-slen//2-(-len(bartag)//2)):15] + bartag + sequence[5:5+(slen//2-len(bartag)//2)]
            us_start, us_end = f.start+(-slen//2-(-len(bartag)//2)), f.start
            ds_start, ds_end = f.end, f.end+(slen//2-len(bartag)//2)
            upstream = contig_sequence[us_start:us_end]
            downstream = contig_sequence[ds_start:ds_end]
            dDNA = upstream + bartag + downstream
            
            # We don't calculate targets
            #new_targets = Target.get_targets(args, dDNA) # [(orientation, start, end, filt_seq, side, filt_spacer, filt_pam), ...]
            
            # We create the new ExcisionDonor object (with no spacers)
            cls(f.name, f.contig, orientation, dDNA, (us_start, us_end), bartag, (ds_start, ds_end), spacer=None)
    
    @classmethod
    def addtag(cls, args, f, contig_sequence, orientation):
        """
        Method for generating dDNAs with complete SPACER+PAM targeting sites
        with flanking homology arms
        """
        # Use this instead of 'args.excise_donor_lengths'
        slen = 100
        
        genome_composition = nucleotides.get_seq_dist(contig_sequence, 1, 4) # step_size=1, kmer_size=4
        
        for m in OnTargetMotif.motifs: # <------ Do I need to add OffTargetMotif.motifs here as well?
            addtag_seqs = [] # holds at most, 10x1000= 10000 sequences
            
            # Keep generating and testing addtags until the score raises above the arbitrary threshold of 0.9
            # Or until the maximum number of batches have been run
            max_batches = 10
            current_batch = 0
            best_score = 0
            while ((best_score < 0.9) and (current_batch < max_batches)):
                current_batch+= 1
                # Generate batches of 1000 addtags:
                addtag_batch = []
                for i in range(1000):
                    # Generate a pseudorandom motif
                    addtag_batch.append(m.generate_sequence(compositions=genome_composition, complement=True, samples=100))
                
                # Write the addtag sequences to a file, and align them to the genome
                #### This happens in 'main()' ####
                
                # Calculate the scores of the addtag sequences
                #### This happens in 'main()' ####
                
                # Add the batched oligos to the list of all generated ones
                addtag_seqs += addtag_batch
                
                best_score = 0
            
            # Only keep the top N best (i.e. 100)
            #for addtag in sorted(addtag_seqs, key=lambda x: x.score, reverse=True)[:100]:
            # Since we haven't scored these yet, then let's just keep all of them
            for addtag in addtag_seqs:
                # Get the contig locations of the homology regions
                us_start, us_end = f.start+(-slen//2-(-len(addtag)//2)), f.start
                ds_start, ds_end = f.end, f.end+(slen//2-len(addtag)//2)
                
                # Get the sequences for the homology regions
                upstream = contig_sequence[us_start:us_end]
                downstream = contig_sequence[ds_start:ds_end]
                
                # Stitch together the homology regions with the addtag (exogenous gRNA target)
                dDNA = upstream + addtag + downstream
                
                ##### Extra condition: ADDTAG must NOT be present in any of the ki-dDNAs #####
                
                # The target must match the 'addtag' sequence completely
                targets = Target.get_targets(args, dDNA) # [(orientation, start, end, filt_seq, side, filt_spacer, filt_pam), ...]
                
                for t in targets:
                    #                     0            1      2    3         4           5         6     7       8    9
                    # Target is a tuple: (orientation, start, end, upstream, downstream, sequence, side, spacer, pam, motif)
                    # The target must match the addtag start and end positions within the dDNA exactly
                    if ((t[1] == len(upstream)) and (t[2] == len(upstream)+len(addtag))):
                        # Add this candidate dDNA to the list of all candidate dDNAs
                        cls(f.name, f.contig, orientation, dDNA, (us_start, us_end), addtag, (ds_start, ds_end), spacer=t)
        #### End 'addtag()' ####
    
    @classmethod
    def unitag(cls, args, f, contig_sequence, orientation, unitag):
        """
        Method for generating dDNAs with complete SPACER+PAM targeting sites
        with flanking homology arms
        
        This function assumes the input 'unitag' has already been calculated
        """
        # Use this instead of 'args.excise_donor_lengths'
        slen = 100
        
        # Get the contig locations of the homology regions
        us_start, us_end = f.start+(-slen//2-(-len(unitag)//2)), f.start
        ds_start, ds_end = f.end, f.end+(slen//2-len(unitag)//2)
        
        # Get the sequences for the homology regions
        upstream = contig_sequence[us_start:us_end]
        downstream = contig_sequence[ds_start:ds_end]
        
        # Stitch together the homology regions with the addtag (exogenous gRNA target)
        dDNA = upstream + unitag + downstream
        
        ##### Extra condition: UNITAG must NOT be present in any of the ki-dDNAs #####
        
        # The target must match the 'addtag' sequence completely
        # This will search all OnTargetMotif motifs
        targets = Target.get_targets(args, dDNA) # [(orientation, start, end, filt_seq, side, filt_spacer, filt_pam), ...]
        
        for t in targets:
            #                     0            1      2    3         4           5         6     7       8    9
            # Target is a tuple: (orientation, start, end, upstream, downstream, sequence, side, spacer, pam, motif)
            # The target must match the addtag start and end positions within the dDNA exactly
            if ((t[1] == len(upstream)) and (t[2] == len(upstream)+len(unitag))):
                # Add this candidate dDNA to the list of all candidate dDNAs
                cls(f.name, f.contig, orientation, dDNA, (us_start, us_end), unitag, (ds_start, ds_end), spacer=t)
        #### End 'unitag()' ####
    
    @classmethod
    def generate_donors(cls, args, contigs):
        """
        Creates the DNA oligo with the structure:
        [upstream homology][unique gRNA][downstream homology]
        that excises the target feature
        """
        
        # Behave differently depending on user input (mintag/addtag/unitag/bartag)
        if (args.ko_dDNA == 'mintag'):
            # Generate the full set of potential dDNAs
            for feature_name, f in Feature.features.items():
            #for feature in features:
                #contig, start, end, strand = features[feature]
                # start & end are 0-based indices, inclusive/exclusive
                
                # assumes start < end
                # DNA 5' of feature is upstream
                # DNA 3' of feature is downstream
                
                contig_sequence = contigs[f.contig]
                orientation = '+' # ????? what is this for? I don't remember
                
                # For each potential dDNA, evaluate how good it is
                #cls.mintag_exhaustive_site_search(args, feature, contig, start, end, strand, contig_sequence, orientation)
                cls.mintag_exhaustive_site_search(args, f, contig_sequence, orientation)
            
        elif (args.ko_dDNA == 'addtag'): ######################## NEED TO ADD FLANKTAG PROCESSING ########################
            # Generate one unique tag for each feature (same take for homologs)
            # using the complement of the DNA sequence composition for the organism
            
            for feature_name, f in Feature.features.items():
                contig_sequence = contigs[f.contig]
                orientation = '+' # I forgot what this is for
                cls.addtag(args, f, contig_sequence, orientation)
            
        elif (args.ko_dDNA == 'unitag'): ######################## NEED TO ADD FLANKTAG PROCESSING ########################
            # Generate a single tag that is the same for all features
            
            # Generate a random tag (that is the complement sequence composition)
            # Evaluate all the tags
            # Pick the best tag
            # Add that best tag to all dDNAs
            
            # Pick an arbitrary motif (this is just the first-defined one)
            m = OnTargetMotif.motifs[0] # <------ Do I need to add OffTargetMotif.motifs here as well?
            
            # Create a single, random unitag according to the chosen motif (for testing)
            unitag = m.generate_sequence(compositions=genome_composition, complement=True, samples=100)
            
            # Create the dDNAs using the chosen unitag
            for feature_name, f in Feature.features.items():
                contig_sequence = contigs[f.contig]
                orientation = '+' # I forgot what this is for
                cls.unitag(args, f, contig_sequence, orientation, unitag)
        
        elif (args.ko_dDNA == 'bartag'): ######################## NEED TO ADD FLANKTAG PROCESSING ########################
            # if number of features x number of bartags > 200, then exit with warning
            if (args.bartag_number * len(Feature.features) > 200):
                raise Exception("Too many features or bartags. Must be <= 200.\nYou have: features × bartags = " + str(len(Feature.features)) + " × " + str(args.bartag_number) + " = " + str(args.bartag_number * len(Feature.features)))
            
            # Generate up to 200 barcodes
            logging.info('Calculating barcodes...')
            barcodes, barcode_fail_n, barcode_success_n = bartag.generate_min_distance(args.bartag_motif, args.bartag_distance, max_successes=200)
            logging.info('Found {} barcodes: {} fails, {} successes'.format(len(barcodes), barcode_fail_n, barcode_success_n))
            
            # Assign a barcode to each feature,
            # Generate the dDNAs
            bartag_assignments = {}
            for feature_name, f in Feature.features.items():
                # Add a barcodes equal to the number specified via the command line
                bartag_assignment[feature_name] = [barcodes.pop() for i in range(args.bartag_number)]
                
                contig_sequence = contigs[f.contig]
                orientation = '+' # I forgot what this is for
                cls.bartag(args, f, contig_sequence, orientation, bartag_assignment[feature_name])
        
        else: ######################## NEED TO ADD FLANKTAG PROCESSING ########################
            # 'args.ko_dDNA' is the path to a FASTA file
            logging.info('Processing ko-dDNA with FASTA as input.')
            
            # Load the FASTA file into memory
            ko_contigs = utils.old_load_fasta_file(args.ko_dDNA)
            
            if (len(ko_contigs) == 1):
                # If only a single sequence appears in FASTA file, then treat it as a unitag
                tag = list(ko_contigs.values())[0]
            
                for feature_name, f in Feature.features.items():
                    contig_sequence = contigs[f.contig]
                    orientation = '+' # I forgot what this is for
                    cls.unitag(args, f, contig_sequence, '+', tag)
                    cls.unitag(args, f, contig_sequence, '+', nucleotides.rc(tag))
            else:
                # If there are multiple sequences in FASTA file, then
                # cross-reference their primary sequence headers with
                # the (parent) feature names and gene names (from homologs file)
                for feature_name, f in Feature.features.items():
                    fp = f.get_expand_parent() # short for 'feature_parent'
                    try:
                        if fp.name in ko_contigs:
                            tag = ko_contigs[fp.name]
                        else:
                            tag = ko_contigs[feature2gene[fp.name]]
                    except KeyError:
                        raise Exception("The ko-dDNA file '{}' has no sequence with a primary header that matches '{}' or '{}'".format(args.ko_dDNA, fp.name, feature2gene[fp.name]))
                    
                    contig_sequence = contigs[f.contig]
                    orientation = '+' # I forgot what this is for
                    cls.unitag(args, f, contig_sequence, '+', tag) # Insert TAG in '+' orientation
                    cls.unitag(args, f, contig_sequence, '+', nucleotides.rc(tag)) # Insert TAG in '-' orientation
        
        #### End 'generate_donors()' ####
    
    @classmethod
    def generate_alignments(cls):
        # Eventually, this function should return something like the following:
        #                                                                  ┌insert
        #                 ┌──────────────────ds homology─────────────────┐┌┴┐┌──────────────────us homology──────────────────┐
        # exDonor-42 dDNA ACCATAACGTTTACTTGTTTAATATGCTATTGATATCTATATTTTTTTCCCTATGTGTAGTGCTTGTATATGCGTGTGTGATGAGAATAAGATGAATAGA
        # pam<spacer gRNA                                                 CCCTATGTGTAGTGCTTGTATAT
        for ind, obj in cls.indices.items():
            logging.info(obj.name)
            segment_string = ''
            loc = next(iter(obj.locations)) # Pull an arbitrary location record
            for segment in loc[3:]:
                if isinstance(segment, str):
                    length = len(segment)
                else:
                    length = segment[1] - segment[0]
                segment_string += cls.make_label(length, '')
            logging.info(segment_string)
            logging.info(obj.sequence)
            for s in obj.spacers:
                orientation = s[0]
                start = s[1]
                end = s[2]
                spacer = s[7]
                pam = s[8]
                if (orientation == '+'):
                    logging.info(' '*start + spacer + pam)
                else:
                    logging.info(' '*start + nucleotides.rc(spacer+pam))
    
    @staticmethod
    def make_label(length, label):
        """Helper function for 'generate_alignments()'"""
        if (length == 0):
            return ''
        if (length == 1):
            return '╥'
        if (length > 1):
            out = ['─'] * length
            out[0] = '┌'
            out[-1] = '┐'
            if (length >= len(label) + 4):
                start = int(length/2 - len(label)/2)
                for i, c in enumerate(label):
                    out[start+i] = c
            return ''.join(out)
    
    def get_inserts(self):
        """
        Returns a list constructed of the mAT inserts for locations
        """
        return [x[4] for x in self.locations]
    
    def get_trims(self):
        """
        Get the length of the upstream and downstream trims for each location
        """
        trims = []
        for loc in self.locations:
            #contig, start, end, strand = features[loc[0]]
            f = Feature.features[loc[0]]
            #mAT = loc[4]
            left_trim = f.start - loc[3][1]
            right_trim = loc[5][0] - f.end
            #trims.append(len(loc[4]) + loc[5][0]-l[3][1])
            trims.append((left_trim, right_trim))
        return trims
    
    def get_inserts_and_trims(self):
        """
        Returns list of insert size, and us/ds trims with the following format:
          [(insert, us-trim, ds-trim), ...]
        """
        return_list = []
        for loc in self.locations:
            f = Feature.features[loc[0]]
            #contig, start, end, strand = features[loc[0]]
            mAT = loc[4]
            left_trim = f.start - loc[3][1]
            right_trim = loc[5][0] - f.end
            return_list.append((mAT, left_trim, right_trim))
        return return_list

class ReversionDonor(Donor):
    prefix = 'reDonor'
    sequences = {}
    indices = {}
    
    @classmethod
    def generate_donors(cls, args, contigs):
        """
        Creates the DNA oligo with the structure:
        [amp-F primer]                              [amp-R primer]
        [upstream homology][original feature][downstream homology]
        """
        tm_max_difference = 4.0
        amplicon_size = (min(args.revert_upstream_homology)+min(args.revert_downstream_homology), max(args.revert_upstream_homology)+max(args.revert_downstream_homology))
        tm_range = (52, 64)
        primer_length_range = (19, 32)
        min_delta_g = -5.0
        
        subset_size = 100 # 500 # Temporarily limit number of Primer objects used as input to the pair() function
        pp_subset_size = 500 # Temporarily limit the number of PrimerPairs that are used to create ReversionDonor objects
        temp_folder = os.path.join('/dev/shm/addtag', os.path.basename(args.folder))
        
        # Make a folder to store the *.gb Genbank flat files
        os.makedirs(os.path.join(args.folder, 'reversion-gb'), exist_ok=True)
        
        # There should be at least one ReversionDonor for each feature
        for feature_name, f in Feature.features.items():
            my_contig = contigs[f.contig]
            
            #feature_length = f.end - f.start
            feature_sequence = my_contig[f.start:f.end]
            
            orientation = '+'
            
            if (f.strand == '-'):
                region_F_start, region_F_stop = f.start-max(args.revert_downstream_homology), f.start-min(args.revert_downstream_homology)
                region_R_start, region_R_stop = f.end+min(args.revert_upstream_homology), f.end+max(args.revert_upstream_homology)
            else:
                region_F_start, region_F_stop = f.start-max(args.revert_upstream_homology), f.start-min(args.revert_upstream_homology)
                region_R_start, region_R_stop = f.end+min(args.revert_downstream_homology), f.end+max(args.revert_downstream_homology)
            
            region_F = my_contig[region_F_start:region_F_stop]
            region_R = my_contig[region_R_start:region_R_stop]
            
            logging.info('Calculating upstream_F primers...')
            upstream_F = args.selected_oligo.scan(region_F,   'left',  primer_size=primer_length_range, tm_range=tm_range, min_delta_g=min_delta_g, folder=temp_folder)
            upstream_F = sorted(upstream_F, key=lambda x: x.weight, reverse=True)
            logging.info('  len(upstream_F) = {}'.format(len(upstream_F)))
            
            logging.info('Calculating downstream_R primers...')
            downstream_R = args.selected_oligo.scan(region_R, 'right', primer_size=primer_length_range, tm_range=tm_range, min_delta_g=min_delta_g, folder=temp_folder)
            downstream_R = sorted(downstream_R, key=lambda x: x.weight, reverse=True)
            logging.info('  len(downstream_R) = {}'.format(len(downstream_R)))
            
            upstream_F = upstream_F[:subset_size]
            downstream_R = downstream_R[:subset_size]
            
            if ((subset_size != None) and (len(upstream_F) > subset_size)):
                logging.info('upstream_F: skipping {}/{} calculated primers'.format(max(0, len(upstream_F)-subset_size), len(upstream_F)))
            if ((subset_size != None) and (len(downstream_R) > subset_size)):
                logging.info('downstream_R: skipping {}/{} calculated primers'.format(max(0, len(downstream_R)-subset_size), len(downstream_R)))
            
            logging.info('Calculating: primer pairs...')
            uf_dr_paired_primers = args.selected_oligo.pair(upstream_F, downstream_R, amplicon_size=amplicon_size, tm_max_difference=tm_max_difference, intervening=0, min_delta_g=min_delta_g, folder=temp_folder)
            uf_dr_paired_primers = sorted(uf_dr_paired_primers, key=lambda x: x.get_joint_weight(), reverse=True)
            logging.info('  len(uf_dr_paired_primers) = {}'.format(len(uf_dr_paired_primers)))
            
            # Ammend the intervening sequence
            for pp in uf_dr_paired_primers:
                pp.intervening = region_R_start - region_F_stop
            
            if ((pp_subset_size != None) and (len(uf_dr_paired_primers) > pp_subset_size)):
                logging.info('Skipping {}/{} primer pairs'.format(max(0, len(uf_dr_paired_primers)-pp_subset_size), len(uf_dr_paired_primers)))
            
            
            ###### alignment ######
            pp_labels_list = []
            ###### alignment ######
            
            
            logging.info('\t'.join(['ReversionDonor', 'feature', 'weight', 'PrimerPair']))
            for pp_count, pp in enumerate(uf_dr_paired_primers[:pp_subset_size]):
                start1, end1 = region_F_start + pp.forward_primer.position, f.start
                start2, end2 = f.end, region_R_start+pp.reverse_primer.position+len(pp.reverse_primer.sequence)
                upstream_seq = my_contig[start1:end1]
                downstream_seq = my_contig[start2:end2]
                
                dDNA = upstream_seq + feature_sequence + downstream_seq
                
                cls(feature_name, f.contig, orientation, dDNA, (start1, end2), spacer=None)
                
                rd = cls.sequences[dDNA]
                
                logging.info('\t'.join([rd.name, feature_name, str(pp.get_joint_weight()), str(pp)]))
                
                ###### alignment ######
                pp_labels_list.append(rd.name)
                ###### alignment ######
                
                ############# Write to Genbank flat file #############
                colors = {
                    'grey': '#DDDDDD',
                    'gray': '#DDDDDD',
                    'pink': '#FF9CCD',
                    'bpink': '#FF9C9A', # salmon
                    'dpink': '#DDB4DD', # grey-pink
                }
                segment_start, segment_end = max(0, f.start-2000), min(len(my_contig), f.end+2000)
                genome_segment = my_contig[segment_start:segment_end]
                gb = utils.GenBankFile(f.contig + ':' + str(segment_start+1) + '-' + str(segment_end), genome_segment) # chr:start-end
                #gb.add_annotation('source', 0, len(genome_segment), '+', organism=args.fasta[0], mol_type='genomic DNA') # only uses the first listed args.fasta
                gb.add_annotation('feature', f.start-segment_start, f.end-segment_start, f.strand, label='feature', ApEinfo_revcolor=colors['grey'], ApEinfo_fwdcolor=colors['grey'])
                gb.add_annotation('region', region_F_start-segment_start, region_F_stop-segment_start, '+', label='upstream_primer_region', ApEinfo_revcolor=colors['dpink'], ApEinfo_fwdcolor=colors['dpink'])
                gb.add_annotation('region', region_R_start-segment_start, region_R_stop-segment_start, '-', label='downstream_primer_region', ApEinfo_revcolor=colors['dpink'], ApEinfo_fwdcolor=colors['dpink'])
                gb.add_annotation('amplicon', start1-segment_start, end2-segment_start, '+', label='ampf_ampr_amplicon', ApEinfo_revcolor=colors['bpink'], ApEinfo_fwdcolor=colors['bpink']) # salmon
                gb.add_annotation('primer', start1 - segment_start, start1 - segment_start + len(pp.forward_primer.sequence), '+', label='ampf_primer', ApEinfo_revcolor=colors['pink'], ApEinfo_fwdcolor=colors['pink']) # pink
                gb.add_annotation('primer', end2 - segment_start-len(pp.reverse_primer.sequence), end2 - segment_start, '-', label='ampr_primer', ApEinfo_revcolor=colors['pink'], ApEinfo_fwdcolor=colors['pink']) # pink
                gb.write(os.path.join(args.folder, 'reversion-gb', 'pp-'+str(pp_count)+'.gb')) # pp_labels_list[-1]
                ############# Write to Genbank flat file #############
                
            
            ###### alignment ######
            label_list = ['us_region', 'us_skipped', 'feature', 'ds_skipped', 'ds_region']
            sequence_list = [region_F, my_contig[region_F_stop:f.start], feature_sequence, my_contig[f.end:region_R_start], region_R]
            aln_out = make_labeled_primer_alignments(label_list, sequence_list, f.contig, pp_labels_list, uf_dr_paired_primers[:pp_subset_size])
            
            for oline in aln_out:
                print(oline)
            ###### alignment ######
            
            
            
            
            
            
    
    @classmethod
    def generate_donors_old(cls, args, contigs):
        """
        Creates the DNA oligo with the structure:
        [upstream homology][original feature][downstream homology]
        """
        # There should be one ReversionDonor for each feature
        for feature_name, f in Feature.features.items():
        #for feature in features:
            # First, get the homology blocks up- and down-stream of the feature
            #contig, start, end, strand = features[feature]
            
            feature_length = f.end - f.start
            
            my_contig = contigs[f.contig]
            orientation = '+'
            
            for us_hom in range(args.revert_upstream_homology[0], args.revert_upstream_homology[1]+1):
                for ds_hom in range(args.revert_downstream_homology[0], args.revert_downstream_homology[1]+1):
                    if (args.revert_donor_lengths[0] <= us_hom+feature_length+ds_hom <= args.revert_donor_lengths[1]):
                        # to implement:
                        #   make revert_upstream_homology and revert_downstream_homology exclude the trim sequences
                        if (f.strand == '+'):
                            start1, end1 = f.start-us_hom, f.start
                            start2, end2 = f.end, f.end+ds_hom
                        else:
                            start1, end1 = f.start-ds_hom, f.start
                            start2, end2 = f.end, f.end+us_hom
                        upstream = my_contig[start1:end1]
                        downstream = my_contig[start2:end2]
                        feature_sequence = my_contig[f.start:f.end]
                        
                        dDNA = upstream + feature_sequence + downstream
                        #for re_seq, re_target in ReversionTarget.indices.items():
                        #    if feature in [x[0] for x in re_target.locations]:
                        #        cls(feature, contig, orientation, (start1, end1), '..'.join(map(str, [start, end])), (start2, end2), dDNA, re_target)
                        #cls(feature, contig, orientation, (start1, end1), '..'.join(map(str, [start, end])), (start2, end2), dDNA, None)
                        cls(feature_name, f.contig, orientation, dDNA, (start1, end2), spacer=None)

class Target(object):
    """Data structure defining a gRNA Target"""
    prefix = 'Target'
    sequences = {} # key = nucleotide sequence, value = ExcisionTarget/ReversionTarget object
    indices = {} # key = exTarget-102, value = ExcisionTarget/ReversionTarget object
    
    @classmethod
    def load_alignment(cls, filename, args, contigs):
        """
        Read in alinment file (SAM/BLASTN).
        sep is the separator for the header. Positions are converted to 0-index
        Creates a list of Sequence objects
        """
        
        # Code to decompress a *.bam file should go here
        with open(filename, 'r') as flo:
            args.selected_aligner.current_file = filename
            record = None
            while True:
                record = args.selected_aligner.load_record(flo)
                #logging.info(record)
                if (record == None):
                    break
                else:
                    # Record(
                    #     query_name, subject_name,
                    #     query_sequence, subject_sequence,
                    #     query_position, subject_position,
                    #     query_length, subject_length,
                    #     flags, cigar, score, evalue, length
                    # )
                    
                    target = cls.indices[record.query_name] # key=exTarget-519
                    
                    # Get orientation
                    alignment_orientation = utils.sam_orientation(record.flags)
                    
                    # Get alignment position
                    alignment_contig = record.subject_name
                    alignment_start, alignment_end = record.subject_position
                    
                    # Reverse-complement if needed
                    alignment_contig_sequence = contigs[alignment_contig]
                    alignment_sequence = alignment_contig_sequence[alignment_start:alignment_end]
                    alignment_upstream = alignment_contig_sequence[alignment_start-10:alignment_start]
                    alignment_downstream = alignment_contig_sequence[alignment_end:alignment_end+10]
                    #actual_sequence = record.query_sequence # Should be column: 9, SEQ, query SEQuence on the same strand as the reference
                    if (alignment_orientation == '-'):
                        alignment_sequence = nucleotides.rc(alignment_sequence)
                        alignment_upstream, alignment_downstream = nucleotides.rc(alignment_downstream), nucleotides.rc(alignment_upstream)
                        #actual_sequence = nucleotides.rc(actual_sequence)
                    
                    # No actual check made here to see if length of query matches length of subject, and they still conform to motif.
                    # This is done inside the 'target.add_alignment()' method
                    
                    target.add_alignment(
                        args,
                        alignment_sequence, # aligned_sequence (as when matched with reference, thus may be revcomp of initial query)
                        alignment_contig, # aligned_contig
                        alignment_start, # aligned_start
                        alignment_end, # aligned_end
                        alignment_orientation, # aligned_orientation (+/-)
                        alignment_upstream,
                        alignment_downstream,
                    )
        
        logging.info(cls.__name__ + ' alignment file parsed: {!r}'.format(filename))
    
    @classmethod
    def generate_query_fasta(cls, filename, sep=':'):
        """
        Creates a FASTA file (That becomes header for SAM file)
        >exTarget-id feature:contig:orientation:start..end ...
        """
        with open(filename, 'w') as flo:
            for sequence, obj in sorted(cls.sequences.items(), key=lambda x: int(x[1].name.split('-')[1])):
                #for location in obj.locations:
                #    print(' '.join(['>'+obj.format_location(location), obj.name]), file=flo)
                #    print(sequence, file=flo)
                print(' '.join(['>'+obj.name] + sorted([obj.format_location(x, sep) for x in obj.locations])), file=flo)
                print(sequence, file=flo)
                
        logging.info(cls.__name__ + ' query FASTA generated: {!r}'.format(filename))
        return filename
    
    @classmethod
    def generate_spacers_fasta(cls, filename, sep=':'):
        """
        Creates a FASTA file
        >exTarget-id motif=STR locations=N alignments=N/N on-target=N off-target=N
        """
        with open(filename, 'w') as flo:
            for sequence, obj in cls.sequences.items():
                print(' '.join([
                    '>' + obj.name,
                    'motif=' + obj.motif,
                    'locations=' + str(len(obj.locations)),
                    'alignments=' + str(len([a for a in obj.alignments if a.postfilter])) + '/' + str(len(obj.alignments)),
                    'on-target=' + str(round(obj.score['Azimuth'], 2)),
                    'off-target=' + str(round(obj.off_targets['Hsu-Zhang'], 2)),
                    'pam=' + obj.pam
                ]), file=flo)
                print(obj.spacer, file=flo)
        logging.info(cls.__name__ + ' spacers FASTA generated: {!r}'.format(filename))
        return filename
    
    def format_location(self, location, sep=':'):
        feature, contig, orientation, start, end, upstream, downstream = location
        return sep.join([feature, contig, orientation, '..'.join([str(start), str(end)])])
    
    def format_sequence(self):
        """
        Returns nucleotides with proper 'SPACER>PAM' or 'PAM<SPACER' formatting.
        """
        if (self.side == '>'):
            return self.spacer + '>' + self.pam
        else:
            return self.pam + '<' + self.spacer
    
    def __init__(self, feature, contig, orientation, start, end, upstream, downstream, sequence, side, spacer, pam, motif, parsed_motif, add=True):
        """Create a structure for holding individual sequence information"""
        location = (feature, contig, orientation, start, end, upstream, downstream)
        self.locations = set()
        
        self.sequence = sequence
        self.side = side
        self.spacer = spacer
        self.pam = pam
        self.motif = motif
        self.parsed_motif = parsed_motif
        
        # List to store alignments
        self.alignments = []
        
        # Get the index number
        self.index = len(self.sequences)
        self.name = self.prefix + '-' + str(self.index)
        
        if add:
            the_target = self.sequences.setdefault(self.sequence, self)
            the_target.locations.add(location)
            self.indices.setdefault(the_target.name, the_target)
        
        # Scores for this sequence only (not PairedSequenceAlgorithm)
        self.score = {}
        
        # We need an off-target score for each algorithm and each potential dDNA
        # Thus it will be structured as a dict nested within a dict:
        # self.off_targets = {'CFD': {'exDonor-0': 10.0, 'exDonor-1': 12.0}}
        self.off_targets = {}
        if (len(self.locations) > 0):
            self.calculate_default_scores()
    
    @classmethod
    def get_targets(cls, args, sequence):
        """
        Tries to match all OnTargetMotif motifs to the input sequence
        ANYWHERE in the sequence
        Returns list of matches (as a tuple) from all OnTargetMotif motifs.
        """
        targets = set()
        #for seq_i, sequence in enumerate(dDNAs):
        for orientation in ['+', '-']:
            if (orientation == '-'):
                sequence = nucleotides.rc(sequence)
            
            #for i in range(len(args.parsed_motifs)):
            for mymotif in OnTargetMotif.motifs: # <------ Do I need to add OffTargetMotif.motifs here as well?
                #spacers, pams, side = args.parsed_motifs[i]
                spacers, pams, side = mymotif.parsed_list
                #compiled_regex = args.compiled_motifs[i]
                #matches = nucleotides.motif_search(sequence, spacers, pams, side)
                matches = nucleotides.motif_search2(sequence, side, mymotif.compiled_regex)
                for seq, start, end, spacer, pam in matches:
                    if (orientation == '-'):
                        start, end = len(sequence) - end, len(sequence) - start
                    upstream = sequence[start-10:start]
                    downstream = sequence[end:end+10]
                    filtered_targets = target_filter(seq, spacer, pam, upstream, downstream, args)
                    for filt_seq, filt_spacer, filt_pam in filtered_targets:
                        targets.add((orientation, start, end, upstream, downstream, filt_seq, side, filt_spacer, filt_pam, mymotif.motif_string, tuple([tuple(x) if isinstance(x, list) else x for x in mymotif.parsed_list])))
        return sorted(targets) # becomes a list
    
    def calculate_default_scores(self):
        """Populate this scores for this Sequence"""
        loc = next(iter(self.locations)) # Pull an arbitrary location record
        #parent = (self.sequence, self.target, self.pam, self.upstream, self.downstream)
        parent = (self.sequence, self.spacer, self.pam, loc[5], loc[6])
        for C in algorithms.single_algorithms:
            if (C.default != None):
                self.score[C.name] = C.default
            else:
                self.score[C.name] = C.calculate(parent, parsed_motif=self.parsed_motif)
        for C in algorithms.paired_algorithms:
            if (C.default != None):
                self.score[C.name] = C.default
            else:
                self.score[C.name] = 0.0
    
    def add_alignment(self, args, aligned_sequence, aligned_contig, aligned_start, aligned_end, aligned_orientation, aligned_upstream, aligned_downstream):
        """Add a genomic position to the list of alignments"""
        loc = next(iter(self.locations)) # Pull an arbitrary location record
        #parent = (self.sequence, self.target, self.pam, self.upstream, self.downstream)
        parent = (self.sequence, self.spacer, self.pam, loc[5], loc[6])
        aligned_target, aligned_pam, aligned_motif = self.split_spacer_pam(aligned_sequence, args)
        #aligned_parsed_motif = args.parsed_motifs[args.motifs.index(aligned_motif)]
        #aligned_parsed_motif = aligned_motif.parsed_list
        if ((aligned_target != None) and (aligned_pam != None)):
            a = Alignment(
                aligned_sequence,
                aligned_target,
                aligned_pam,
                aligned_motif,
                aligned_contig,
                aligned_start,
                aligned_end,
                aligned_orientation,
                upstream=aligned_upstream,
                downstream=aligned_downstream
            )
            a.calculate_scores(parent)
            self.alignments.append(a)
        else:
            print('Cannot add alignment:', aligned_sequence, aligned_contig, aligned_start, aligned_end, aligned_orientation, file=sys.stderr)
    
    def split_spacer_pam(self, sequence, args):
        #for i in range(len(args.parsed_motifs)):
        for motif in OnTargetMotif.motifs:
            #spacers, pams, side = args.parsed_motifs[i]
            spacers, pams, side = motif.parsed_list
            #compiled_regex = args.compiled_motifs[i]
            m = nucleotides.motif_conformation2(sequence, side, motif.compiled_regex)
            
            # If the motif matches, then return it immediately
            if m:
                return m[0], m[1], motif # spacer, pam, motif
        
        for motif in OffTargetMotif.motifs:
            spacers, pams, side = motif.parsed_list
            m = nucleotides.motif_conformation2(sequence, side, motif.compiled_regex)
            if m:
                return m[0], m[1], motif
    
        # If the motif does not match, then fudge it for an off-target/on-target motif
        # This will return the last-fudged motif, and not necessarily the best-fudged motif
        # *MAY* need to improve this code later
        for motif in (OnTargetMotif.motifs + OffTargetMotif.motifs)[::-1]:
            spacers, pams, side = motif.parsed_list
            if (side == '>'):
                l = max(map(len, pams))
                return sequence[:-l], sequence[-l:], motif # spacer, pam, motif
            else: #elif (side == '<'):
                l = max(map(len, pams))
                return seq[:l], seq[l:], motif # spacer, pam, motif
    
    def get_features(self):
        """Return (sorted) list of all feature names this Target maps to"""
        return sorted(set(x[0] for x in self.locations))
    
    def get_location_features(self):
        """Return (sorted) list of all feature names this Target maps to"""
        return sorted(set(x[0] for x in self.locations))
    
    def get_contigs(self):
        """Return (sorted) list of all contig names this Target maps to"""
        return sorted(set(x[1] for x in self.locations))
    
    def get_parent(self):
        """Returns a parent tuple for an arbitrary location record"""
        loc = next(iter(self.locations)) # Pull an arbitrary location record
        #parent = (self.sequence, self.target, self.pam, self.upstream, self.downstream)
        parent = (self.sequence, self.spacer, self.pam, loc[5], loc[6])
        return parent
    
    @classmethod
    def score_batch(cls):
        """Performs BatchedSingleSequenceAlgorithm calculations on all sequences"""
        # Get immutable list of dict keys
        t_sorted = sorted(cls.indices.keys(), key=lambda x: int(x.split('-')[1]))
        
        # Make list to populate with values to pass into the calculate() methods
        queries = []
        for i, t_index in enumerate(t_sorted):
            t_obj = cls.indices[t_index]
            loc = next(iter(t_obj.locations)) # Pull an arbitrary location record
            parent = (t_obj.sequence, t_obj.spacer, t_obj.pam, loc[5], loc[6]) # upstream, downstream
            queries.append(parent)
        
        # Loop through the Algorithms
        for C in algorithms.batched_single_algorithms:
            # Calculate
            batch_scores = C.calculate(queries) # motif=cls.motif
            
            # Assign the score to the appropriate Target
            for i, t_index in enumerate(t_sorted):
                t_obj = cls.indices[t_index]
                t_obj.score[C.name] = batch_scores[i]
    
    def score_off_targets(self, args, homologs):
        """
        Calculate Guide Score (off-target score) for all single/paired
        algorithms with the 'off_target=True' attribute.
        """
        # Get list of algorithms whose scores should be used for off-target
        # calculations
        calculators = []
        for C in algorithms.single_algorithms:
            if C.off_target:
                calculators.append(C)
        for C in algorithms.paired_algorithms:
            if C.off_target:
                calculators.append(C)
        
        # Make empty lists for scores
        # It should have this format
        #  dict[algorithm][genome/exDonor/reDonor] = [score1, score2, score3, ...]
        # This is for keeping track of which scores should get the ratio multiplier
        on_targets = {}
        off_targets = {}
        for C in calculators:
            on_targets[C.name] = {'gDNA': [], 'dDNA': []}
            off_targets[C.name] = {'gDNA': [], 'dDNA': []}
        
        # Get list of on-target features
        on_target_features = set([x[0] for x in self.locations])
        if isinstance(self, ReversionTarget):
            on_target_features = set(utils.flatten([x.split(',') for x in on_target_features]))
        
        if homologs:
            for f in list(on_target_features):
                on_target_features.update(homologs.get(f, set()))
        
        logging.info('on_target_features' + str(on_target_features))
        # Check each alignment
        for a in self.alignments:
            a_features = None
            temp_string = ''
            if a.postfilter:
                if isinstance(self, ExcisionTarget):
                    if a.contig.startswith('exDonor-'):
                        a_features = ExcisionDonor.indices[a.contig].get_features()
                        temp_string += 'case exT-exD'.rjust(18) + ' '
                        
                        for i, C in enumerate(calculators):
                            c_score = a.score[C.name]
                            if (C.minimum <= c_score <= C.maximum):
                                if ((len(a_features) == 0) or (len(on_target_features.intersection(a_features)) == 0)):
                                    # if alignment is NOT a target:
                                    # And is also an exDonor, then we will ignore it
                                    # As we are assuming mono-plex gRNA
                                    a.action = 'skip'
                                # If the features match, then this is an off-target
                                else:
                                    # Scale this by the ratio, as they are dDNA alignments
                                    off_targets[C.name]['dDNA'].append(c_score)
                                    a.action = 'off'
                                    
                    else: # This is a genomic contig and not a dDNA
                        a_features = Feature.get_overlapping_features(a.contig, a.start, a.end)
                        #a_features = self.get_features(features, a.contig, a.start, a.end) # contig, start, end
                        temp_string += 'case exT-gDNA'.rjust(18) + ' '
                        
                        for i, C in enumerate(calculators):
                            c_score = a.score[C.name]
                            if (C.minimum <= c_score <= C.maximum):
                                # if alignment is NOT a target:
                                if ((len(a_features) == 0) or (len(on_target_features.intersection(a_features)) == 0)):
                                    off_targets[C.name]['gDNA'].append(c_score)
                                    a.action = 'off'
                                #if (len(on_target_features.intersection(a_features)) > 0):
                                # if alignment is an intended target
                                else:
                                    on_targets[C.name]['gDNA'].append(c_score)
                                    a.action = 'on'
                
                elif isinstance(self, ReversionTarget):
                    if a.contig.startswith('exDonor-'):
                        # reTarget vs exDonor:
                        #   if multiplex:
                        #     if exDonor == intended exDonor (for this feature), then on-target
                        #     if exDonor is the same feature, but NOT intended, then ignore
                        #     if exDonor is for another feature, then it is off-target
                        #   if monoplex:
                        #     if exDonor == intended exDonor (for this feature), then on-target
                        #     if exDonor is the same feature, but NOT intended, then ignore
                        #     if exDonor is for another feature, then ignore
                        #a_features = ExcisionDonor.indices[a.contig].get_features()
                        intended_contigs = self.get_contigs()
                        temp_string += 'case reT-exD'.rjust(18) + ' '
                        
                        for i, C in enumerate(calculators):
                            c_score = a.score[C.name]
                            if (C.minimum <= c_score <= C.maximum):
                                if a.contig in intended_contigs:
                                    on_targets[C.name]['gDNA'].append(c_score)
                                    a.action = 'on'
                                else:
                                    a.action = 'skip'
                                
                                #if ((len(a_features) == 0) or (len(on_target_features.intersection(a_features)) == 0)):
                                #    pass # do nothing for monoplex
                                #else:
                                #    on_targets[C.name]['gDNA'].append(c_score)
                    
                    elif a.contig.startswith('reDonor-'):
                        # reTarget vs reDonor:
                        #   multiplex: these are all off-target (use ratio)
                        #   monoplex: if reDonor has same feature as reTarget, then it is an off-target (use ratio)
                        a_features = ReversionDonor.indices[a.contig].get_features()
                        temp_string += 'case reT-reD'.rjust(18) + ' '
                        
                        for i, C in enumerate(calculators):
                            c_score = a.score[C.name]
                            if (C.minimum <= c_score <= C.maximum):
                                if ((len(a_features) == 0) or (len(on_target_features.intersection(a_features)) == 0)):
                                    a.action = 'skip'
                                else:
                                    # use ratio
                                    off_targets[C.name]['dDNA'].append(c_score)
                                    a.action = 'off'
                    else:
                        # reTarget vs genome: these are all off-target, regardless of the feature
                        temp_string += 'case reT-gDNA'.rjust(18) + ' '
                    
                        for i, C in enumerate(calculators):
                            c_score = a.score[C.name]
                            if (C.minimum <= c_score <= C.maximum):
                                off_targets[C.name]['gDNA'].append(c_score)
                                a.action = 'off'
            else:
                temp_string += 'failed post-filter' + ' '
            temp_string += ' '.join([a.action.rjust(4), str(a), str(a_features)])
            logging.info(temp_string)
        
        for i, C in enumerate(calculators):
            on_str = str(len(on_targets[C.name]['gDNA'])) + ' + (' + str(args.dDNA_gDNA_ratio)+')'+str(len(on_targets[C.name]['dDNA']))
            off_str = str(len(off_targets[C.name]['gDNA'])) + ' + (' + str(args.dDNA_gDNA_ratio)+')'+str(len(off_targets[C.name]['dDNA']))
            logging.info(' '.join([self.name, C.name, '('+on_str+')/(('+on_str+') + ('+off_str+'))']))
        
        # Perform off-target calculations
        for C in calculators:
            try:
                on_list = on_targets[C.name]['gDNA'] + args.dDNA_gDNA_ratio * on_targets[C.name]['dDNA']
                off_list = off_targets[C.name]['gDNA'] + args.dDNA_gDNA_ratio * off_targets[C.name]['dDNA']
                self.off_targets[C.name] = scores.off_target_score(off_list, on_list)
            except ZeroDivisionError:
                self.off_targets[C.name] = 0.0 # This should never happen
    
    def __repr__(self):
        """Return a string containing a printable representation of the Target object."""
        return self.__class__.__name__ + '(' + ' '.join([
            self.name,
            self.format_sequence(),
            'motif=' + self.motif,
            'locations=' + str(len(self.locations)),
            'alignments=' + str(len([a for a in self.alignments if a.postfilter])) + '/' + str(len(self.alignments))] +
            [x + '=' + str(round(self.score[x], 2)) for x in self.score] + 
            ['OT:' + x + '=' + str(round(self.off_targets[x], 2)) for x in self.off_targets]
            ) + ')'

class ExcisionTarget(Target):
    prefix = 'exTarget'
    # Non-redundant dicts
    sequences = {} # key = nucleotide sequence, value = ExcisionTarget object
    indices = {} # key = exTarget-102, value = ExcisionTarget object
    
    @classmethod
    def old_search_all_features(cls, args, features, contigs):
        for feature_name in features:
            feature_contig, feature_start, feature_end, feature_strand = features[feature_name]
            
            if feature_contig in contigs:
                cls.expand_feature(args, (feature_name, feature_contig, feature_start, feature_end, feature_strand), (feature_contig, contigs[feature_contig]))
            else:
                logging.info("The contig '{}' for feature '{}' is not in the input FASTA.".format(feature_contig, feature_name))
        
    
    @classmethod
    def search_all_features(cls, args, contigs):
        for feature_name, f in Feature.features.items():
            contig_sequence = contigs[f.contig]
            ef = f.expand_feature(args, contig_sequence) # Need to improve this: expanded feature should have a "parent" feature
            
            # Search for targets in the feature
            feature_sequence = contig_sequence[ef.start:ef.end]
            targets = Target.get_targets(args, feature_sequence) # Does both orientations (+/-)
            
            # Create ExcisionTarget objects for each found target
            for t in targets:
                t_orientation = t[0]
                t_start = t[1]
                t_end = t[2]
                t_upstream = t[3]
                t_downstream = t[4]
                t_sequence = t[5]
                t_side = t[6]
                t_spacer = t[7]
                t_pam = t[8]
                t_motif_string = t[9]
                t_motif_parsed_list = t[10]
                
                if (t_orientation == '+'):
                    real_start = f.start + t_start
                    real_end = f.start + t_end
                else:
                    contig_length = len(contig_sequence)
                    real_start = contig_length - t_end + f.start
                    real_end = contig_length - t_start + f.start
                
                #OLD cls(feature, feature_contig, orientation, real_start, real_end, t_upstream, t_downstream, filt_seq, side, filt_spacer, filt_pam, motif.motif_string, motif.parsed_list)
                #REFERENCE CODE targets.add((orientation, start, end, upstream, downstream, filt_seq, side, filt_spacer, filt_pam, mymotif.motif_string, tuple([tuple(x) if isinstance(x, list) else x for x in mymotif.parsed_list])))
                cls(f.name, f.contig, t_orientation, real_start, real_end, t_upstream, t_downstream, t_sequence, t_side, t_spacer, t_pam, t_motif_string, t_motif_parsed_list)
                # t = (orientation, start, end, t_upstream, t_downstream, filt_seq, side, filt_spacer, filt_pam, args.motifs[i], args.parsed_motifs[i])
                #SIMILAR ReversionTarget(obj_features, obj.name, t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], t[10])
            
            
            
    
    @classmethod
    def old_expand_feature(cls, args, f, contig_sequence, expansion_size=100, minimum_targets_per_feature=10):
        
        feature_start = f.start
        feature_end = f.end
        
        contig_length = len(contig_sequence)
        
        if (feature_end == None):
            feature_end = len(contig_sequence)
        
        targets = []
        while ((len(targets) < minimum_targets_per_feature) and ((feature_start, feature_end) != (0, contig_length))):
            feature_sequence = contig_sequence[feature_start:feature_end]
            targets = Target.get_targets(args, feature_sequence) # Does both orientations (+/-)
            
            feature_start = max(0, feature_start-expansion_size)
            feature_end = min(contig_length, feature_end + expansion_size)
        
        # Create ExcisionTarget objects for each found target
        for t in targets:
            t_orientation = t[0]
            t_start = t[1]
            t_end = t[2]
            t_upstream = t[3]
            t_downstream = t[4]
            t_sequence = t[5]
            t_side = t[6]
            t_spacer = t[7]
            t_pam = t[8]
            t_motif_string = t[9]
            t_motif_parsed_list = t[10]
            
            if (t_orientation == '+'):
                real_start = feature_start + t_start
                real_end = feature_start + t_end
            else:
                real_start = contig_length - t_end + feature_start
                real_end = contig_length - t_start + feature_start
            
            #OLD cls(feature, feature_contig, orientation, real_start, real_end, t_upstream, t_downstream, filt_seq, side, filt_spacer, filt_pam, motif.motif_string, motif.parsed_list)
            #REFERENCE CODE targets.add((orientation, start, end, upstream, downstream, filt_seq, side, filt_spacer, filt_pam, mymotif.motif_string, tuple([tuple(x) if isinstance(x, list) else x for x in mymotif.parsed_list])))
            cls(f.name, f.contig, t_orientation, real_start, real_end, t_upstream, t_downstream, t_sequence, t_side, t_spacer, t_pam, t_motif_string, t_motif_parsed_list)
            # t = (orientation, start, end, t_upstream, t_downstream, filt_seq, side, filt_spacer, filt_pam, args.motifs[i], args.parsed_motifs[i])
            #SIMILAR ReversionTarget(obj_features, obj.name, t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], t[10])
    
    # @classmethod
    # def get_targets(cls, args, contigs, features):
    #     """
    #     Searches within the annotated features on all contigs for targets
    #     that match the args.motifs criteria, then filters them.
    #     
    #     The constructor call adds the valid gRNA sites to 'ExcisionTarget.sequences' and 'ExcisionTarget.indices'.
    #     """
    #     # Find unique gRNA sites within each feature
    #     # Use a sliding window to make a list of queries
    #     for feature in features: # This is the filtered list of features, so we want to find all of these
    #         
    #         feature_contig, feature_start, feature_end, feature_strand = features[feature]
    #         if (feature_end == None):
    #             feature_end = len(contigs[feature_contig])
    #         
    #         # Make sure the contig the feature is on is present in the FASTA
    #         if feature_contig in contigs:
    #             # Find a site within this feature that will serve as a unique gRNA
    #             
    #             # for each orientation:
    #             # if (args.strands in ['+', 'both']):
    #             #  targets.extend...
    #             # if (args.strands in ['-', 'both']):
    #             #  targets.extend...
    #             
    #             # Code to enable automatic feature expansion
    #             # expansion_size = 10 # number of nt to expand features in 5' and 3' direction
    #             # minimum_spacers_per_feature = 100
    #             # found_spacers = 0
    #             # while ((found_spacers < minimum_spacers_per_feature) and ((feature_start, feature_end) != (0, contig_length)):
    #             #     for orientation in ['+', '-']:
    #             #         for motif in OnTargetMotif.motifs:
    #             #             for m in matches:
    #             #                 for target in filtered_targets:
    #             #                     cls() # need to add the novel (feature_start, feature_end) so it can be stored in this object
    #             #                     found_spacers += 1.0/len(filtered_targets) # Penalize for targets at the exact same location (IUPAC disambiguation)
    #             #     # Expand feature before next iteration
    #             #     feature_start = max(0, feature_start-expansion_size)
    #             #     feature_end = min(contig_length, feature_end + expansion_size)
    #             
    #             # Search both the '+' and '-' strands
    #             for orientation in ['+', '-']:
    #                 if (orientation == '+'):
    #                     sequence = contigs[feature_contig][feature_start:feature_end]
    #                 else:
    #                     sequence = nucleotides.rc(contigs[feature_contig][feature_start:feature_end])
    #                 
    #                 #for i in range(len(args.parsed_motifs)):
    #                 for motif in OnTargetMotif.motifs:
    #                     #spacers, pams, side = args.parsed_motifs[i]
    #                     spacers, pams, side = motif.parsed_list
    #                     #compiled_regex = args.compiled_motifs[i]
    #                     #matches = nucleotides.motif_search(sequence, spacers, pams, side)
    #                     matches = nucleotides.motif_search2(sequence, side, motif.compiled_regex)
    #                     for seq, start, end, spacer, pam in matches:
    #                         if (orientation == '+'):
    #                             real_start = feature_start + start
    #                             real_end = feature_start + end
    #                         else:
    #                             real_start = len(sequence) - end + feature_start
    #                             real_end = len(sequence) - start + feature_start
    #                         t_upstream = sequence[start-10:start]
    #                         t_downstream = sequence[end:end+10]
    #                         filtered_targets = target_filter(seq, spacer, pam, t_upstream, t_downstream, args)
    #                         for filt_seq, filt_spacer, filt_pam in filtered_targets:
    #                             cls(feature, feature_contig, orientation, real_start, real_end, t_upstream, t_downstream, filt_seq, side, filt_spacer, filt_pam, motif.motif_string, motif.parsed_list)
    #                             # Maybe add to ReversionDonor here
    #         else:
    #             logging.info("The contig '{}' for feature '{}' is not in the input FASTA.".format(feature_contig, feature))
    
    #def get_donors(self):
    #    return [ReversionDonor.indices[x] for x in self.get_contigs()]

class ReversionTarget(Target):
    prefix = 'reTarget'
    # Non-redundant dicts
    sequences = {} # key = nucleotide sequence, value = ReversionTarget object
    indices = {} # key = reTarget-234, value = ReversionTarget object
    
    @classmethod
    def get_targets(cls):
        """
        Extracts the spacer element from each ExcisionDonor object, and uses it
        to create the associated ReversionTarget objects.
        
        The constructor call adds the valid gRNA sites to 'ReversionTarget.sequences' and 'ReversionTarget.indices'.
        """
        for dDNA, obj in ExcisionDonor.sequences.items():
            #obj_features = ','.join([x[0] for x in list(obj.locations)])
            obj_features = ','.join(sorted(set([x[0] for x in obj.locations])))
            
            # Populate the ReversionTarget sequences indices dicts
            for t in obj.spacers: # [(orientation, start, end, filt_seq, side, filt_spacer, filt_pam), ...]
                #final_targets.append(tuple([obj.name, obj_feature, obj_contig] + list(t)))
                # t = (orientation, start, end, t_upstream, t_downstream, filt_seq, side, filt_spacer, filt_pam, args.motifs[i], args.parsed_motifs[i])
                ReversionTarget(obj_features, obj.name, t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], t[10])
    
    def get_donors(self):
        return [ExcisionDonor.indices[x] for x in self.get_contigs()]

class Motif(object):
    motifs = []
    def __init__(self, motif_string):
        # Save the input string: 'N{3,5}>NRG'
        self.motif_string = motif_string
        
        # Create [spacers, pams, side] list: (['NNN', 'NNNN', 'NNNNN'], ['NRG'], '>')
        self.parsed_list, self.spacer_sense_cuts, self.spacer_antisense_cuts, self.pam_sense_cuts, self.pam_antisense_cuts = self.parse_motif(self.motif_string)
        
        # Build the regex string: '((?:[ACGT][ACGT][ACGT])|(?:[ACGT][ACGT][ACGT][ACGT])|(?:[ACGT][ACGT][ACGT][ACGT][ACGT]))((?:[ACGT][AG]G))'
        self.regex_string = self.build_motif_regex(self.parsed_list[0], self.parsed_list[1], self.parsed_list[2])
        
        # Compile the regex
        self.compiled_regex = self.compile_regex(self.regex_string, ignore_case=True)
        
        # Save object to publicly-accessible, class list
        self.motifs.append(self)
        
        # Example:
        # motif          TTTN<N{19}/.{4}\
        # genome      ...CCAAAGCCAACGACTTTAGCTAGCTAAAGGACCTATGCCCATTACATGCCGCCAA...
        # sequence                     TTTAGCTAGCTAAAGGACCTATGCCCA
        # upstream           AGCCAACGAC
        # pam                          TTTA
        # target                           GCTAGCTAAAGGACCTATG
        # cut sites                                          ><  ><
        # sense cuts                                         ><
        # antisense cuts                                         ><
        # double-strand cuts
        # downstream                                              TTACATGCCG
        
    def parse_motif(self, motif):
        """
        Takes string of SPACER>PAM and returns ([spacer, spacer, ...], [pam, pam ...], '>')
        
        Example inputs: 'G{,2}N{19,20}>NGG', 'N{20,21}>NNGRRT', 'TTTN<N{20,23}'
        """
        
        # Still need to implement:
        #  Identify gRNA pairs on opposite strands appropriate for "nickase" Cas9
        #  allow for PAM-out and PAM-in orientations
        #  --motifs, where '.' is any non-spacer nucleotide (from the "Dimeric CRISPR RNA-guided FokI nucleases for highly specific genome editing" paper)
        #    'CCN<N{20}.{13,18}N{20}>NGG'     # PAM-out
        #    'N{20}>NGG.{13,18}CCN<N{20}'     # PAM-in
        #
        # The FokI-dCas9 restriction site can cut anywhere between either PAM
        # site (for PAM-out). It typically cuts with a 3' overhang
        # (Targeting individual subunits of the FokI restriction endonuclease to specific DNA strands):
        #   5'-.../.... ...-3'
        #   3'-... ....\...-5'
        #
        # Maybe this? It isn't perfect, but it lets the user specify the cut sites.
        #   'CCN<N{20}.{5,9}/....\.{9,5}N{20}>NGG' # One quantifier counts up while the other counts down?
        
        gt_count = motif.count('>')
        lt_count = motif.count('<')
        lb_count = motif.count('{')
        rb_count = motif.count('}')
        # Make sure motif does not violate basic rules
        if (gt_count + lt_count < 1):
            raise Exception("Motif lacks distinction between spacer and PAM sequences ('>' or '<' character)")
        elif (gt_count + lt_count > 1):
            raise Exception("Motif has too many '>' or '<' characters")
        if (lb_count != rb_count):
            raise Exception("Motif braces '{' and '}' do not match")
        if (motif.count(' ') > 0):
            raise Exception("Motif contains invalid space ' ' characters")
        if (motif.count('{}') > 0):
            raise Exception("Motif contains invalid quantifier '{}'")
        if (motif.count('{,}') > 0):
            raise Exception("Motif contains invalid quantifier '{,}'")
        
        if (gt_count == 1):
            spacer_motif, pam_motif = motif.split('>')
            side = '>'
        elif (lt_count == 1):
            pam_motif, spacer_motif = motif.split('<')
            side = '<'
        
        spacer_sequences, spacer_sense_cuts, spacer_antisense_cuts = self.parse_motif_helper(spacer_motif)
        pam_sequences, pam_sense_cuts, pam_antisense_cuts = self.parse_motif_helper(pam_motif)
        
        return (spacer_sequences, pam_sequences, side), spacer_sense_cuts, spacer_antisense_cuts, pam_sense_cuts, pam_antisense_cuts
    
    # List of common motifs obtained from
    #  https://github.com/maximilianh/crisporWebsite/crispor.py
    
    # from (https://benchling.com/pub/cpf1):
    # Cpf1 is an RNA-guided nuclease, similar to Cas9. It recognizes a T-rich
    # PAM, TTTN, but on the 5' side of the guide. This makes it distinct from
    # Cas9, which uses an NGG PAM on the 3' side. The cut Cpf1 makes is
    # staggered. In AsCpf1 and LbCpf1, it occurs 19 bp after the PAM on the
    # targeted (+) strand and 23 bp on the other strand, as shown here:
    #                        Cpf1
    #   TTTC GAGAAGTCATCTAATAAGG|CCAC TGTTA
    #   AAAG CTCTTCAGTAGATTATTCC GGTG|ACAAT
    #   =PAM =========gRNA====== =
    #
    # Benchling suggests 20 nt guides can be used for now, as there is not real
    # suggestion for optimal guide length. Robust guide scores for Cpf1 are
    # still in development, but simple scoring based on the number of off-target
    # sites is available on Benchling.
    # 
    # Cpf1 requires only a crRNA for activity and does not need a tracrRNA to
    # also be present.
    #
    # Two Cp1-family proteins, AsCpf1 (from Acidaminococcus)
    # and LbCpf1 (from Lachnospiraceae), have been shown to perform efficient
    # genome editing in human cells.
    #
    # Why use Cpf1 over Cas9?
    #  see https://benchling.com/pub/cpf1
    
    def parse_motif_helper(self, submotif):
        """
        Helper function that parses either the SPACER or PAM motif.
        Decodes quantifiers and returns a list
        """
        # Keep track of expanded sequences
        sequences = ['']
        
        # Keep track if a quantifier is being parsed
        quantifier = None
        
        # Iterate through the characters
        for i, c in enumerate(submotif):
            if (c == '{'): # Start the quantifier
                quantifier = ''
            elif (c == '}'): # End the quantifier
                quantifier_list = quantifier.split(',')
                if (len(quantifier_list) == 1):
                    min_length = int(quantifier)
                    max_length = min_length
                    
                elif (len(quantifier_list) == 2):
                    if (quantifier_list[0] == ''):
                        min_length = 0
                    else:
                        min_length = int(quantifier_list[0])
                    if (quantifier_list[1] == ''):
                        raise Exception("Motif quantifier '{" + quantifier + "}' contains no maximum value")
                    else:
                        max_length = int(quantifier_list[1])
                    if (min_length > max_length):
                        raise Exception("Motif quantifier '{" + quantifier + "}' minimum and maximum lengths are invalid")
                
                last_chars = [ x[-1] for x in sequences ]
                
                sequences = [ x[:-1] for x in sequences ]
                
                new_sequences = []
                for j, s in enumerate(sequences):
                    for length in range(min_length, max_length+1):
                        new_sequences.append(s + last_chars[j]*length)
                sequences = new_sequences
                quantifier = None
            elif (quantifier != None): # add current character to quantifier if it is open
                # We strip out any cut sites
                #if (c not in ['/', '\\', '|']):
                quantifier += c
            else: # add the current character to the expanded sequences
                for j in range(len(sequences)):
                    sequences[j] = sequences[j] + c
        
        # Define empty list holding all cut sites
        sense_cuts = []
        antisense_cuts = []
        
        # Iterate through sequences, and record the location of any cut sites
        # Then strip out any cut sites
        for i, seq in enumerate(sequences):
            sense_cuts.append([])
            antisense_cuts.append([])
            found = 0
            
            nseq = seq
            for j, c in enumerate(seq):
                if (c  == '/'):
                    sense_cuts[-1].append(j-found)
                    nseq = nseq.replace('/', '', 1)
                    found += 1
                elif (c == '\\'):
                    antisense_cuts[-1].append(j-found)
                    nseq = nseq.replace('\\', '', 1)
                    found += 1
                elif (c == '|'):
                    sense_cuts[-1].append(j-found)
                    antisense_cuts[-1].append(j-found)
                    nseq = nseq.replace('|', '', 1)
                    found += 1
            sequences[i] = nseq
        
        return sequences, sense_cuts, antisense_cuts
    
    def compile_regex(self, re_pattern, ignore_case=False, enhance_match=False, best_match=False):
        """
        Returns compiled regex
        """
        # Choose the regex flags
        myflags = 0
        if ignore_case:
            myflags |= regex.IGNORECASE
        if enhance_match:
            myflags |= regex.ENHANCEMATCH
        if best_match:
            myflags |= regex.BESTMATCH
            
        c = regex.compile(re_pattern, flags=myflags) # regex.ENHANCEMATCH|regex.IGNORECASE
        return c
    
    def build_motif_regex(self, spacers, pams, side, anchored=False):
        """
        Returns regex pattern string for motif
        """
        spacer_pattern = '(' + '|'.join([self.build_regex_pattern(x, capture=False) for x in spacers]) + ')'
        pam_pattern = '(' + '|'.join([self.build_regex_pattern(x, capture=False) for x in pams]) + ')'
        if (side == '>'): # SPACER>PAM
            if anchored:
                re_pattern = '^'+ spacer_pattern + pam_pattern + '$'
            else:
                re_pattern = spacer_pattern + pam_pattern
        elif (side == '<'): # PAM<SPACER
            if anchored:
                re_pattern = '^'+ pam_pattern + spacer_pattern + '$'
            else:
                re_pattern = pam_pattern + spacer_pattern
        return re_pattern
    
    def build_regex_pattern(self, iupac_sequence, max_substitutions=0, max_insertions=0, max_deletions=0, max_errors=0, capture=True):
        """
        Build a regular expression pattern for the nucleotide search, taking IUPAC
        ambiguities into account.
        
        So far for DNA only
        Expects string as input (NOT a list)
        """
        # Format the fuzzy matching restrictions
        fuzzy = []
        if (max_substitutions > 0):
            fuzzy.append('s<={}'.format(max_substitutions))
        if (max_insertions > 0):
            fuzzy.append('i<={}'.format(max_insertions))
        if (max_deletions > 0):
            fuzzy.append('d<={}'.format(max_deletions))
        if (max_errors > 0):
            fuzzy.append('e<={}'.format(max_errors))
        if (len(fuzzy) > 0):
            fuzzy = '{' + ','.join(fuzzy) + '}'
        else:
            fuzzy = ''
        
        # Convert the IUPAC sequence to an equivalent regex
        iupac = {
            'a': 'a',
            'c': 'c',
            'g': 'g',
            't': 't',
            'r': '[ag]',
            'y': '[ct]',
            'm': '[ac]',
            'k': '[gt]',
            'w': '[at]',
            's': '[cg]',
            'b': '[cgt]',
            'd': '[agt]',
            'h': '[act]',
            'v': '[acg]',
            'n': '[acgt]',
            
            'A': 'A',
            'C': 'C',
            'G': 'G',
            'T': 'T',
            'R': '[AG]',
            'Y': '[CT]',
            'M': '[AC]',
            'K': '[GT]',
            'W': '[AT]',
            'S': '[CG]',
            'B': '[CGT]',
            'D': '[AGT]',
            'H': '[ACT]',
            'V': '[ACG]',
            'N': '[ACGT]',
            
            '.': '.', # Any character
            
            '|': r'\|', # Double-stranded cut
            '/': r'\/', # Sense cut
            '\\': r'\\', # anti-sense cut
        }
        sequence = ''.join(map(lambda x: iupac[x], iupac_sequence))
        if capture:
            pattern = '(' + sequence + ')' + fuzzy
        else:
            pattern = '(?:' + sequence + ')' + fuzzy
        #logger.info('Built regex string: {!r}'.format(pattern))
        return pattern
    
    def generate_sequence(self, compositions=None, complement=False, samples=100):
        """
        Generate a random sequence according to input sequence composition
        or with uniform composition of None.
        """
        
        if compositions:
            seqs = []
            scores = []
            for i in range(samples):
                s = bartag.random_motif_sequence(self.motif_string)
                seqs.append(s)
                scores.append(nucleotides.sequence_likelihood(s, compositions)/len(s))
            
            #sorted_seqs = sorted(seqs, key=lambda x: seqs[x], reverse=not complement)
            #return sorted_seqs[0]
            if complement:
                return min(zip(scores, seqs))[1]
            else:
                return max(zip(scores, seqs))[1]
        else:
            return bartag.random_motif_sequence(self.motif_string)
        

class Feature(object):
    features = {}
    excluded_features = {}
    # Origin flags
    NONE=0
    INPUT=1
    DERIVED=2
    def __init__(self, contig, start, end, strand, name=None, attributes=None, source=None, feature_type=None, score=None, frame=None, origin=NONE, sep=';', expand_parent=None):
        self.contig = contig
        self.source = source
        self.feature_type = feature_type
        self.start = start
        self.end = end
        if (score != '.'):
            self.score = score
        else:
            self.score = None
        self.strand = strand
        if (frame != '.'):
            self.frame = frame
        else:
            self.frame = None
        self.attributes = {}
        if isinstance(attributes, dict):
            self.attributes = attributes
        elif isinstance(attributes, str): # "ID=12;Parent=nope"
            alist = regex.split(sep+'\s*', attributes) # ['ID=12', 'Parent=nope']
            self.attributes = dict([regex.split('\s*=\s*', x) for x in alist]) # {'ID': '12', 'Parent': 'nope'}
        self.name = name
        self.origin = origin
        self.expand_parent = expand_parent
    
    def get_expand_parent(self):
        if (self.expand_parent == None):
            return self
        else:
            return self.expand_parent.get_expand_parent()
    
    @classmethod
    def parse_gff_line(cls, text):
        obj = None
        
        line = text.rstrip()
        if not line.startswith('#'):
            sline = line.split('\t')
            if (len(sline) > 6):
                contig = sline[0]
                source = sline[1]
                feature_type = sline[2]
                start = int(sline[3])-1 # Should always be true: start <= end 
                end = int(sline[4]) # Should always be true: start <= end
                score = sline[5]
                strand = sline[6]
                frame = None
                attributes = None
                if (len(sline) > 7):
                    frame = sline[7]
                if (len(sline) > 8):
                    attributes = sline[8]
                obj = Feature(contig, start, end, strand, source=source, feature_type=feature_type, score=score, frame=frame, attributes=attributes, origin=Feature.INPUT)
        return obj
    
    @classmethod
    def load_gff_file(cls, filename, feature_types, excluded_feature_types, selected_features, tag):
        """
        Load General Feature Format (GFF) file into dict()
        One line per feature, each containing 9 columns of data, plus optional
        track definition lines.
        
        Converts positions to 0-based index.
        """
        # Fields must be tab-separated. Also, all but the final field in each
        # feature line must contain a value; "empty" columns should be denoted
        # with a '.'
        #  1) seqid - name of the chromosome or scaffold; chromosome names can
        #     be given with or without the 'chr' prefix. Important note: the
        #     seqname must be a standard chromosome name or an identifier such as
        #     a scaffold ID, without any additional content such as species or
        #     assembly.
        #  2) source - name of the program that generated this feature, or the
        #     data source (database or project name)
        #  3) feature - feature type name, e.g. Gene, Variation, Similarity
        #  4) start - Start position of the feature, with sequence numbering
        #     starting at 1.
        #  5) end - End position of the feature, with sequence numbering
        #     starting at 1.
        #  6) score - A floating point value. As in earlier versions of the format,
        #     the semantics of the score are ill-defined. It is strongly
        #     recommended that E-values be used for sequence similarity features,
        #     and that P-values be used for ab initio gene prediction features.
        #     If there is no score, put a '.' (a period) in this field.
        #  7) strand - defined as '+' (forward), '-' (reverse), '.' (unstranded),
        #     '?' (relevant, but unknown).
        #  8) frame - for CDS features, '0', '1' or '2'. '0' indicates that the first base of
        #     the feature is the first base of a codon, '1' that the second base
        #     is the first base of a codon, and so on. Other features can use '.'.
        #  9) attribute - A semicolon-separated list of tag-value pairs, providing
        #     additional information about each feature. A list of feature
        #     attributes in the format tag=value. Multiple tag=value pairs are
        #     separated by semicolons. URL escaping rules are used for tags or
        #     values containing the following characters: ",=;". Spaces are allowed
        #     in this field, but tabs must be replaced with the %09 URL escape.
        #     This field is not required.
        #     Column 9 tags have predefined meanings:
        #       ID - Indicates the unique identifier of the feature. IDs must be
        #            unique within the scope of the GFF file.
        #       Name - Display name for the feature. This is the name to be displayed to the user. Unlike IDs, there is no requirement that the Name be unique within the file.
        #       Alias - A secondary name for the feature. It is suggested that this tag be used whenever a secondary identifier for the feature is needed, such as locus names and accession numbers. Unlike ID, there is no requirement that Alias be unique within the file.
        #       Parent - Indicates the parent of the feature. A parent ID can be used to group exons into transcripts, transcripts into genes, and so forth. A feature may have multiple parents. Parent can *only* be used to indicate a partof relationship.
        #       Target - Indicates the target of a nucleotide-to-nucleotide or protein-to-nucleotide alignment. The format of the value is "target_id start end [strand]", where strand is optional and may be "+" or "-". If the target_id contains spaces, they must be escaped as hex escape %20.
        #       Gap - The alignment of the feature to the target if the two are not collinear (e.g. contain gaps). The alignment format is taken from the CIGAR format described in the Exonerate documentation. http://cvsweb.sanger.ac.uk/cgi-bin/cvsweb.cgi/exonerate?cvsroot=Ensembl). See the GFF3 specification for more information.
        #       Derives_from - Used to disambiguate the relationship between one feature and another when the relationship is a temporal one rather than a purely structural "part of" one. This is needed for polycistronic genes. See the GFF3 specification for more information.
        #       Note - A free text note.
        #       Dbxref - A database cross reference. See the GFF3 specification for more information.
        #       Ontology_term - A cross reference to an ontology term. See the GFF3 specification for more information.
        #     Multiple attributes of the same type are indicated by separating the
        #     values with the comma "," character, as in: 'Parent=AF2312,AB2812,abc-3'
        #     Note that attribute names are case sensitive. "Parent" is not the
        #     same as "parent". All attributes that begin with an uppercase letter
        #     are reserved for later use. Attributes that begin with a lowercase
        #     letter can be used freely by applications. You can stash any
        #     semi-structured data into the database by using one or more
        #     unreserved (lowercase) tags.
        
        with open(filename, 'r') as flo:
            for line in flo:
                line = line.rstrip()
                obj = cls.parse_gff_line(line)
                if obj:
                    if ((obj.feature_type in feature_types) or ('all' in feature_types)):
                        obj.name = obj.attributes[tag]
                        
                        # Reduce the total number of features to just the ones indicated in the selection
                        if selected_features:
                            if (obj.name in selected_features):
                                cls.features[obj.name] = obj
                        else:
                            cls.features[obj.name] = obj
                    elif ((obj.feature_type in excluded_feature_types) or ('all' in excluded_feature_types)):
                        obj.name = obj.attributes[tag]
                        cls.excluded_features[obj.name] = obj
        
        #logger.info('GFF file parsed: {!r}'.format(filename))
    
    @classmethod
    def assert_features(cls, selection, contigs):
        # Current implementation limitations
        # GFF file only has 'Gene=' tag on ONE of the homologs, and not the other
        # User will have to specify the feature twice if he wants to target both homologs
        
        # Require at least one in order to filter
        # Feature has this format: key=gene/tag, value = (contig, start(bp), end(bp), strand)
        if selection:
            for s in selection:
                f = cls.features.get(s)
                if not f:
                    raise Exception("Feature specified as '--selection "+s+"' does not exist in input '--gff'.")
        
        for feature_name, f in cls.features.items():
            if not f.contig in contigs:
                raise Exception("Feature '{}' lies on contig '{}' which does not exist in input '--fasta'.".format(feature_name, f.contig))
        
        if (len(cls.features) == 0):
            raise Exception("Input '--gff' and '--selection' combination returns no valid features.")
    
    @classmethod
    def expand_all_features(cls, args, contigs):
        checked_keys = []
        # Some make-shift code to get this working when the 'Feature.features' dict changes size during the loop process
        while(len(checked_keys) < len(Feature.features)):
            for k, v in Feature.features.items():
                if (k not in checked_keys):
                    feature_name = k
                    f = v
                    break
            checked_keys.append(feature_name)
            
            contig_sequence = contigs[f.contig]
            if (f.origin == Feature.INPUT):
                #ef = f.expand_feature(args, contig_sequence)
                logging.info('Expanding feature: {}'.format(feature_name))
                f.expand_feature(args, contig_sequence)
    
    def expand_feature(self, args, contig_sequence):
        """
        New method of expanding a feature.
        Will generate all possible expanded features, which will be
        evaluated later to determine which is the best.
        """
        # Option 1 (Can't do because target evaluation happens later)
        #  keep expanding the feature until a target of minimum quality is found
        #  or until the maximum size 'args.feature_expansion_lengths[1]' is reached
        
        # Option 2 (Can't do because target evaluation happens later)
        #  scan the widest region for targets
        #  evaluate target scores/weights
        #  pick the best target
        #  center/justify the FEATURE and TARGET
        #  create the derived feature object
        
        # Option 3
        #  create all derived features within 'args.feature_expansion_lengths' limits
        #  these are centered/justified already around each target/feature (as specified by user command line options)
        #  later:
        #   evaluate the target scores/weights
        #   pick the feature with the best target
        
        contig_length = len(contig_sequence)
        
        # First we identify all the features from -max to +max
        feature_start = self.start
        feature_end = self.end or contig_length # if (self.end == None), then it will be len(contig_sequence)
        
        # Find the maximum distance the feature can be expanded both up- and down-stream:
        max_upstream_coord = 0
        max_downstream_coord = contig_length
        for exf_name, exf_obj in Feature.excluded_features.items():
            if (self.contig == exf_obj.contig):
                if (exf_obj.end < self.start):
                    max_upstream_coord = max(max_upstream_coord, exf_obj.end)
                if (exf_obj.start > self.end):
                    max_downstream_coord = min(max_downstream_coord, exf_obj.start)
                if (Feature.overlap_coverage(self.start, self.end, exf_obj.start, exf_obj.end) > 0):
                    logging.info("WARNING: Selected feature '{}' overlaps with excluded feature '{}'".format(self.name, exf_obj.name))
        
        logging.info('    feature: {}'.format(self.name))
        logging.info('     bounds: {}..{}'.format(self.start, self.end))
        logging.info('     limits: {}..{}'.format(max_upstream_coord, max_downstream_coord))
        
        
        
        #feature_sequence = contig_sequence[feature_start:feature_end]
        # The '10' is for upstream/downstream adjacent sequences
        #ex_feature_start = max(0, feature_start-(10+args.feature_expansion_lengths[1]))
        #ex_feature_end = min(feature_end+(10+args.feature_expansion_lengths[1]), contig_length)
        ex_feature_start = max(0, feature_start-args.feature_expansion_lengths[1])
        ex_feature_end = min(feature_end+args.feature_expansion_lengths[1], contig_length)
        
        #max_feature_sequence = contig_sequence[ex_feature_start:ex_feature_end]
        
        # We scan for targets only once
        targets = Target.get_targets(args, contig_sequence, start=ex_feature_start, end=ex_feature_end) # Does both orientations (+/-)
        logging.info('max targets: {}'.format(len(targets)))
        
        # The relative coordinates of FEATURE within EX_FEATURE:
        #relative_feature_start = feature_start - ex_feature_start
        #relative_feature_end = feature_end - ex_feature_start
        
        derived_sets = set()
        
        # We generate a derived feature for each identified target
        for t in targets:
            # t = (orientation, start, end, upstream, downstream, filt_seq, side, filt_spacer, filt_pam, mymotif.motif_string, tuple([tuple(x) if isinstance(x, list) else x for x in mymotif.parsed_list])))
            
            target_start = t[1]
            target_end = t[2]
            
            #  center_feature: --------HHHH[...............FEATURE.........TARGET]HHHH------------------------
            #                  -----HHHH[..................FEATURE.........TARGET...]HHHH--------------------- pad=3
            #   center_target: -----------------------HHHH[FEATURE.........TARGET................]HHHH--------
            #                  --------------------HHHH[...FEATURE.........TARGET...................]HHHH----- pad=3
            #     center_both: -----------------------HHHH[FEATURE.........TARGET]HHHH------------------------
            #                  --------------------HHHH[...FEATURE.........TARGET...]HHHH--------------------- pad=3
            # justify_feature: -----------------------HHHH[FEATURE.........TARGET]HHHH------------------------
            #                  -----------------------HHHH[FEATURE.........TARGET...]HHHH--------------------- pad=3
            #  justify_target: -----------------------HHHH[FEATURE.........TARGET]HHHH------------------------
            #                  --------------------HHHH[...FEATURE.........TARGET]HHHH------------------------ pad=3
            
            # max(end2-start1, end1-start2) # Full distance including X and Y = 19
            # max(start2-end1, start1-end2) # Distance in between between X and Y = 3
            #     1 XXXXXXXXXX
            #     2              YYYYYY
            #    
            #     1          XXXXXXXXXX
            #     2 YYYYYY
            
            if (args.feature_expansion_method == 'center_feature'):
                # ----[........FFF....TTTT]----  Terminal 3'
                # ----[TTTT....FFF........]----  Terminal 5'
                # ----------[..fffTT]----------  Overlapping 3'
                # ----------[TTfff..]----------  Overlapping 5'
                # ----------[.TfffTT]----------  Complete overlap 3' longer
                # ----------[TTfffT.]----------  Complete overlap 5' longer
                full_dist = max(target_end-feature_start, feature_end-target_start)
                derived_start = feature_end - full_dist - args.feature_expansion_pad
                derived_end = feature_start + full_dist + args.feature_expansion_pad
            
            elif (args.feature_expansion_method == 'center_target'):
                # ----[....TTTT.FFF]----
                # ----[FFF.TTTT....]----
                full_dist = max(target_end-feature_start, feature_end-target_start)
                derived_start = target_end - full_dist - args.feature_expansion_pad
                derived_end = target_start + full_dist + args.feature_expansion_pad
            
            elif (args.feature_expansion_method == 'center_both'):
                # ----[TTTT...FFF]----
                # ----[FFF...TTTT]----
                derived_start = min(feature_start, target_start) - args.feature_expansion_pad
                derived_end = max(target_end, feature_end) + args.feature_expansion_pad
            
            elif (args.feature_expansion_method == 'justify_feature'):
                # ----[FFF....TTTTxxx]----
                # ----[xxxTTTT....FFF]----
                derived_start = min(feature_start, target_start)
                derived_end = max(target_end, feature_end)
                if (derived_start == feature_start):
                    derived_end += args.feature_expansion_pad
                else:
                    derived_start -= args.feature_expansion_pad
                
            elif (args.feature_expansion_method == 'justify_target'):
                # ----[TTTT....FFFxxx]----
                # ----[xxxFFF....TTTT]----
                derived_start = min(feature_start, target_start)
                derived_end = max(target_end, feature_end)
                if (derived_start == target_start):
                    derived_end += args.feature_expansion_pad
                else:
                    derived_start -= args.feature_expansion_pad
            
            # If derived feature is too small, then we automatically pad it
            # so it reaches the minimum length
            if (derived_end - derived_start < args.feature_expansion_lengths[0]):
                size_difference = args.feature_expansion_lengths[0] - (derived_end - derived_start)
                if (args.feature_expansion_method in ['center_feature', 'center_target', 'center_both']):
                    new_pad = (size_difference+1)//2
                    derived_start -= new_pad
                    derived_end += new_pad
                    logging.info("Added {} nt of padding bases to either side of DERIVED FEATURE".format(new_pad))
                
                elif (args.feature_expansion_method == 'justify_feature'):
                    if (derived_start == feature_start):
                        derived_end += size_difference
                    else:
                        derived_start -= size_difference
                    logging.info("Added {} nt of padding bases to one side of DERIVED FEATURE".format(size_difference))
                
                elif (args.feature_expansion_method == 'justify_target'):
                    if (derived_start == target_start):
                        derived_end += size_difference
                    else:
                        derived_start += size_difference
                    logging.info("Added {} nt of padding bases to one side of DERIVED FEATURE".format(size_difference))
            
            logging.info("derived_start = {}".format(derived_start))
            logging.info("  derived_end = {}".format(derived_end))
            logging.info("derived_end - derived_start = {}".format(derived_end - derived_start))
            derived_sets.add((derived_start, derived_end))
        
        # Create a counter for the derived feature
        count = 0
        
        for derived_start, derived_end in derived_sets:
            # If derived feature length is within the desired size range, then create the derived feature
            if (args.feature_expansion_lengths[0] <= derived_end - derived_start <= args.feature_expansion_lengths[1]):
                if (max_upstream_coord <= derived_start <= derived_end <= max_downstream_coord):
                    new_name = self.name + '_derived-' + str(count)
                    new_attributes = self.attributes.copy()
                    new_attributes[args.tag] = new_name
                    new_feature = Feature(self.contig, derived_start, derived_end, self.strand, name=new_name, source=self.source, feature_type=self.feature_type, score=self.score, frame=self.frame, attributes=new_attributes, origin=Feature.DERIVED, expand_parent=self)
                    Feature.features[new_name] = new_feature
                    count += 1
                    logging.info("DERIVED FEATURE '{}' created.".format(new_name))
                else:
                    logging.info("DERIVED FEATURE would overlap with excluded feature.")
        # END expand_feature()
    
    def previous_expand_feature(self, args, contig_sequence, expansion_size=20, minimum_targets_per_feature=5):
        feature_start = self.start
        feature_end = self.end
        
        contig_length = len(contig_sequence)
        
        if (feature_end == None):
            feature_end = len(contig_sequence)
        
        # Find the maximum distance the feature can be expanded both up- and down-stream:
        max_upstream_coord = 0
        max_downstream_coord = contig_length
        for exf_name, exf_obj in Feature.excluded_features.items():
            if (self.contig == exf_obj.contig):
                if (exf_obj.end < self.start):
                    max_upstream_coord = max(max_upstream_coord, exf_obj.end)
                if (exf_obj.start > self.end):
                    max_downstream_coord = min(max_downstream_coord, exf_obj.start)
                if (Feature.overlap_coverage(self.start, self.end, exf_obj.start, exf_obj.end) > 0):
                    logging.info("WARNING: Selected feature '{}' overlaps with excluded feature '{}'".format(self.name, exf_obj.name))
        
        logging.info('feature: {}'.format(self.name))
        logging.info( 'bounds: {}..{}'.format(self.start, self.end))
        logging.info(' limits: {}..{}'.format(max_upstream_coord, max_downstream_coord))
        
        # Gradually expand feature size until the minimum number of targets is found
        feature_sequence = contig_sequence[feature_start:feature_end]
        targets = Target.get_targets(args, feature_sequence) # Does both orientations (+/-)
        logging.info('targets: {}'.format(len(targets)))
        count = 0
        while ((len(targets) < minimum_targets_per_feature) and (((feature_start, feature_end) != (0, contig_length)) or ((feature_start, feature_end) != (max_upstream_coord, max_downstream_coord)))):
            feature_start = max(0, feature_start-expansion_size, max_upstream_coord)
            feature_end = min(contig_length, feature_end + expansion_size, max_downstream_coord)
            
            feature_sequence = contig_sequence[feature_start:feature_end]
            targets = Target.get_targets(args, feature_sequence) # Does both orientations (+/-)
            logging.info('targets: {}'.format(len(targets)))
            count += 1
        
        # Save the new feature as a Feature object
        if (count > 0):
            new_name = self.name + '_derived'
            new_attributes = self.attributes.copy()
            new_attributes[args.tag] = new_name
            new_feature = Feature(self.contig, feature_start, feature_end, self.strand, name=new_name, source=self.source, feature_type=self.feature_type, score=self.score, frame=self.frame, attributes=new_attributes, origin=Feature.DERIVED, expand_parent=self)
            Feature.features[new_name] = new_feature
            return new_feature
        else:
            return self
    
#    def overlap_distance(self, start1, end1, start2, end2):
#        coverage = self.overlap_coverage(start1, end1, start2, end2)
#        if (coverage > 0):
#            return -coverage
#        else:
#            #return min(abs(start2-end1), abs(start1-end2))
#            return max(start2-end1, start1-end2)
#            # 1 ..........
#            # 2              XXXXXX
#            
#            # 1          ........
#            # 2 XXXXXX
    
#    def overlap_percent_coverage(self, start1, end1, start2, end2):
#        coverage = self.overlap_coverage(start1, end1, start2, end2)
#        return float(coverage) / max(abs(end1 - start1), abs(end2 - start2))
    
#    def percent_full_length_coverage(self, start1, end1, len1, start2, end2, len2):
#        coverage = self.overlap_coverage(start1, end1, start2, end2)
#        return float(coverage) / max(len1, len2)
    
    @staticmethod
    def overlap_coverage(start1, end1, start2, end2, index_base=0):
        coverage = 0
        
        # Expects 0-based, left-inclusive, right-exclusive indexing
        if (start2 <= start1 <= end1 <= end2):
            # 1      .............
            # 2   XXXXXXXXXXXXXXXXXX
            coverage = abs(end1 - start1)
        elif (start1 <= start2 <= end2 <= end1):
            # 1      .....................
            # 2            XXXXXXXXXX
            coverage = abs(end2 - start2)
        elif (start1 <= start2 <= end1 <= end2):
            # 1      ........
            # 2        XXXXXXXX
            coverage = abs(end1 - start2)
        elif (start2 <= start1 <= end2 <= end1):
            # 1      .............
            # 2  XXXXXXXXXX
            coverage = abs(end2 - start1)
        
        # if 1-based, left-inclusive, right-inclusive indexing,
        # then give +1 to coverage
        if index_base == 1:
            coverage += 1
        
        return coverage
    
    @classmethod
    def get_overlapping_features(cls, contig, start, end):
        """Returns list (sorted) of features that overlap with arguments"""
        the_features = set()
        for feature_name, f in cls.features.items():
            if (contig == f.contig):
                if (cls.overlap_coverage(start, end, f.start, f.end) > 0):
                    the_features.add(feature_name)
        
        return sorted(the_features)

class OnTargetMotif(Motif):
    motifs = []

class OffTargetMotif(Motif):
    motifs = []

class CustomHelpFormatter(argparse.HelpFormatter):
    """Help message formatter which retains any formatting in descriptions
    and adds default values to argument help.
    
    Only the name of this class is considered a public API. All the methods
    provided by the class are considered an implementation detail.
    """
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

class main(object):
    def __init__(self):
        """Create argument parser, then run the selected sub-program"""
        
        # Obtain command line arguments and parse them
        args = self.parse_arguments()
        
        self.process_general_arguments(args)
        
        # Get timestamp for analysis beginning
        start_time = time.time()
        
        # Echo the command line parameters to STDOUT
        #print(args, flush=True)
        
        # Print command line parameters (arguments) to log file
        if hasattr(args, 'folder'):
            logging.info(args)
        
        # call the function for which action was used
        args.func(args)
        
        if hasattr(args, 'folder'):
            # Print time taken for program to complete
            logging.info('{} finished'.format(__program__))
            logging.info('Runtime: {}s'.format(time.time()-start_time))
        
    def _glossary(self, args):
        """Print the Glossary"""
        utils.print_local_file('glossary.txt')
    
    def _motifs(self, args):
        """Print the list of common CRISPR/Cas motifs"""
        utils.print_local_file('motifs.txt')
    
    def _algorithms(self, args):
        """Print information on the algorithms"""
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
    
    def _aligners(self, args):
        """Print information about the supported aligners"""
        # Aligners intended to implement = ['addtag', 'blast+', 'blat', 'bowtie', 'bowtie2', 'bwa', 'cas-offinder']
        # Other aligners to consider: 'rmap', 'maq', 'shrimp2', 'soap2', 'star', 'rhat', 'mrsfast', 'stampy'
        for x in aligners.aligners:
            print('==========', x.name, '==========')
            print('     Author:', x.author)
            print('       Year:', x.year)
            print('   Citation:', x.citation)
            print('      Input:', x.input)
            print('     Output:', x.output)
            print('  Truncated:', x.truncated)
            print('')
    
    def _oligos(self, args):
        """Print information about the supported oligos"""
        for x in oligos.oligos:
            print('==========', x.name, '==========')
            print('     Author:', x.author)
            print('       Year:', x.year)
            print('   Citation:', x.citation)
            print('')
    
    def _search(self, args):
        """
        Search input FASTA for arbitrary IUPAC sequence
        Print out a GFF3 that can be used as an input to "generate"
        """
        if not args.identifier:
            args.identifier = ["search_"+str(x) for x in range(len(args.query))]
        if (len(args.query) != len(args.identifier)):
            raise Exception("The number of sequences specified in '--query' does not match the number of features in '--identifier'.")
        
        myflags = regex.ENHANCEMATCH | regex.IGNORECASE # regex.ENHANCEMATCH|regex.IGNORECASE|regex.BESTMATCH
        
        contig_sequences = utils.load_multiple_fasta_files(args.fasta)
        
        for i, qseq in enumerate(args.query):
            qpat = nucleotides.build_regex_pattern(qseq) #max_substitutions=0, max_insertions=0, max_deletions=0, max_errors=0, capture=True
            c = regex.compile(qpat, flags=myflags)
            
            rc_qseq = nucleotides.rc(qseq)
            rc_qpat = nucleotides.build_regex_pattern(rc_qseq)
            rc_c = regex.compile(rc_qpat, flags=myflags)
            
            j = 0
            
            for ctg_name, ctg_seq in contig_sequences.items():
                for m in c.finditer(ctg_seq):
                    # contig_id, source, feature_type, start, end, score, strand, frame, attributes
                    contig_id = ctg_name
                    source = "addtag"
                    feature_type = "search"
                    start = m.start()+1
                    end = m.end()
                    score = '.' # <------ replace with E-value or percent identity
                    strand = '+'
                    frame = '.'
                    attributes = args.tag+'='+args.identifier[i] + '_' + str(j)
                    print("\t".join(str(x) for x in [contig_id, source, feature_type, start, end, score, strand, frame, attributes]))
                    j += 1
                
                for m in rc_c.finditer(ctg_seq):
                    contig_id = ctg_name
                    source = "addtag"
                    feature_type = "search"
                    start = m.start()+1
                    end = m.end()
                    score = '.' # <------ replace with E-value or percent identity
                    strand = '-'
                    frame = '.'
                    attributes = args.tag+'='+args.identifier[i] + '_' + str(j)
                    print("\t".join(str(x) for x in [contig_id, source, feature_type, start, end, score, strand, frame, attributes]))
                    j += 1
        # End _search()
    
    def _feature(self, args):
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
                obj = Feature.parse_gff_line(line)
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
                                m = regex.search(q_pattern, obj.attributes[k], flags=regex.IGNORECASE)
                                if m:
                                    #print(m)
                                    #matched_lines.add((line, 0))
                                    matched_lines.add(line)
                                    includes.append({ tag_key: obj.attributes[tag_key] for tag_key in args.linked_tags if tag_key in obj.attributes })
                                    break
        #for line in sorted(matched_lines):
        #    print(line)
        #for inc in includes:
        #    print(inc)
        # search linked attributes
        with open(args.gff, 'r') as flo:
            for i, line in enumerate(flo):
                line = line.rstrip()
                obj = Feature.parse_gff_line(line)
                if obj:
                    for inc in includes:
                        for k, v in inc.items():
                            #errors = len(q)//4
                            #inc_pattern = '(?:'+v+'){e<='+str(errors)'}'
                            inc_pattern = v
                            if k in obj.attributes:
                                m = regex.search(inc_pattern, obj.attributes[k], flags=regex.IGNORECASE)
                                if m:
                                    #matched_lines.add((line, 1))
                                    matched_lines.add(line)
                                    break
        
        if args.header:
            print('# '+'\t'.join(['seqid', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']))
        
        for line in sorted(matched_lines):
            print(line)
        # End _feature()
    
    def _extract(self, args):
        """
        Search input FASTA headers for arbitrary text
        Print out headers+sequences that match
        """
        
        # Build the regex pattern
        q_patterns = []
        for q in args.query:
            if args.allow_errors:
                errors = len(q)//4
                q_pattern = '(?:'+q+'){e<='+str(errors)+'}'
            else:
                q_pattern = q
            q_patterns.append(q_pattern)
        
        pattern = '|'.join(q_patterns)
        
        print_on = False
        for fn in args.fasta:
            with open(fn, 'r') as flo:
                for line in flo:
                    line = line.rstrip()
                    if (len(line) > 0):
                        if line.startswith('>'):
                            m = regex.search(pattern, line, flags=regex.IGNORECASE | regex.ENHANCEMATCH) # regex.BESTMATCH
                            if m:
                                print_on = True
                            else:
                                print_on = False
                        if print_on:
                            print(line)
        
        # End _extract()
    
    def _evaluate(self, args):
        '''UNDER DEVELOPMENT'''
        print("Perform the evaluation here.")
        
        contig_sequences = utils.load_multiple_fasta_files(args.fasta)
        
        dDNA_file = ExcisionDonor.generate_fasta(os.path.join(args.folder, 'excision-dDNAs.fasta'))
        genome_fasta_file = utils.write_merged_fasta(contig_sequences, os.path.join(args.folder, 'genome.fasta'))
        
        # search spacer FASTA against genome+dDNA FASTA
        genome_index_file = args.selected_aligner.index(genome_fasta_file, os.path.basename(genome_fasta_file), args.folder, args.processors)
        dDNA_index_file = args.selected_aligner.index(dDNA_file, os.path.basename(dDNA_file), args.folder, args.processors)
        
        q2gDNA_align_file = args.selected_aligner.align(ex_query_file, genome_index_file, 'excision-query-2-gDNA.'+args.selected_aligner.output, args.folder, args.processors)
        q2exdDNA_align_file = args.selected_aligner.align(ex_query_file, ex_dDNA_index_file, 'excision-query-2-excision-dDNA.'+args.selected_aligner.output, args.folder, args.processors)
        
        ExcisionTarget.load_alignment(q2gDNA_align_file, args, contig_sequences)
        ExcisionTarget.load_alignment(q2exdDNA_align_file, args, ExcisionDonor.get_contig_dict())
        
        # Calculate off-target/guide scores for each algorithm
        logging.info("ExcisionTarget after SAM parsing and off-target scoring")
        ##### Add short-circuit/heuristic #####
        for et_seq, et_obj in ExcisionTarget.sequences.items():
            et_obj.score_off_targets(args, homologs)
            logging.info(et_obj)
            for a in et_obj.alignments:
                logging.info('  ' + str(a))
        
        # Batch calculate with new ExcisionTarget class
        ExcisionTarget.score_batch()
        
        logging.info("ExcisionTarget after Azimuth calculation")
        for et_seq, et_obj in ExcisionTarget.sequences.items():
            logging.info(et_obj)
        
    def get_gene_from_feature(self, feature_name, feature2gene):
        parent = Feature.features[feature_name].get_expand_parent().name
        
        return feature2gene[parent]
    
    def _generate(self, args):
        """Perform complete CRISPR/Cas analysis for input"""
        #contigs = utils.load_fasta_file(args.fasta[0]) # To do --> Do this for all FASTA files, then merge them into the same dictionary
        
        # Load the FASTA file specified on the command line
        # Merge all sequence information into the same dictionary
        #fasta_index, contig_index, contig_sequences = utils.load_indexed_fasta_files(args.fasta)
        contig_sequences = utils.load_multiple_fasta_files(args.fasta)
        
        # Open and parse the GFF file specified on the command line
        #features = utils.load_gff_file(args.gff, args.features, args.tag)
        # Filter features by what is selected
        #features = self.filter_features(features, args.selection)
        for gff_file in args.gff:
            Feature.load_gff_file(gff_file, args.features, args.excluded_features, args.selection, args.tag)
        Feature.assert_features(args.selection, contig_sequences)
        
        # Make index of homologs
        if args.homologs:
            homologs, feature2gene = utils.load_homologs(args.homologs)
        else:
            homologs, feature2gene = None, None
        
        logging.info('Feature.features')
        for f_name, f in sorted(Feature.features.items()):
            logging.info("  {}:{}:{}:{}..{}".format(f.name, f.contig, f.strand, f.start, f.end))
        logging.info('Feature.excluded_features')
        for exf_name, f in sorted(Feature.excluded_features.items()):
            logging.info("  {}:{}:{}:{}..{}".format(f.name, f.contig, f.strand, f.start, f.end))
        
        # Merge features?
        #features = merge_features(features)
        
        
        #### Some checks that need to be added ####
        # The feature being targeted MUST not go up to the edge of the contigs
        # Otherwise there will be no junction. i.e. the upstream and downstream
        # regions won't exist.
        # Actually, these regions must be at minimum, 50 nt
        ###########################################
        
        
        if (args.ko_gRNA): # Takes the values: (True/False)
            # Search for good targets within specified features
            # If no good target is found, then expand the feature
            pass
        
        # This code is performed in 'ExcisionDonor.generate_donors()'
        # if (args.ko_dDNA): # (mintag/addtag/unitag/bartag)
        #     if (args.ko_dDNA == 'mintag'):
        #         pass
        #     elif (args.ko_dDNA == 'addtag'):
        #         pass
        #     elif (args.ko_dDNA == 'unitag'):
        #         pass
        #     elif (args.ko_dDNA == 'bartag'):
        #         pass
        
        
        # Design gRNAs to target the ko-dDNA.
        if (args.ki_gRNA): # Takes the values: (True/False)
            pass
        
        
        
        
        # Code here (maybe as part of ExcisionTarget.search_all_features()), should expand features
        # if necessary. How?
        #   Look to see if the us/ds homology regions flanking the feature are identical.
        #     If they are identical, then this would be a homozygous dDNA
        #     If they are different, then this would be an allele-specific (heterozygous) dDNA
        
        if (args.ko_dDNA or args.ko_gRNA):
            # Expand features if necessary
            if (args.feature_expansion_method != None):
                Feature.expand_all_features(args, contig_sequences)
                
                # Print the set of new features
                logging.info('Feature.features')
                for f_name, f in sorted(Feature.features.items()):
                    logging.info("  {}:{}:{}:{}..{} PARENT={}".format(f.name, f.contig, f.strand, f.start, f.end, f.get_expand_parent().name))
        
        if (args.ko_gRNA):
            # Search features within contigs for targets that match the motifs
            # Old code (without feature expansion) ExcisionTarget.get_targets(args, contig_sequences, features)
            ExcisionTarget.search_all_features(args, contig_sequences)
            
            # Write the query list to FASTA
            ex_query_file = ExcisionTarget.generate_query_fasta(os.path.join(args.folder, 'excision-query.fasta'))
        
        # Generate excision dDNAs and their associated reversion gRNA spacers
        ExcisionDonor.generate_donors(args, contig_sequences)
        
        if (args.ki_gRNA):
            ReversionTarget.get_targets()
        ex_dDNA_file = ExcisionDonor.generate_fasta(os.path.join(args.folder, 'excision-dDNAs.fasta')) # Program will fail with error if this file is empty...
        
        if (args.ki_gRNA):
            re_query_file = ReversionTarget.generate_query_fasta(os.path.join(args.folder, 'reversion-query.fasta'))
        
        if (args.ki_dDNA != None): # args.ki_dDNA can take one of 3 possible values: (None/True/'*.fasta')
            # Generate reversion dDNAs and write them to FASTA
            ReversionDonor.generate_donors(args, contig_sequences, feature2gene)
            re_dDNA_file = ReversionDonor.generate_fasta(os.path.join(args.folder, 'reversion-dDNAs.fasta'))
        
        # Merge input FASTA files into a single one
        genome_fasta_file = utils.write_merged_fasta(contig_sequences, os.path.join(args.folder, 'genome.fasta'))
        
        # Index args.fasta for alignment
        #index_file = index_reference(args)
        genome_index_file = args.selected_aligner.index(genome_fasta_file, os.path.basename(genome_fasta_file), args.folder, args.processors)
        ex_dDNA_index_file = args.selected_aligner.index(ex_dDNA_file, os.path.basename(ex_dDNA_file), args.folder, args.processors)
        if (args.ki_dDNA == True):
            re_dDNA_index_file = args.selected_aligner.index(re_dDNA_file, os.path.basename(re_dDNA_file), args.folder, args.processors)
        
        if (args.ko_gRNA):
            # Use selected alignment program to find all matches in the genome and dDNAs
            #ex_genome_align_file = align(ex_query_file, genome_index_file, args)
            exq2gDNA_align_file = args.selected_aligner.align(ex_query_file, genome_index_file, 'excision-query-2-gDNA.'+args.selected_aligner.output, args.folder, args.processors)
            exq2exdDNA_align_file = args.selected_aligner.align(ex_query_file, ex_dDNA_index_file, 'excision-query-2-excision-dDNA.'+args.selected_aligner.output, args.folder, args.processors)
            
            #print("ExcisionTarget before SAM parsing")
            #for et_seq, et_obj in ExcisionTarget.sequences.items():
            #    print(et_obj)
            
            # Load the SAM files and add Alignments to ExcisionTarget sequences
            ExcisionTarget.load_alignment(exq2gDNA_align_file, args, contig_sequences)
            ExcisionTarget.load_alignment(exq2exdDNA_align_file, args, ExcisionDonor.get_contig_dict())
            
            # Calculate off-target/guide scores for each algorithm
            logging.info("ExcisionTarget after SAM parsing and off-target scoring")
            ##### Add short-circuit/heuristic #####
            for et_seq, et_obj in ExcisionTarget.sequences.items():
                et_obj.score_off_targets(args, homologs)
                logging.info(et_obj)
                for a in et_obj.alignments:
                    logging.info('  ' + str(a))
            
            # Batch calculate with new ExcisionTarget class
            ExcisionTarget.score_batch()
            
            logging.info("ExcisionTarget after Azimuth calculation")
            for et_seq, et_obj in ExcisionTarget.sequences.items():
                logging.info(et_obj)
        
        # Generate the FASTA with the final scores
        excision_spacers_file = ExcisionTarget.generate_spacers_fasta(os.path.join(args.folder, 'excision-spacers.fasta'))
        
        # Use selected alignment program to find all matches in the genome and dDNAs
        #re_align_file = align(re_query_file, genome_index_file, args)
        if (args.ki_gRNA):
            req2gDNA_align_file = args.selected_aligner.align(re_query_file, genome_index_file, 'reversion-query-2-gDNA.'+args.selected_aligner.output, args.folder, args.processors)
            req2exdDNA_align_file = args.selected_aligner.align(re_query_file, ex_dDNA_index_file, 'reversion-query-2-excision-dDNA.'+args.selected_aligner.output, args.folder, args.processors)
        if (args.ki_dDNA == True):
            req2redDNA_align_file = args.selected_aligner.align(re_query_file, re_dDNA_index_file, 'reversion-query-2-reversion-dDNA.'+args.selected_aligner.output, args.folder, args.processors)
        
        # Load the SAM files and add Alignments to ReversionTarget sequences
        if (args.ki_gRNA):
            ReversionTarget.load_alignment(req2gDNA_align_file, args, contig_sequences)
            ReversionTarget.load_alignment(req2exdDNA_align_file, args, ExcisionDonor.get_contig_dict())
        if (args.ki_dDNA == True):
            ReversionTarget.load_alignment(req2redDNA_align_file, args, ReversionDonor.get_contig_dict())
        
        
        if (args.ki_gRNA):
            # Calculate off-target/guide scores for each algorithm
            logging.info("ReversionTarget after SAM parsing and off-target scoring")
            # Somehow need to prioritize sequences for scoring.
            # Sequences with best diversity should be scored first.
            # If calculation time runs out, then just stop.
            # Should time each feature separately?
            ##### Subroutine #####
            #start_time = time.time()
            #if (time.time() - start_time >= args.max_time):
            #    logging.info('Site search terminated due to time constraints.')
            #    return
            ##### Subroutine #####
            for re_seq, re_obj in ReversionTarget.sequences.items():
                re_obj.score_off_targets(args, homologs)
                logging.info(re_obj)
                for a in re_obj.alignments:
                    logging.info('  ' + str(a))
            
            # Batch calculate with new ReversionTarget class
            ReversionTarget.score_batch()
            
            logging.info("ReversionTarget after Azimuth calculation")
            for rt_seq, rt_obj in ReversionTarget.sequences.items():
                logging.info(rt_obj)
            
            # Generate the FASTA with the final scores
            reversion_spacers_file = ReversionTarget.generate_spacers_fasta(os.path.join(args.folder, 'reversion-spacers.fasta'))
        
        # Test code to generate alignments
        ExcisionDonor.generate_alignments()
        
        # Pick out the best ones and print them out
        self.get_best_table(args, homologs, feature2gene)
        self.log_results(args, homologs, n=5)
        
        # Existing pipeline - Want to mutate feature, and cut site within feature
        # ,,,,,,,,,,,,,,,,NNNNNNNNNNNNNNNNNNNN111NNN................... Feature to be mutated contains cut site 1
        # ,,,,,,,,,,,,,,,,aaa222aaa.................................... Feature is knocked-out, and unique cut site 2 added
        # ,,,,,,,,,,,,,,,,mmmmmmmmmmmmmmmmmmmmmmmmmm................... Cut site 2 targeted and replaced by mutant
        
        # New pipeline - Want to mutate feature, but cut site outside of feature
        # ,,,,,,,,,,,,,,,,NNNNNNNNNNNNNN---------------111---.......... Feature to be mutated does not contain closest cut site 1
        # ,,,,,,,,,,,,,,,,nnnnnnnnnnnnnnnnnnnnnnnnnnnnn111nnn.......... Feature expanded to include cut site
        # ,,,,,,,,,,,,,,,,aaa222aaa.................................... Feature is knocked-out, and unique cut site 2 added
        # ,,,,,,,,,,,,,,,,mmmmmmmmmmmmmm---------------111---.......... Cut site 2 targeted and replaced by mutant (input) plus wild type cut site 1
        
        # 1) input: feature start and end
        # 2) if closest/best cut site is outside the feature
        #    Then expand the feature to include the closest/best cut site
        #    Also, concatenate the wild-type version of the expanded site to the mutant
        #    to create the new ki-dDNA
        
        # What if multiple features are adjacent to each other?
        # Then they should be combined into a single transformation.
        # ,,,,,,,,,NNNNNNNNNNN111NNN........NNNNNNNNN111NNN............ Since features 1 and 2 are independent, they can be kept as separate reactions
        # ,,,,,,,,,NNNNNNNNNN--------111-------NNNNNNNNNNNN............ Since features 1 and 2 both have the closest/best cut site, they are combined
        # ,,,,,,,,,nnnnnnnnnnnnnnnnnn111nnnnnnnnnnnnnnnnnnn............
        # ,,,,,,,,,aaa222aaa...........................................
        # ,,,,,,,,,mmmmmmmmmm--------111-------mmmmmmmmmmmm............ the dDNA is expanded to include the wt region in between both of these features
        
    def make_primer_set(self, args, qname, us_seq, ds_seq, insert_seq, feature_seq, q_hih_seq, s_ush_seq, q_ush_seq, s_dsh_seq, q_dsh_seq):
        tm_max_difference = 4.0
        amplicon_size = (200, 900)
        tm_range = (52, 64)
        primer_length_range = (19,32)
        min_delta_g = -5.0
        
        subset_size = 100 # 200 # Temporary size limit for the number of primers that go into the pair() function
        # primer group:
        #          feature: Fo      + Rf
        #          feature:      Ff +      Ro
        # optional-feature:      Ff + Rf
        #           insert: Fo      + Ri
        #           insert:      Fi +      Ro
        #           insert: Fo      +      Ro
        # optional-insert:       Fi + Ri
        
        # Save the sequence regions
        sseq_upstream = us_seq[-max(amplicon_size):]
        sseq_downstream = ds_seq[:max(amplicon_size)]
        sseq_insert = insert_seq
        sseq_feature = feature_seq
        
        # Make lists of primers for each region
        logging.info('Scanning the 4 regions (upstream, downstream, feature, insert) for 6 sets of primers')
        
        # try:
        #     logging.info('Loading primers from region+direction...')
        #     
        #     upstream_F   = load_object('upstream_F_'+qname+'.pickle')
        #     downstream_R = load_object('downstream_R_'+qname+'.pickle')
        #     feature_F    = load_object('feature_F_'+qname+'.pickle')
        #     feature_R    = load_object('feature_R_'+qname+'.pickle')
        #     insert_F     = load_object('insert_F_'+qname+'.pickle')
        #     insert_R     = load_object('insert_R_'+qname+'.pickle')
        #     
        # except FileNotFoundError:
        logging.info('Calculating primers found in each region+direction...')
        temp_folder = os.path.join('/dev/shm/addtag', os.path.basename(args.folder))
        
        logging.info('Number of primers found in each region+direction:')
        upstream_F   = args.selected_oligo.scan(sseq_upstream,   'left',  primer_size=primer_length_range, tm_range=tm_range, min_delta_g=min_delta_g, folder=temp_folder)
        logging.info('  len(upstream_F) = {}'.format(len(upstream_F)))
        
        downstream_R = args.selected_oligo.scan(sseq_downstream, 'right', primer_size=primer_length_range, tm_range=tm_range, min_delta_g=min_delta_g, folder=temp_folder)
        logging.info('  len(downstream_R) = {}'.format(len(downstream_R)))
        
        feature_F    = args.selected_oligo.scan(sseq_feature,    'left',  primer_size=primer_length_range, tm_range=tm_range, min_delta_g=min_delta_g, folder=temp_folder)
        logging.info('  len(feature_F) = {}'.format(len(feature_F)))
        
        feature_R    = args.selected_oligo.scan(sseq_feature,    'right', primer_size=primer_length_range, tm_range=tm_range, min_delta_g=min_delta_g, folder=temp_folder)
        logging.info('  len(feature_R) = {}'.format(len(feature_R)))
        
        insert_F     = args.selected_oligo.scan(sseq_insert,     'left',  primer_size=primer_length_range, tm_range=tm_range, min_delta_g=min_delta_g, folder=temp_folder)
        logging.info('  len(insert_F) = {}'.format(len(insert_F)))
        
        insert_R     = args.selected_oligo.scan(sseq_insert,     'right', primer_size=primer_length_range, tm_range=tm_range, min_delta_g=min_delta_g, folder=temp_folder)
        logging.info('  len(insert_R) = {}'.format(len(insert_R)))
        #     
        #     logging.info('Saving: primers found in each region+direction')
        #     save_object(upstream_F,   'upstream_F_'+qname+'.pickle')
        #     save_object(downstream_R, 'downstream_R_'+qname+'.pickle')
        #     save_object(feature_F,    'feature_F_'+qname+'.pickle')
        #     save_object(feature_R,    'feature_R_'+qname+'.pickle')
        #     save_object(insert_F,     'insert_F_'+qname+'.pickle')
        #     save_object(insert_R,     'insert_R_'+qname+'.pickle')
        
        # Sort each set of region+direction primers
        upstream_F = sorted(upstream_F, key=lambda x: x.weight, reverse=True)
        downstream_R = sorted(downstream_R, key=lambda x: x.weight, reverse=True)
        feature_F = sorted(feature_F, key=lambda x: x.weight, reverse=True)
        feature_R = sorted(feature_R, key=lambda x: x.weight, reverse=True)
        insert_F = sorted(insert_F, key=lambda x: x.weight, reverse=True)
        insert_R = sorted(insert_R, key=lambda x: x.weight, reverse=True)
        
        # Select the best subset of primers so the pairing doesn't take too long
        # This will cap the maximum number of primer pairs for each region span
        # to be 100 x 100 = 10,000
        upstream_F = upstream_F[:subset_size]
        downstream_R = downstream_R[:subset_size]
        feature_F = feature_F[:subset_size]
        feature_R = feature_R[:subset_size]
        insert_F = insert_F[:subset_size]
        insert_R = insert_R[:subset_size]
        logging.info('upstream_F: skipping {}/{} calculated primers'.format(max(0, len(upstream_F)-subset_size), len(upstream_F)))
        logging.info('downstream_R: skipping {}/{} calculated primers'.format(max(0, len(downstream_R)-subset_size), len(downstream_R)))
        logging.info('feature_F: skipping {}/{} calculated primers'.format(max(0, len(feature_F)-subset_size), len(feature_F)))
        logging.info('feature_R: skipping {}/{} calculated primers'.format(max(0, len(feature_R)-subset_size), len(feature_R)))
        logging.info('insert_F: skipping {}/{} calculated primers'.format(max(0, len(insert_F)-subset_size), len(insert_F)))
        logging.info('insert_R: skipping {}/{} calculated primers'.format(max(0, len(insert_R)-subset_size), len(insert_R)))
        
        # Pseudocode:
        #   create list to hold 6-set of primers, and their joint weight
        #   make all upstream_F downstream_R pairs, and weigh them
        #   starting with the best upstream_F and downstream_R pairs:
        #       pair upstream_F with all feature_R
        #       pair upstream_F with all insert_R
        #       pair all feature_F with downstream_R
        #       pair all insert_F with downstream_R
        #       calculate joint weight (product of all pair weights)
        #       add 6-set with joint weight to list
        #   select the best-weighted 6-set
        
        # Make all upstream_F downstream_R pairs, and then weigh them
        # Sort all primer pairs such that the ones with greatest weight are first in the list
        
        # try:
        #     logging.info('Loading: uf_dr_paired_primers')
        #     uf_dr_paired_primers = load_object('uf_dr_paired_primers_'+qname+'.pickle')
        # except FileNotFoundError:
        
        ###### Old version ######
        #logging.info('Calculating: uf_dr_paired_primers (amplicon across insert)...')
        #uf_dr_paired_primers = sorted(
        #    args.selected_oligo.pair(upstream_F, downstream_R, amplicon_size=amplicon_size, tm_max_difference=tm_max_difference, intervening=len(q_hih_seq), min_delta_g=min_delta_g, folder=temp_folder),
        #    key=lambda x: x.get_joint_weight(), reverse=True)
        #logging.info('         len(uf_dr_paired_primers) = {}'.format(len(uf_dr_paired_primers)))
        ###### End old version ######
        
        logging.info('Calculating: uf_dr_paired_primers...')
        uf_dr_paired_primers = sorted(
            args.selected_oligo.pair(upstream_F, downstream_R, amplicon_size=(amplicon_size[0]-len(s_ush_seq)-len(s_dsh_seq), amplicon_size[1]-len(s_ush_seq)-len(s_dsh_seq)), tm_max_difference=tm_max_difference, intervening=0, min_delta_g=min_delta_g, folder=temp_folder),
            key=lambda x: x.get_joint_weight(), reverse=True)
        logging.info('         len(uf_dr_paired_primers) = {}'.format(len(uf_dr_paired_primers)))
        #     
        #     logging.info('Saving: uf_dr_paired_primers')
        #     save_object(uf_dr_paired_primers, 'uf_dr_paired_primers_'+qname+'.pickle')
        # 
        # try:
        #     logging.info('Loading: uf_fr_paired_primers')
        #     uf_fr_paired_primers = load_object('uf_fr_paired_primers_'+qname+'.pickle')
        # except FileNotFoundError:
        logging.info('Calculating: uf_fr_paired_primers...')
        uf_fr_paired_primers = sorted(
            args.selected_oligo.pair(upstream_F, feature_R, amplicon_size=amplicon_size, tm_max_difference=tm_max_difference, intervening=len(s_ush_seq), min_delta_g=min_delta_g, folder=temp_folder),
            key=lambda x: x.get_joint_weight(), reverse=True)
        logging.info('         len(uf_fr_paired_primers) = {}'.format(len(uf_fr_paired_primers)))
        #     
        #     logging.info('Saving: uf_fr_paired_primers')
        #     save_object(uf_fr_paired_primers, 'uf_fr_paired_primers_'+qname+'.pickle')
        # 
        # try:
        #     logging.info('Loading: uf_ir_paired_primers')
        #     uf_ir_paired_primers = load_object('uf_ir_paired_primers_'+qname+'.pickle')
        # except:
        logging.info('Calculating: uf_ir_paired_primers...')
        uf_ir_paired_primers = sorted(
            args.selected_oligo.pair(upstream_F, insert_R, amplicon_size=amplicon_size, tm_max_difference=tm_max_difference, intervening=len(q_ush_seq), min_delta_g=min_delta_g, folder=temp_folder),
            key=lambda x: x.get_joint_weight(), reverse=True)
        logging.info('         len(uf_ir_paired_primers) = {}'.format(len(uf_ir_paired_primers)))
        #     
        #     logging.info('Saving: uf_ir_paired_primers')
        #     save_object(uf_ir_paired_primers, 'uf_ir_paired_primers_'+qname+'.pickle')
        # 
        # try:
        #     logging.info('Loading: ff_dr_paired_primers')
        #     ff_dr_paired_primers = load_object('ff_dr_paired_primers_'+qname+'.pickle')
        # except:
        logging.info('Calculating: ff_dr_paired_primers...')
        ff_dr_paired_primers = sorted(
            args.selected_oligo.pair(feature_F, downstream_R, amplicon_size=amplicon_size, tm_max_difference=tm_max_difference, intervening=len(s_dsh_seq), min_delta_g=min_delta_g, folder=temp_folder),
            key=lambda x: x.get_joint_weight(), reverse=True)
        logging.info('         len(ff_dr_paired_primers) = {}'.format(len(ff_dr_paired_primers)))
        #     
        #     logging.info('Saving: ff_dr_paired_primers')
        #     save_object(ff_dr_paired_primers, 'ff_dr_paired_primers_'+qname+'.pickle')
        # 
        # try:
        #     logging.info('Loading: if_dr_paired_primers')
        #     if_dr_paired_primers = load_object('if_dr_paired_primers_'+qname+'.pickle')
        # except:
        logging.info('Calculating: if_dr_paired_primers...')
        if_dr_paired_primers = sorted(
            args.selected_oligo.pair(insert_F, downstream_R, amplicon_size=amplicon_size, tm_max_difference=tm_max_difference, intervening=len(q_dsh_seq), min_delta_g=min_delta_g, folder=temp_folder),
            key=lambda x: x.get_joint_weight(), reverse=True)
        logging.info('         len(if_dr_paired_primers) = {}'.format(len(if_dr_paired_primers)))
        #     
        #     logging.info('Saving: if_dr_paired_primers')
        #     save_object(if_dr_paired_primers, 'if_dr_paired_primers_'+qname+'.pickle')
        
        def nsum(n):
            return n*(n+1)/2
        
        def rank_order(x, reverse=False, shift=0):
            return [y+shift for y in sorted(range(len(x)), key=x.__getitem__, reverse=reverse)]
        
        def rank(x, reverse=False, shift=0):
            return [z[0]+shift for z in sorted(enumerate(sorted(enumerate(x), key=lambda w: w[1], reverse=reverse)), key=lambda y: y[1][0])]
        
        def rank_dist(x):
            s = nsum(len(x))
            return [y/s for y in rank(x, reverse=False, shift=1)]
        
        def product(x):
            z = 1
            for y in x:
                z *= y
            return z
        
        def random_choices(population, weights, k=1):
            """
            Return a k sized list of population elements chosen with replacement.
            """
            weight_sum = sum(weights)
            choices = zip(population, weights)
            values = []
            for i in range(k):
                r = random.uniform(0, weight_sum)
                upto = 0
                for c, w in choices:
                    if upto + w >= r:
                        values.append(c)
                        break
                    upto += w
                else:
                    values.append(random.choice(population))
            return values
        
        def filter_primer_pairs(pairs, forward=None, reverse=None):
            ooo = []
            
            if (forward and reverse):
                for pp in pairs:
                    if ((pp.forward_primer.sequence == forward.sequence) and (pp.reverse_primer.sequence == reverse.sequence)):
                        ooo.append(pp)
            elif forward:
                for pp in pairs:
                    if (pp.forward_primer.sequence == forward.sequence):
                        ooo.append(pp)
            elif reverse:
                for pp in pairs:
                    if (pp.reverse_primer.sequence == reverse.sequence):
                        ooo.append(pp)
            else:
                ooo = pairs
            
            return ooo
        
        def random_primer_by_weight_old(pairs, forward=None, reverse=None):
            ooo = filter_primer_pairs(pairs, forward, reverse)
            
            if (len(ooo) == 0):
                return None
            elif ('choices' in random.__all__):
                return random.choices(ooo, [x.get_joint_weight() for x in ooo])[0]
            else:
                return random_choices(ooo, [x.get_joint_weight() for x in ooo])[0]
        
        def random_primer_by_weight(pairs):
            if (len(pairs) == 0):
                return None
            elif ('choices' in random.__all__):
                return random.choices(pairs, [x.get_joint_weight() for x in pairs])[0]
            else:
                return random_choices(pairs, [x.get_joint_weight() for x in pairs])[0]
        
        # Build a set of random, best ones
        starting_max = 1000
        starting_count = 0
        
        # Create lists to hold 5-set of primer pairs, and their joint weight
        starting_set = [] # initial set
        finished_set = [] # final set
        
        uf_dr_iterator = iter(uf_dr_paired_primers)
        
        if ((starting_max != None) and (len(uf_dr_paired_primers) > starting_max)):
            logging.info('uf_dr_paired_primers: skipping {}/{} calculated primer pairs'.format(max(0, len(uf_dr_paired_primers)-starting_max), len(uf_dr_paired_primers)))
        
        while (starting_count < starting_max):
            # Increment the count
            starting_count += 1
            
            # Create random 6-sets, with each pair's probability determined by its joint weight
            #uf_dr_pair = random_primer_by_weight(uf_dr_paired_primers) # comment this out to prevent random sampling of this one
            try:
                uf_dr_pair = next(uf_dr_iterator)
            except StopIteration:
                break
            
            sources = [
                filter_primer_pairs(uf_fr_paired_primers, forward=uf_dr_pair.forward_primer),
                filter_primer_pairs(uf_ir_paired_primers, forward=uf_dr_pair.forward_primer),
                filter_primer_pairs(ff_dr_paired_primers, reverse=uf_dr_pair.reverse_primer),
                filter_primer_pairs(if_dr_paired_primers, reverse=uf_dr_pair.reverse_primer)
            ]
            
            logging.info('loop {}: length of sources: '.format(starting_count) + str([len(x) for x in sources]))
            
            uf_fr_pair = random_primer_by_weight(sources[0])
            uf_ir_pair = random_primer_by_weight(sources[1])
            ff_dr_pair = random_primer_by_weight(sources[2])
            if_dr_pair = random_primer_by_weight(sources[3])
            
            #if None in [uf_dr_pair, uf_fr_pair, ff_dr_pair]:
            #    logging.info("'None' in [uf_dr_pair, uf_fr_pair, ff_dr_pair] --> skipping to next...")
            #    continue
            
            set4 = [uf_fr_pair, uf_ir_pair, ff_dr_pair, if_dr_pair]
            
            group = [
                uf_dr_pair.forward_primer,
                uf_dr_pair.reverse_primer,
                uf_fr_pair.reverse_primer if uf_fr_pair else None,
                uf_ir_pair.reverse_primer if uf_ir_pair else None,
                ff_dr_pair.forward_primer if ff_dr_pair else None,
                if_dr_pair.forward_primer if if_dr_pair else None
            ]
            
            group_weight = args.selected_oligo.group_weight(group)
            joint_weight = group_weight * product(x.get_joint_weight() for x in [uf_dr_pair]+set4 if x)
            starting_set.append((joint_weight, [uf_dr_pair]+set4))
            
            logging.info('loop {}: starting_set: '.format(starting_count) + str(starting_set[-1]))
            
            # iteratively improve until a local maxima is found
            # By continually swapping out the least-weighted component primer
            min_weight_delta = 0.0000000001
            weight_delta = 1
            while(weight_delta > min_weight_delta):
                # Select the component primer with the worst weight
                # excluding uf and dr
                swappables = [
                    set4[0].reverse_primer if set4[0] else None, # uf_fr_pair
                    set4[1].reverse_primer if set4[1] else None, # uf_ir_pair
                    set4[2].forward_primer if set4[2] else None, # ff_dr_pair
                    set4[3].forward_primer if set4[3] else None  # if_dr_pair
                ]
                #swappables = [s for s in swappables if s]
                
                # Get list of indices from smallest weight to largest weight
                wi_order = rank_order([x.weight if (x != None) else math.inf for x in swappables])
                #worst_index = wi_order[0]
                
                #worst = swappables[worst_index]
                #worst.weight
                
                # start at the worst, and go to the next-worst, then next, then next
                for wi in wi_order:
                    new_set4 = set4[:] #[uf_fr_pair, uf_ir_pair, ff_dr_pair, if_dr_pair]
                    
                    # Swap the worst-performing primer with its best alternative
                    #for source_pp in sources[worst_index]:
                    for source_pp in sources[wi]:
                        
                        #new_set4[worst_index] = source_pp
                        new_set4[wi] = source_pp
                        
                        new_group = [
                            uf_dr_pair.forward_primer,
                            uf_dr_pair.reverse_primer,
                            new_set4[0].reverse_primer if new_set4[0] else None,
                            new_set4[1].reverse_primer if new_set4[1] else None,
                            new_set4[2].forward_primer if new_set4[2] else None,
                            new_set4[3].forward_primer if new_set4[3] else None
                        ]
                        
                        new_group_weight = args.selected_oligo.group_weight(new_group)
                        new_joint_weight = new_group_weight * product(x.get_joint_weight() for x in [uf_dr_pair]+new_set4 if x)
                        
                        if (new_joint_weight > joint_weight):
                            logging.info('                 set: ' + str((new_joint_weight, [uf_dr_pair]+new_set4)))
                            
                            weight_delta = new_joint_weight - joint_weight
                            
                            # replace values for next iteration
                            set4 = new_set4
                            joint_weight = new_joint_weight
                            break
                    else:
                        weight_delta = 0
                    
                    if (weight_delta > 0):
                        break
            
            # Does not re-calculate 'pp.weight', thus the amplicon_size will not be further penalized
            # However, the 'pp.get_amplicon_size()' will accurately reflect the updated 'pp.intervening' length
            uf_insert_dr_pair = copy.deepcopy(uf_dr_pair)
            uf_insert_dr_pair.intervening = len(q_hih_seq)
            uf_feature_dr_pair = copy.deepcopy(uf_dr_pair)
            uf_feature_dr_pair.intervening = len(s_ush_seq)+len(sseq_feature)+len(s_dsh_seq)
            
            #finished_set.append((joint_weight, [uf_dr_pair]+set4)) # Originally, only a single copy of 'uf_dr_pair'
            finished_set.append((joint_weight, [uf_insert_dr_pair, uf_feature_dr_pair]+set4)) # two copies of 'uf_dr_pair'
        
        # Suzanne suggests
        # 1) start with best from sorted list of uf_dr_pairs
        # 2) run the optimization a finite number of times, or until the delta drops below a certain amount
        #    by swapping worst component
        #
        
        return starting_set, finished_set
    
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
    
    def filter_alignment_records_for_cPCR(self, args, alignment_filename, dDNA_contigs):
        """
        Parse the alignment to identify the exogenous DNA, as well as pairs of left/right flanking homology regions
        Note: This will only allow for each dDNA to align to a single place on each contig.
        If a single dDNA could target multiple loci per contig, then this function will not suffice
        Thus, this code needs modifications for broader applications.
        """
        # Allow for alignments to omit up to 3 nt at edges (of homology regions)
        permitted_edge = 3
        
        # Make a set of invalid records that will be removed because they failed some quality threshold
        invalid_records = set()
        
        # Make dict of all decent records
        records = {}
        
        # Read the alignment file
        with open(alignment_filename, 'r') as flo:
            record = None
            while True:
                record = args.selected_aligner.load_record(flo)
                if (record == None):
                    break
                else:
                    # Process the record
                    # Require certain e-value and length for a "significant" alignment
                    if ((record.evalue <= args.max_evalue) and (record.length >= args.min_length)):
                        # Require alignment to occur at the termini of the query (either at the extreme beginning, or extreme end)
                        #if ((record.query_position[0] < permitted_edge) or (record.query_position[1] > len(dDNA_contigs[record.query_name])-permitted_edge)):
                        #    records.setdefault((record.query_name, record.subject_name), []).append(record)
                        qs_pair = (record.query_name, record.subject_name)
                        
                        # If the alignment is on the edges of the query
                        if ((record.query_position[0] < permitted_edge) or (record.query_position[1] > len(dDNA_contigs[record.query_name])-permitted_edge)):
                            if (qs_pair not in records):
                                records[qs_pair] = [None, None]
                        
                        # Arrange the records according to their flanking: records[(qname, sname)] = [left_record, right_record]
                        if (record.query_position[0] < permitted_edge): # left
                            if (records[qs_pair][0] == None):
                                records[qs_pair][0] = record
                            else:
                                logging.info(str(qs_pair) + ' has too many valid alignments')
                                invalid_records.add(qs_pair)
                        
                        elif (record.query_position[1] > len(dDNA_contigs[record.query_name])-permitted_edge): # right
                            if (records[qs_pair][1] == None):
                                records[qs_pair][1] = record
                            else:
                                logging.info(str(qs_pair) + ' has too many valid alignments')
                                invalid_records.add(qs_pair)
        
        logging.info("before filtering:")
        for kkk, vvv in records.items():
            logging.info("  "+str(kkk) + " " + str(vvv))
        
        # Queue for removal the (query, subject) pairs that only have 1 of the 2 required alignments
        for qs_pair, record_pair in records.items():
            if (None in record_pair):
                invalid_records.add(qs_pair)
        
        # Remove the (query, subject) pairs that had too many valid alignments
        for bad_key in invalid_records:
            records.pop(bad_key)
            logging.info('Removing invalid record: ' + str(bad_key))
        
        logging.info("after filtering:")
        for kkk, vvv in records.items():
            logging.info("  "+str(kkk) + " " + str(vvv))
        
        return records
    
    def filter_alignment_records2(self, args, alignment_filename, dDNA_contigs):
        records = {}
        
        # Allow for alignments to omit up to 3 nt at edges (of homology regions)
        permitted_edge = 3
        
        # Make a set of invalid records that will be removed because they failed some quality threshold
        invalid_records = set()
        
        # Read the alignment file
        with open(alignment_filename, 'r') as flo:
            record = None
            while True:
                record = args.selected_aligner.load_record(flo)
                if (record == None):
                    break
                else:
                    # Process the record
                    # Require certain e-value and length for a "significant" alignment
                    if ((record.evalue <= args.max_evalue) and (record.length >= args.min_length)):
                        # Require alignment to occur at the termini of the query (either at the extreme beginning, or extreme end)
                        #if ((record.query_position[0] < permitted_edge) or (record.query_position[1] > len(dDNA_contigs[record.query_name])-permitted_edge)):
                        #    records.setdefault((record.query_name, record.subject_name), []).append(record)
                        qs_pair = (record.query_name, record.subject_name)
                        
                        # The entire query should align to the subject
                        if (record.query_position[0] < permitted_edge < len(dDNA_contigs[record.query_name])-permitted_edge < record.query_position[1]):
                            if (qs_pair not in records):
                                records[qs_pair] = record
                            else:
                                logging.info(str(qs_pair) + ' has too many valid alignments')
                                invalid_records.add(qs_pair)
        for bad_qs_pair in invalid_records:
            records.pop(bad_qs_pair)
            logging.info('Removing invalid record: ' + str(bad_qs_pair))
        
        return records
        
    
    def _confirm(self, args):
        print("# Design confirmation primers here.")
        
        
        # Set intended amplicon size range
        #amplicon_size = (400, 700) # definition moved to make_primer_set() function
        
        # Hard-coded melting temp range (should be adjustable via command line in the future)
        #tm_range = (53, 57) # definition moved to make_primer_set() function
        
        # Require primer pairs to be within 2 degrees Celcius from each other
        #tm_max_difference = 3.0 # definition moved to make_primer_set() function
            
        # Order in which primer lengths should be interrogated, starting from the left
        # (18 doesn't work reliably with UNAFold)
        #primer_sizes = [20, 21, 19, 22, 23, 18, 24, 25, 26]
        
        #args.number_pcr_conditions
        
        # Define variables to hold slines for output
        output1 = []
        output2 = []
        
        # First, create the engineered genomes:
        #   FILENAME          DESCRIPTION
        #   genome-r0.fasta   wild type genomic DNA
        #   genome-r1.fasta   genomic DNA with dDNA1 incorporated
        #   genome-r2.fasta   genomic DNA with dDNA1 incorporated, then dDNA2 incorporated
        
        # For each round, create certain data structures or files
        genome_contigs_list = []
        dDNA_contigs_list = []
        dDNA_fasta_file_list = []
        dDNA_alignment_file_list = []
        genome_fasta_file_list = []
        
        #Datum = namedtuple('Datum', ['r', 'qname', 'sname', 'us_seq', 'ds_seq', 'insert_seq', 'feature_seq', 'q_hih_seq', 'q_ush_seq', 's_ush_seq', 'q_dsh_seq', 's_dsh_seq'])
        #Datum = namedtuple('Datum', ['r', 'sname', 'ush_start', 'ush_end', 'dsh_start', 'dsh_end'])
        Datum = namedtuple('Datum', ['dDNA_r', 'dDNA_contig', 'genome_r', 'genome_contig', 'ush_start', 'ush_end', 'dsh_start', 'dsh_end', 'ins_start', 'ins_end'])
        
        # Create 'r0'
        logging.info('Working on round r{}'.format(0))
        # Parse the input FASTA
        genome_contigs_list.append(utils.load_multiple_fasta_files(args.fasta))
        dDNA_contigs_list.append(None) # The 'r0' round of genome engineering has no dDNA
        dDNA_fasta_file_list.append(None)
        dDNA_alignment_file_list.append(None)
        
        # Merge input FASTA files into a single one
        genome_fasta_file_list.append(utils.write_merged_fasta(genome_contigs_list[0], os.path.join(args.folder, 'genome-r0.fasta')))
        
        # We need the groups of contig names that correspond to the loci in each round
        # For example:
        #   contig_groups = [ ['chr1', 'chr1-r1[gene1A,gene2]', 'chr1-r1[gene1A,gene2]-r2[gene1B]'] ]
        # Populate contig_groups with the initial, 'r0' contig names
        contig_groups = [[c] for c in genome_contigs_list[0]]
        
        # We link the dDNA contigs together
        # For example:
        #   dDNA_groups = [ ['gene1A', 'gene1B'], ['gene2'] ]
        dDNA_groups = []
        
        # We link the actual genome modifications together
        # For example:
        #   mod_groups = [
        #     [('chr1', 1022, 2000, 2300, 4500), ('chr1-r1[gene1A,gene2]', 1322, 2300, 2600, 4800), ('chr1-r1[gene1A,gene2]-r2[gene1B]', 1055, 2000, 2301, 4000)],
        #     [('chr1', 200, 220, 240, 300), ('chr1-r1[gene1A,gene2]', 200, 240, 240, 280)]
        #  ]
        mod_groups = []
        
        # group_links = {
        #     # (dDNA round, dDNA contig name): gDNA r0                                          gDNA r1
        #     (1, 'gene1A'):                    ((0, 'chr1', 100, 200, 300, 400,),              (1, 'chr1-r1[gene1A,gene2]', 200,300,400,500)),
        #     (2, 'gene1B'):                    ((1, 'chr1-r1[gene1A,gene2]', 200,300,400,500), (2, 'chr1-r1[gene1A,gene2]-r2[gene1B]', 100,200,300,400))
        # }
        # Logic: Because the first tuple of 'gene 1B' overlaps with the second tuple of 'gene1A', these are linked
        group_links = [] # {}
        
        # Record calculated dDNA homology regions
        #dDNA_homology_seqs = {} # key=(dDNA_r, dDNA_contig, genome_contig), value=(s_ush_seq, s_dsh_seq)
        
        # For each dDNA, do some operations
        for j, dDNA_filename in enumerate(args.dDNAs):
            # Variable to count the round of genome engineering
            r = j+1
            
            # Let the user know which round is being processed
            logging.info('Working on round r{}'.format(r))
            
            # Load input dDNA FASTA into the list of dicts
            dDNA_contigs_list.append(utils.old_load_fasta_file(dDNA_filename))
            
            # Build an index of the previous round's genome
            genome_index_file = args.selected_aligner.index(genome_fasta_file_list[r-1], os.path.basename(genome_fasta_file_list[r-1]), args.folder, args.processors)
            
            # Write current dDNA to FASTA
            dDNA_fasta_file_list.append(utils.write_merged_fasta(dDNA_contigs_list[r], os.path.join(args.folder, 'dDNA-r'+str(r)+'.fasta')))
            
            # We align dDNA against the genome+dDNA(prev)
            dDNA_alignment_file_list.append(args.selected_aligner.align(dDNA_fasta_file_list[r], genome_index_file, 'dDNA-r'+str(r)+'-alignment.'+args.selected_aligner.output, args.folder, args.processors))
            
            # Parse the alignment to make an intermediary record data structure
            records = self.filter_alignment_records_for_cPCR(args, dDNA_alignment_file_list[r], dDNA_contigs_list[r])
            
            # Iterate through all detected paired cross-over events (i.e. only simulate the "proper" cross-overs for significant alignment pairs)
            # Starting with the right-most position in each sequence, make changes to the genome contig sequence
            # We assume no record pairs overlap
            
            # Keep track of how contig names will change
            name_changes = {}
            
            # Seed the results of this round's engineering with the contigs from the previous round
            # (Will still need to change their names)
            genome_contigs_list.append(copy.deepcopy(genome_contigs_list[r-1]))
            
            # Sort the 'qs_pair' tuples from greatest to least based on position in each chromosome
            record_order = sorted(records, key=lambda x: records[x][0].subject_position, reverse=True)
            
            # Iterate through the sorted records
            for i, qs_pair in enumerate(record_order):
                qname, sname = qs_pair
                record_list = records[qs_pair]
                # Each should have only one upstream and one downstream homology region
                # Should be able to comment this condition out as it is already accounted for in the record population/filtering step
                if (len(record_list) == 2):
                    if (record_list[0].flags & 16 == record_list[1].flags & 16): # Each homology region should have identical strand orientation
                        logging.info('Working on {} vs {}'.format(qname, sname))
                        logging.info("crossover: " + str((i, qname, sname, record_list)))
                        
                        # Identify which record is upstream, and which is downstream (by query)
                        us_record, ds_record = sorted(record_list, key=lambda x: x.query_position)
                        
                        #new_name = sname+'-r'+str(r) # Add round number to the contig header
                        name_changes.setdefault(sname, []).append(qname)
                        
                        # get the entire contig sequences
                        qcontig = dDNA_contigs_list[r][qname]
                        scontig = genome_contigs_list[r][sname]
                        
                        # Get the QUERY homology coordinates
                        qush_start, qush_end = us_record.query_position[0], us_record.query_position[1]
                        qdsh_start, qdsh_end = ds_record.query_position[0], ds_record.query_position[1]
                        
                        # Get the SUBJECT homology coordinates
                        if (us_record.flags & 16): # Reverse complement means upstream/downstream are swapped
                            logging.info('  us_record is reverse complemented')
                            sush_start, sush_end = ds_record.subject_position[0], ds_record.subject_position[1]
                            sdsh_start, sdsh_end = us_record.subject_position[0], us_record.subject_position[1]
                        else: # not reverse-complemented
                            sush_start, sush_end = us_record.subject_position[0], us_record.subject_position[1]
                            sdsh_start, sdsh_end = ds_record.subject_position[0], ds_record.subject_position[1]
                        
                        
                        
                        ######## Get the sequences ########
                        
                        for seqloop in range(2):
                            ######## Manage QUERY sequences ########
                            
                            # Calculate new insert sequence (with homology arms) from dDNA
                            q_hih_seq = qcontig[qush_start:qdsh_end] # us-homology + insert + ds-homology (should be entire contig)
                            
                            # Calculate new insert sequence (excluding homology arms) from dDNA
                            insert_seq = qcontig[qush_end:qdsh_start] # insert
                            
                            # Get the up/down-stream homology regions for the dDNA (not the chromosome)
                            q_ush_seq = qcontig[qush_start:qush_end]
                            q_dsh_seq = qcontig[qdsh_start:qdsh_end]
                            
                            if (us_record.flags & 16): # Reverse complement if necessary
                                # We must reverse-complement these dDNA substrings so their orientation matches the genome
                                q_hih_seq = nucleotides.rc(q_hih_seq)
                                insert_seq = nucleotides.rc(insert_seq)
                                q_ush_seq, q_dsh_seq = nucleotides.rc(q_dsh_seq), nucleotides.rc(q_ush_seq)
                            
                            ######## Manage SUBJECT sequences ########
                            
                            # Extract the feature sequence (excluding homology arms)
                            feature_seq = scontig[sush_end:sdsh_start]
                            
                            # Calculate new upstream & downstream chromosome arms (swap us/ds for subject) (excludes homology arms)
                            us_seq = scontig[:sush_start]
                            ds_seq = scontig[sdsh_end:]
                            
                            # Get the up/down-stream homology regions for the chromosome (not the dDNA)
                            s_ush_seq = scontig[sush_start:sush_end]
                            s_dsh_seq = scontig[sdsh_start:sdsh_end]
                            
                            if (seqloop == 0):
                                # Overview for making homology regions mutually exclusive
                                #   No overlap
                                #     q          ......111111111111......2222222222.....           none
                                #     s      ....11111111111................222222222222.......    none
                                #   query overlap
                                #     q          ........1111111111XXX22222222.......              overlap
                                #     s    .......11111111111111............222222222222.......    none
                                #      q2        ........111111111122222222222.......
                                #      s2  .......11111111111'''............222222222222.......
                                #   subject overlap
                                #     q          ......111111111111......2222222222.....           none
                                #     s    ........11111111111111XXXXX2222222222222222....         overlap
                                #      q2        ......1111111'''''......2222222222.....  
                                #      s2  ........11111111111111222222222222222222222....
                                #   query & subject overlap
                                #     q          ........1111111111XXX22222222.......              overlap
                                #     s    ........11111111111111XXXXX2222222222222222....         overlap
                                #      q2        ........111111111122222222222.......
                                #      s2  ........11111111111111222222222222222222222....
                                
                                # if both query and subject overlap
                                if ((qdsh_start < qush_end) and (sdsh_start < sush_end)):
                                    # No fancy calculations required. Just set the new coordinates
                                    qush_end = qdsh_start
                                    sush_end = sdsh_start
                                # If query overlaps
                                elif (qdsh_start < qush_end):
                                    # Calculate the subject trim
                                    q_overlap = qcontig[qdsh_start:qush_end] # Just the overlapping subsequence of the query alignments
                                    if (us_record.flags & 16):
                                        q_overlap = nucleotides.rc(q_overlap)
                                    m = regex.search('(?:'+q_overlap+'){e<'+str(len(q_overlap)//2)+'}$', s_ush_seq, flags=regex.IGNORECASE|regex.BESTMATCH)
                                    if m:
                                        sush_end = sush_start + m.start()
                                    
                                    # Set the query upstream homology end
                                    qush_end = qdsh_start
                                # If subject overlaps
                                elif (sdsh_start < sush_end):
                                    # Calculate the query trim
                                    s_overlap = scontig[sdsh_start:sush_end]
                                    m = regex.search('(?:'+s_overlap+'){e<'+str(len(s_overlap)//2)+'}$', q_ush_seq, flags=regex.IGNORECASE|regex.BESTMATCH)
                                    if m:
                                        qush_end = qush_start + m.start()
                                    
                                    # Set the subject upstream homology end
                                    sush_end = sdsh_start
                                # If there is no overlap
                                else:
                                    # Exit this for loop because no coordinates were modified
                                    # and there is no need to re-calculate the sequences
                                    break
                                
                                # Otherwise, re-calculate all the sequences
                        
                        # Record the dDNA homology regions for later reference
                        # Data structure limitation that each dDNA can only engineer a single locus per contig (i.e. this code needs improvement)
                        #dDNA_homology_seqs[(r, qname, sname)] = (s_ush_seq, s_dsh_seq)
                        
                        #########
                        # Replace entry in 'genome_contigs' dict with the new cross-over contig
                        # This probably doesn't work with multiple targeted engineering events on the same chromosome
                        #genome_contigs_list[r][sname] = us_seq + q_hih_seq + ds_seq # <---------------------- this should be done in the 'r' outer loop (earlier draft)
                        
                        # Keep the subject up/down-stream homology regions in case 2 engineered loci have overlapping homology arms
                        # like this:
                        #    site1            site2
                        #   ...-----iii------...........
                        #   .............------iii---...
                        # The result would be:
                        #   ...-----iii--------iii---...
                        #genome_contigs_list[r][sname] = us_seq + s_ush_seq + insert_seq + s_dsh_seq + ds_seq # For some reason, this code doesn't properly incorporate dDNA into the genome
                        genome_contigs_list[r][sname] = us_seq + q_hih_seq + ds_seq # Giving this a try
                        
                        # this messes with calling self.calculate_amplicons() for each query-subject pair.
                        
                        ##### alternate #####
                        # This code only works for a single-locus per chromosome! <------ need to update this so it works for multiple loci per chromosome
                        #updated_contigs[sname] = us_seq + q_hih_seq + ds_seq
                        ### end alternate ###
                        
                        # Add this parsed record data to the 'group_links'
                        # Key refers to dDNA file and contig name, value referes to gDNA file, contig name, and homology regions
                        #group_links.append(Datum(r, qname, r-1, sname, sush_start, sush_end, sdsh_start, sdsh_end, us_record.query_position[1], ds_record.query_position[0]))
                        group_links.append(Datum(r, qname, r-1, sname, sush_start, sush_end, sdsh_start, sdsh_end, qush_end, qdsh_start))
                        
                else:
                    logging.info(qname + ' vs ' + sname + ' has ' + str(len(record_list)) + ' regions of homology (2 needed)')
            
            # Rename the contigs based on the modifications
            # This code only works for a single-locus per chromosome! <------ need to update this so it works for multiple loci per chromosome
            for sname in name_changes:
                new_name = sname+'-r'+str(r)+'['+','.join(name_changes[sname])+']'
                genome_contigs_list[r][new_name] = genome_contigs_list[r].pop(sname)
                ##### alternate #####
                #genome_contigs.pop(sname)
                #genome_contigs[new_name] = updated_contigs[sname]
                ### end alternate ###
                
                # Add the new contig name to the contig name group
                for cg in contig_groups:
                    if sname in cg:
                        cg.append(new_name)
                        break
                else: # This 'else' statement should never run because we've already seeded 'contig_groups' with 'r0' contig names
                    contig_groups.append([sname, new_name])
                
            
            # Write the new modified genome to FASTA file
            genome_fasta_file_list.append(utils.write_merged_fasta(genome_contigs_list[r], os.path.join(args.folder, 'genome-r'+str(r)+'.fasta')))
        
        # Now that the modified genomes have been calculated, we need to look for shared DNA outside of the homology regions
        
        
        
        
        # This aligns each dDNA-rN against the engineered genome-rN[].
        # That is, it tells us the new location of each dDNA contig in the new genome
        # This is important because multiple edits (insertions/deletions) may have occurred
        # in the same contig, and the locations may have shifted.
        # We loop through the dDNAs as before
        for j, dDNA_filename in enumerate(args.dDNAs):
            # Variable to count the round of genome engineering
            r = j+1
            
            # Let the user know which round is being processed
            logging.info('Working on round r{}'.format(r))
            
            # Align each dDNA to the genome it produced
            genome_index_file = args.selected_aligner.index(genome_fasta_file_list[r], os.path.basename(genome_fasta_file_list[r]), args.folder, args.processors)
            temp_alignment = args.selected_aligner.align(dDNA_fasta_file_list[r], genome_index_file, 'dDNA-temp-r'+str(r)+'-alignment.'+args.selected_aligner.output, args.folder, args.processors)
            
            # Parse the alignment to make an intermediary record data structure
            records = self.filter_alignment_records2(args, temp_alignment, dDNA_contigs_list[r])
            
            # Find the exact locations of each of the dDNAs within the engineered genomes
            # Iterate through the sorted records
            for i, qs_pair in enumerate(records):
                qname, sname = qs_pair
                record = records[qs_pair]
                
                logging.info('Working on {} vs {}'.format(qname, sname))
                
                scontig = genome_contigs_list[r][sname]
                qcontig = dDNA_contigs_list[r][qname]
                
                
                # Calculate new insert sequence (with homology arms) from dDNA
                q_hih_seq = qcontig[record.query_position[0]:record.query_position[1]] # us-homology + insert + ds-homology (should be entire contig)
                
                # Calculate new upstream & downstream chromosome arms (excludes homology arms)
                us_seq = scontig[:record.subject_position[0]]
                ds_seq = scontig[record.subject_position[1]:]
                
                s_hih_seq = scontig[record.subject_position[0]:record.subject_position[1]]
                
                sush_start, sush_end = record.subject_position[0], None
                sdsh_start, sdsh_end = None, record.subject_position[1]
                ins_start, ins_end = None, None
                
                if (record.flags & 16): # Reverse complement if necessary
                    logging.info('  us_record is reverse complemented')
                    # We must reverse-complement these dDNA substrings so their orientation matches the genome
                    q_hih_seq = nucleotides.rc(q_hih_seq)
                else: # not reverse-complemented
                    pass
            
                #group_links.setdefault((r, qname), []).append(Datum(r, sname, sush_start, sush_end, sdsh_start, sdsh_end))
                group_links.append(Datum(r, qname, r, sname, sush_start, sush_end, sdsh_start, sdsh_end, ins_start, ins_end))
        
        ######## Begin for linking the loci ########
        
        logging.info('contig_groups:')
        for cg in contig_groups:
            logging.info('  ' + str(cg))
        
        def in_same_contig_group(contig1, contig2, cgroups=contig_groups):
            for g in cgroups:
                if ((contig1 in g) and (contig2 in g)):
                    return True
            else:
                return False
        
        # Make a copy of 'group_links' that we can pop stuff from
        datum_list = copy.deepcopy(group_links)
        
        logging.info('datum_list:')
        for datum in datum_list:
            logging.info('  ' + str(datum))
        
        
        # For each (round, locus), we need the (contig name, upstream_homology_length, downstream_homology_length)
        # loop through these:
        
        # For each cross-over event, there will be a datum
        # Typically, there should be 1 Datum object for haploid genomes,
        # and 2 Datum objects for diploid genomes
        
        
        # Populate initial groups based on if there is crossing over in the genome-r0
        datum_groups = []
        to_pop = []
        for di, datum in enumerate(datum_list):
            if datum.genome_contig in genome_contigs_list[0].keys():
                datum_groups.append([datum])
                to_pop.append(di)
        for di in sorted(to_pop, reverse=True):
            datum_list.pop(di)
        to_pop = []
        
        logging.info('datum_groups (r0):')
        for dg in datum_groups:
            logging.info('  ' + str(dg))
        
        logging.info('datum_list:')
        for datum in datum_list:
            logging.info('  ' + str(datum))
        
        # Populate the intermediary groups
        # New try:
        # If there is ANY OVERLAP AT ALL between
        #      round n,    AFTER: ush_start-dsh_end
        # and  round n+1, BEFORE: ush_start-dsh_end
        # Then these loci should be linked!
        no_more_after = False
        #no_more_before = False
        groups = []
        wcount = 0
        while (not no_more_after):
            # Specify an error if for some reason this doesn't work
            wcount += 1
            if (wcount > 10000):
                raise Exception('Too many iterations taken to link loci between engineered genomes.')
            
            for target_genome_round in range(1, len(genome_contigs_list)-1):
                to_pop = []
                # Select the first datum ("after")
                datum_after = None
                for di, datum in enumerate(datum_list):
                    if (datum.genome_r == target_genome_round):
                        # If the rounds are the same then the alignment is "after engineering"
                        if (datum.dDNA_r == datum.genome_r):
                            datum_after = datum
                            to_pop.append(di)
                            break
                else:
                    no_more_after = True
                
                # Select the second datum ("before")
                datum_before = None
                if (datum_after != None):
                    for di, datum in enumerate(datum_list):
                        if (datum.genome_r == target_genome_round):
                            # If the rounds are different then the alignment is "before engineering"
                            if (datum.dDNA_r != datum.genome_r):
                                if ((datum.genome_contig == datum_after.genome_contig) and
                                    (datum.dDNA_contig != datum_after.dDNA_contig)):
                                    datum_before = datum
                                    to_pop.append(di)
                                    break
                    #else:
                    #    no_more_before = True
                
                else:
                    # If there is no "before" to pair with selected "after"
                    # then this is a terminal genome engineering
                    # Not sure if this code would work if a locus was "skipped" between rounds as follows:
                    #    genome-r0: locus A wt
                    #    genome-r1: locus A ko
                    #    genome-r2: locus A ko (no change)
                    #    genome-r3: locus A ki
                    
                    # Since this is a terminal, then we need to select the correct "before" using homology
                    pass
                if ((datum_after != None) and (datum_before != None)):
                    # Now that the two datum are selected, determine whether or not they overlap
                    # Feature.overlap_coverage(start1, end1, start2, end2, index_base=0)
                    overlap = Feature.overlap_coverage(datum_after.ush_start, datum_after.dsh_end, datum_before.ush_start, datum_before.dsh_end)
                    
                    # If they overlap
                    if (overlap > 0):
                        # Then add them to the same group
                        group = [datum_after, datum_before]
                        groups.append(group)
                        
                        logging.info('group:')
                        for g in group:
                            logging.info('  ' + str(g))
                        
                        for di in sorted(to_pop, reverse=True):
                            datum_list.pop(di)
                    # Otherwise, they are non-overlapping engineering events on the same genome contig
                    else:
                        # Do nothing
                        pass
        
        logging.info('datum_list:')
        for datum in datum_list:
            logging.info('  ' + str(datum))
        
        # Assuming we have 'group' with the [datum_after, datum_before]
        # Now we match homology regions with r0
        for dg in datum_groups:
            # group = [
            #     Datum(dDNA_r=1, dDNA_contig='exDonor-55', genome_r=1, genome_contig='Ca22chr3A_C_albicans_SC5314-r1[exDonor-55]', ush_start=943366, ush_end=None, dsh_start=None, dsh_end=943468)
            #     Datum(dDNA_r=2, dDNA_contig='reDonor-0', genome_r=1, genome_contig='Ca22chr3A_C_albicans_SC5314-r1[exDonor-55]', ush_start=942963, ush_end=943420, dsh_start=943414, dsh_end=943996)
            # ]
            ush_len = abs(dg[-1].ush_end - dg[-1].ush_start)
            dsh_len = abs(dg[-1].dsh_end - dg[-1].dsh_start)
            
            dg_ush_seq = genome_contigs_list[dg[-1].genome_r][dg[-1].genome_contig][dg[-1].ush_start:dg[-1].ush_start+ush_len]
            dg_dsh_seq = genome_contigs_list[dg[-1].genome_r][dg[-1].genome_contig][dg[-1].dsh_end-dsh_len:dg[-1].dsh_end]
            
            for group in groups:
                g_ush_seq = genome_contigs_list[group[0].genome_r][group[0].genome_contig][group[0].ush_start:group[0].ush_start+ush_len]
                g_dsh_seq = genome_contigs_list[group[0].genome_r][group[0].genome_contig][group[0].dsh_end-dsh_len:group[0].dsh_end]
                
                # If the rounds match, the dDNA_contigs match, and the genome_contigs match
                # and the homology regions match
                if ((dg[-1].dDNA_r == group[0].dDNA_r) and
                    (dg[-1].dDNA_contig == group[0].dDNA_contig) and
                    in_same_contig_group(group[0].genome_contig, dg[-1].genome_contig) and
                    (nucleotides.lcs(dg_ush_seq, g_ush_seq).size > 0.9*ush_len) and
                    (nucleotides.lcs(dg_dsh_seq, g_dsh_seq).size > 0.9*dsh_len)):
                    # Then populate their missing fields
                    #if (group[0].ush_end == None):
                    #    group[0].ush_end = group[-1].ush_start+ush_len
                    #if (group[0].dsh_start == None):
                    #    group[0].dsh_start = group[-1].dsh_end-dsh_len
                    
                    # Then add them to the same group
                    for g in group:
                        dg.append(g)
        
        logging.info('datum_groups:')
        for dg in datum_groups:
            logging.info('  ' + str(dg))
        
        # Populate the final/terminal group 'rN'
        # Using homology! (like we did with r0 above)
        for dg in datum_groups:
            # group = [
            #     Datum(dDNA_r=1, dDNA_contig='exDonor-55', genome_r=1, genome_contig='Ca22chr3A_C_albicans_SC5314-r1[exDonor-55]', ush_start=943366, ush_end=None, dsh_start=None, dsh_end=943468)
            #     Datum(dDNA_r=2, dDNA_contig='reDonor-0', genome_r=1, genome_contig='Ca22chr3A_C_albicans_SC5314-r1[exDonor-55]', ush_start=942963, ush_end=943420, dsh_start=943414, dsh_end=943996)
            # ]
            ush_len = abs(dg[-1].ush_end - dg[-1].ush_start)
            dsh_len = abs(dg[-1].dsh_end - dg[-1].dsh_start)
            
            dg_ush_seq = genome_contigs_list[dg[-1].genome_r][dg[-1].genome_contig][dg[-1].ush_start:dg[-1].ush_start+ush_len]
            dg_dsh_seq = genome_contigs_list[dg[-1].genome_r][dg[-1].genome_contig][dg[-1].dsh_end-dsh_len:dg[-1].dsh_end]
            
            to_pop = []
            for di, datum in enumerate(datum_list):
                d_ush_seq = genome_contigs_list[datum.genome_r][datum.genome_contig][datum.ush_start:datum.ush_start+ush_len]
                d_dsh_seq = genome_contigs_list[datum.genome_r][datum.genome_contig][datum.dsh_end-dsh_len:datum.dsh_end]
                
                # If the rounds match, the dDNA_contigs match, and the genome_contigs match
                # and the homology regions match
                if ((dg[-1].dDNA_r == datum.dDNA_r) and
                    (dg[-1].dDNA_contig == datum.dDNA_contig) and
                    in_same_contig_group(datum.genome_contig, dg[-1].genome_contig) and
                    #group[0].genome_contig.startswith(dg[-1].genome_contig) and
                    (nucleotides.lcs(dg_ush_seq, d_ush_seq).size > 0.9*ush_len) and
                    (nucleotides.lcs(dg_dsh_seq, d_dsh_seq).size > 0.9*dsh_len)):
                    
                    # Then add them to the same group
                    dg.append(datum)
                    to_pop.append(di)
            
            for di in sorted(to_pop, reverse=True):
                datum_list.pop(di)
        
        # Please note:
        #   This "finished" 'datum_groups' could use some work.
        #   The data structure of these linked 'Datum' objects should allow for branching
        #   For instance:
        #     wt1(1,0)──ko1(1,1)──ko1(2,1)──ki1(2,2)
        #      └────────ko1(1,1)──ko1(2,1)──ki1(2,2)
        #   I forgot why, but it is important
        #    ------A-----A-----  contig has 2 locations that should be targeted
        #                        Should PCR be allele/paralog-specific?
        #                        Or should the PCR try to encompass both alleles/paralogs?
        logging.info('datum_groups (finished):')
        for dg in datum_groups:
            logging.info('  ' + str(dg))
        
        logging.info('datum_list (finished):')
        if (len(datum_list) == 0):
            logging.info('  EMPTY')
        else:
            for datum in datum_list:
                logging.info('  ' + str(datum))
        
        ######## End code for linking the loci ########
        
        # 'datum_groups' contains the linked loci
        # Let's find the longest common substring in the far_upstream and far_downstream regions of each
        # We also need the distance (in nt) between this LCS and the feature/insert
        # (for both the upstream and downstream)
        # (so we can do proper amplicon size calculations later)
        pcr_regions = []
        
        # pcr_region_positions[datum_group index][datum index] = [fus_start, fus_end, fds_start, fds_end]
        pcr_region_positions = []
        
        for dg in datum_groups:
            # Find the upstream LCS
            us_dist = 600
            us_done = False
            fus0 = None
            while (us_done == False):
                fus0 = genome_contigs_list[dg[0].genome_r][dg[0].genome_contig][max(0, dg[0].ush_start-us_dist):dg[0].ush_start]
                for di, datum in enumerate(dg[1:]):
                    if ((us_dist > dg[0].ush_start) or (us_dist > datum.ush_start)):
                        logging.info('too far')
                        us_done = True
                    fus = genome_contigs_list[datum.genome_r][datum.genome_contig][max(0, datum.ush_start-us_dist):datum.ush_start]
                    m = nucleotides.lcs(fus0, fus)
                    logging.info('us_dist = ' + str(us_dist))
                    logging.info('   fus0 = ' + fus0)
                    logging.info(('fus'+str(di+1)).rjust(7) +' = ' + fus)
                    logging.info('      m = ' + str(m))
                    if (m.size < 200):
                        us_dist += 100
                        break
                    else:
                        fus0 = fus0[m.a:m.a+m.size]
                else:
                    us_done = True
            
            logging.info('far_upstream_dist: ' + str(us_dist))
            logging.info('far_upstream_seq >= 200: ' + fus0)
            
            # Find the downstream LCS
            ds_dist = 600
            ds_done = False
            fds0 = None
            while (ds_done == False):
                fds0 = genome_contigs_list[dg[0].genome_r][dg[0].genome_contig][dg[0].dsh_end:dg[0].dsh_end+ds_dist]
                for di, datum in enumerate(dg[1:]):
                    if ((dg[0].dsh_end+ds_dist > len(genome_contigs_list[dg[0].genome_r][dg[0].genome_contig])) or (datum.dsh_end+ds_dist > len(genome_contigs_list[datum.genome_r][datum.genome_contig]))):
                        logging.info('too far')
                        ds_done = True
                    fds = genome_contigs_list[datum.genome_r][datum.genome_contig][datum.dsh_end:datum.dsh_end+ds_dist]
                    m = nucleotides.lcs(fds0, fds)
                    logging.info('us_dist = ' + str(ds_dist))
                    logging.info('   fus0 = ' + fds0)
                    logging.info(('fus'+str(di+1)).rjust(7) +' = ' + fds)
                    logging.info('      m = ' + str(m))
                    if (m.size < 200):
                        ds_dist += 100
                        break
                    else:
                        fds0 = fds0[m.a:m.a+m.size]
                else:
                    ds_done = True
            
            logging.info('far_downstream_dist: ' + str(ds_dist))
            logging.info('far_downstream_seq >= 200: ' + fds0)
            
            # Add PCR regions to data structure
            pcr_regions.append((fus0, fds0))
            
            # Need to calculate the position of 'sF' and 'sR' within each genome associated with the Datum
            pcr_region_positions.append([])
            for datum in dg:
                pcr_region_positions[-1].append([])
                
                # Add sF_start, sF_end
                temp_var = max(0, datum.ush_start-us_dist)
                fus_start = genome_contigs_list[datum.genome_r][datum.genome_contig][temp_var:datum.ush_start].rindex(fus0) + temp_var
                fus_end = fus_start + len(fus0)
                pcr_region_positions[-1][-1].append(fus_start)
                pcr_region_positions[-1][-1].append(fus_end)
                
                # Add sR_start, sR_end
                fds_start = genome_contigs_list[datum.genome_r][datum.genome_contig][datum.dsh_end:datum.dsh_end+ds_dist].index(fds0) + datum.dsh_end
                fds_end = fds_start + len(fds0)
                pcr_region_positions[-1][-1].append(fds_start)
                pcr_region_positions[-1][-1].append(fds_end)
            
            logging.info("checking 'sF' region:")
            logging.info("              fus0: " + fus0)
            for prpi2, prp2 in enumerate(pcr_region_positions[-1]):
                logging.info(str(prpi2) + ' ' + str(prp2[:2]).rjust(16) + ": " + genome_contigs_list[dg[prpi2].genome_r][dg[prpi2].genome_contig][prp2[0]:prp2[1]] + ' ' + str(dg[prpi2]))
            logging.info("checking 'sR' region:")
            logging.info("              fds0: " + fds0)
            for prpi2, prp2 in enumerate(pcr_region_positions[-1]):
                logging.info(str(prpi2) + ' ' + str(prp2[2:]).rjust(16) + ": " + genome_contigs_list[dg[prpi2].genome_r][dg[prpi2].genome_contig][prp2[2]:prp2[3]] + ' ' + str(dg[prpi2]))
        
        # Identify the feature/insert sequence where the 'rN-oF', 'rN-oR', 'rN-iF', 'rN-iR' primers should be located
        # and store them in 'insert_seqs'
        Insert = namedtuple('Insert', ['genome_r', 'genome_contig', 'seq', 'us_seq', 'ds_seq', 'fus_dist', 'fds_dist', 'type'])
        max_primer_length = 35
        for i, dg in enumerate(datum_groups):
            fus_seq, fds_seq = pcr_regions[i]
            
            # insert_seqs = [
            #     Insert(qname='ko-dDNA', genome_r=0, genome_contig='chr1',               seq='ACGTAACA') 
            #     Insert(qname='ki-dDNA', genome_r=1, genome_contig='chr1-r1[ko]',        seq='ACGTAACA')
            #     Insert(qname='ki-dDNA', genome_r=2, genome_contig='chr1-r1[ko]-r2[ki]', seq='CGATAAGC')
            # ]
            insert_seqs = []
            
            # Go through every datum, and select all the 'before' ones
            # (because they have the full homology regions specified)
            # And use those to calculate the 'before'=feature, and 'after'=insert
            # sequences (with their up/downstream flanking regions)
            for di, datum in enumerate(dg):
                if (datum.ins_start != None):
                    # Add feature
                    insert_seqs.append(Insert(
                        datum.genome_r, # genome_r
                        datum.genome_contig, # genome_contig
                        genome_contigs_list[datum.genome_r][datum.genome_contig][datum.ush_end:datum.dsh_start], # seq
                        genome_contigs_list[datum.genome_r][datum.genome_contig][max(0, datum.ush_end-(max_primer_length-1)):datum.ush_end], # us_seq
                        genome_contigs_list[datum.genome_r][datum.genome_contig][datum.dsh_start:datum.dsh_start+(max_primer_length-1)], # ds_seq
                        datum.ush_end - pcr_region_positions[i][di][1], # fus_dist
                        pcr_region_positions[i][di][2] - datum.dsh_start, # fds_dist
                        'b' # type
                    ))
                    
                    # Add insert
                    insert_seqs.append(Insert(
                        datum.genome_r, # genome_r
                        datum.genome_contig, # genome_contig
                        #datum.dDNA_r, # genome_r <----------------------- may need to change this to be the genome (not the dDNA)
                        #datum.dDNA_contig, # genome_contig <------------- may need to change this to be the genome (not the dDNA)
                        dDNA_contigs_list[datum.dDNA_r][datum.dDNA_contig][datum.ins_start:datum.ins_end], # seq
                        dDNA_contigs_list[datum.dDNA_r][datum.dDNA_contig][max(0, datum.ins_start-(max_primer_length-1)):datum.ins_start], # us_seq
                        dDNA_contigs_list[datum.dDNA_r][datum.dDNA_contig][datum.ins_end:datum.ins_end+(max_primer_length-1)], # ds_seq
                        datum.ush_end - pcr_region_positions[i][di][1], # fus_dist
                        pcr_region_positions[i][di][2] - datum.dsh_start, # fds_dist
                        'a' # type
                    ))
                # Essentially, this current loop does the following:
                #   if (di == 0): (a 'before')
                #     Then get the original feature seq (for genome_r=0)
                #     And get the insert seq (for genome_r=1)
                #   if (di == 1), (an 'after')
                #     then skip
                #   if (di == 2): (a 'before')
                #     Then get the feature seq (for genome_r=1)
                #     And get the insert seq (for genome_r=2)
                #   if (di == 3): (an 'after)
                #     then skip
            
            # Now that we've identified the sequences where the shared PCR primers should reside ('sF' & 'sR'),
            # As well as the sequences of the feature/insert and their up/downstream flanking sequences
            # (which are used to find primers that span the junctions),
            # we do the actual cPCR calculations
            
            # Old code for 6-set:
            #initial_pair_list, final_pair_list = self.make_primer_set(args, qname, us_seq, ds_seq, insert_seq, feature_seq, q_hih_seq, s_ush_seq, q_ush_seq, s_dsh_seq, q_dsh_seq)
            
            #                                             shared_forward  shared_reverse  features/inserts
            pair_list, insert_pair_list = self.calculate_them_primers(args, fus_seq,        fds_seq,        insert_seqs)
            
            # Filter 'pair_list' to get the top 10
            pair_list = sorted(pair_list, reverse=True)[:10]
            
            # Print the primers for table output 1
            for ppli, (w, pp_list) in enumerate(pair_list):
                for pp_i, pp in enumerate(pp_list):
                    amp_name = '-'
                    if pp: # primer_set_index, primer_pair_index, locus
                        if (pp.forward_primer.name == 'sF'):
                            if (pp.reverse_primer.name == 'sR'):
                                amp_name = 'A'
                            else:
                                amp_name = 'B'
                        else:
                            if (pp.reverse_primer.name == 'sR'):
                                amp_name = 'C'
                            else:
                                amp_name = 'D'
                        output1.append([ppli, pp_i, i, amp_name, 'Round', pp.forward_primer.name, pp.reverse_primer.name, pp.get_amplicon_size(), pp.get_tms(), 'Template'])
                    else:
                        output1.append([ppli, pp_i, i, amp_name, '-', '-', '-', '-', '-', '-'])
            
            # Filter 'insert_pair_list' to get the top 10 (They are already sorted)
            insert_pair_list = [x[:10] for x in insert_pair_list]
            for ip_r, iF_iR_paired_primers in enumerate(insert_pair_list):
                for ppli, pp in enumerate(iF_iR_paired_primers):
                    amp_name = 'D'
                    if pp:
                        output1.append([ppli, ip_r, i, amp_name, 'Round', pp.forward_primer.name, pp.reverse_primer.name, pp.get_amplicon_size(), pp.get_tms(), 'Template'])
                    else:
                        output1.append([ppli, ip_r, i, amp_name, '-', '-', '-', '-', '-', '-'])
            
            # Print the primers for table output 2
            for ppli, (w, pp_list) in enumerate(pair_list):
                for pp_i, pp in enumerate(pp_list):
                    if pp:
                        # 'f_locations' and 'r_locations' should have the same length
                        f_locations = self.get_primer_location(genome_fasta_file_list, genome_contigs_list, pp.forward_primer.sequence)
                        r_locations = self.get_primer_location(genome_fasta_file_list, genome_contigs_list, pp.reverse_primer.sequence)
                        for loc_i in range(len(f_locations)):
                            for loc in f_locations[loc_i]:
                                output2.append([ppli, pp_i, i, pp.forward_primer.name, pp.forward_primer.sequence, loc[0], loc[1], loc[2], loc[3], loc[4]])
                            for loc in r_locations[loc_i]:
                                output2.append([ppli, pp_i, i, pp.reverse_primer.name, pp.reverse_primer.sequence, loc[0], loc[1], loc[2], loc[3], loc[4]])
                        #output2.append([ppli, pp_i, i, pp.reverse_primer.name, pp.reverse_primer.sequence, 'File', 'Contig', 'Start', 'End', '-'])
                    else:
                        output2.append([ppli, pp_i, i, '-', '-', '-', '-', '-', '-', '+'])
                        output2.append([ppli, pp_i, i, '-', '-', '-', '-', '-', '-', '-'])
            
            for ip_r, iF_iR_paired_primers in enumerate(insert_pair_list):
                for ppli, pp in enumerate(iF_iR_paired_primers):
                    if pp:
                        # 'f_locations' and 'r_locations' should have the same length
                        f_locations = self.get_primer_location(genome_fasta_file_list, genome_contigs_list, pp.forward_primer.sequence)
                        r_locations = self.get_primer_location(genome_fasta_file_list, genome_contigs_list, pp.reverse_primer.sequence)
                        for loc_i in range(len(f_locations)):
                            for loc in f_locations[loc_i]:
                                output2.append([ppli, ip_r, i, pp.forward_primer.name, pp.forward_primer.sequence, loc[0], loc[1], loc[2], loc[3], loc[4]])
                            for loc in r_locations[loc_i]:
                                output2.append([ppli, ip_r, i, pp.reverse_primer.name, pp.reverse_primer.sequence, loc[0], loc[1], loc[2], loc[3], loc[4]])
                    else:
                        output2.append([ppli, ip_r, i, '-', '-', '-', '-', '-', '-', '+'])
                        output2.append([ppli, ip_r, i, '-', '-', '-', '-', '-', '-', '-'])
        
        # Output header information for first output table
        print('#                                          Genome')
        print('# Amplicon  ──upstream─┐┌─homology─┐┌──insert/feature──┐┌─homology─┐┌─downstream──')
        print('#        A   sF ===>····················································<=== sR')
        print('#        B   sF ===>···················<=== rN-oR')
        print('#        C                             rN-oF ===>·······················<=== sR')
        print('#        D                        rN-iF ===>······<=== rN-iR')
        print('#')
        print('# F=forward, R=reverse')
        print('# s=shared, o=outer, i=inner')
        print('# rN=round number')
        print('#')
        print('# '+ '\t'.join(['Set', 'Index', 'Locus', 'Amplicon', 'Round', 'F', 'R', 'Size', 'Tm', 'Template']))
        for line in output1:
            print('\t'.join(map(str, line)))
        
        # Output header information for second output table
        print('# ' + '\t'.join(['Set', 'Index', 'Locus', 'Primer', 'Sequence', 'File', 'Contig', 'Start', 'End', 'Strand']))
        for line in output2:
            print('\t'.join(map(str, line)))
        
                        #####################
                        # Cut code was here #
                        #####################
                        
                        # This code output GenBank '*.gb' files
            
    def calculate_them_primers(self, args, far_us_seq, far_ds_seq, insert_seqs):
        """
        """
        tm_max_difference = 4.0
        amplicon_size = (200, 900)
        tm_range = (52, 64)
        primer_length_range = (19,32)
        min_delta_g = -5.0
        
        subset_size = 1000 #35 #200 # Temporary size limit for the number of primers that go into the pair() function
        
        logging.info('Scanning the regions (shared upstream (sF), shared downstream (sR), feature/insert (rN-oF,rN-oR,rN-iF,rN-iR) for all decent primers')
        
        # Make the 'temp_folder' path
        temp_folder = os.path.join('/dev/shm/addtag', os.path.basename(args.folder))
        
        logging.info("Scanning far upstream for 'sF' primers:")
        sF_list = sorted(
            args.selected_oligo.scan(far_us_seq, 'left', primer_size=primer_length_range, tm_range=tm_range, min_delta_g=min_delta_g, folder=temp_folder, time_limit=args.primer_scan_limit),
            key=lambda x: x.weight,
            reverse=True
        )
        logging.info('  len(sF_list) = {}'.format(len(sF_list)))
        logging.info('  sF: skipping {}/{} calculated primers'.format(max(0, len(sF_list)-subset_size), len(sF_list)))
        sF_list = sF_list[:subset_size]
        # Need to add a condition that if (len(sF_list) < N), then it should re-evaluate the far-upstream region to try to get
        # more sequence to search.
        
        logging.info("Scanning far downstream for 'sR' primers:")
        sR_list = sorted(
            args.selected_oligo.scan(far_ds_seq, 'right', primer_size=primer_length_range, tm_range=tm_range, min_delta_g=min_delta_g, folder=temp_folder, time_limit=args.primer_scan_limit),
            key=lambda x: x.weight,
            reverse=True
        )
        logging.info('  len(sR_list) = {}'.format(len(sR_list)))
        logging.info('  sR: skipping {}/{} calculated primers'.format(max(0, len(sR_list)-subset_size), len(sR_list)))
        sR_list = sR_list[:subset_size]
        # Need to add condition that if there aren't enough putative 'sR' primers, then the sR region is expanded
        # Otherwise, the script can just end prematurely
        
        insert_list = []
        for ins in insert_seqs:
            logging.info("Scanning feature/insert for 'rN-oF', 'rN-oR', 'rN-iF', 'rN-iR' primers:")
            iF_list = sorted(
                args.selected_oligo.scan(ins.seq, 'left',  primer_size=primer_length_range, tm_range=tm_range, min_delta_g=min_delta_g, folder=temp_folder, us_seq=ins.us_seq, ds_seq=ins.ds_seq, time_limit=args.primer_scan_limit),
                key=lambda x: x.weight,
                reverse=True
            )
            logging.info('  len(iF_list) = {}'.format(len(iF_list)))
            logging.info('  insert_F: skipping {}/{} calculated primers'.format(max(0, len(iF_list)-subset_size), len(iF_list)))
            iF_list = iF_list[:subset_size]
            
            iR_list = sorted(
                args.selected_oligo.scan(ins.seq, 'right', primer_size=primer_length_range, tm_range=tm_range, min_delta_g=min_delta_g, folder=temp_folder, us_seq=ins.us_seq, ds_seq=ins.ds_seq, time_limit=args.primer_scan_limit),
                key=lambda x: x.weight,
                reverse=True
            )
            logging.info('  len(iR_list) = {}'.format(len(iR_list)))
            logging.info('  insert_R: skipping {}/{} calculated primers'.format(max(0, len(iR_list)-subset_size), len(iR_list)))
            iR_list = iR_list[:subset_size]
            
            insert_list.append([iF_list, iR_list])
        
        # Select the best subset of primers so the pairing doesn't take too long
        # This will cap the maximum number of primer pairs for each region span
        # to be 100 x 100 = 10,000
        
        # Do the pair calculations that involve 'sF' and 'sR'
        logging.info("Calculating: 'sF' 'sR' paired primers...")
        sF_sR_paired_primers = args.selected_oligo.pair(sF_list, sR_list, amplicon_size=(0, amplicon_size[1]), tm_max_difference=tm_max_difference, intervening=0, min_delta_g=min_delta_g, folder=temp_folder, time_limit=args.primer_pair_limit)
        
        for pp in sF_sR_paired_primers:
            # Re-calculate the PrimerPair weights to prefer the smallest amplicon sizes
            pp.weight = pp.get_weight(minimize=True)
            
            # Re-name the primers so when they are printed, they are easy to distinguish
            pp.forward_primer.name = 'sF'
            pp.reverse_primer.name = 'sR'
        
        # Sort by weight
        sF_sR_paired_primers = sorted(
            sF_sR_paired_primers,
            key=lambda x: x.get_joint_weight(),
            reverse=True
        )
        logging.info('  len(sF_sR_paired_primers) = {}'.format(len(sF_sR_paired_primers)))
        
        # pair_list[ins index] = [sF_oR_paired_primers, oF_sR_paired_primers, iF_iR_paired_primers]
        pair_list = []
        insert_pair_list = []
        
        #for i, (iF_list, iR_list) in enumerate(insert_list):
        for i in range(len(insert_seqs)):
            iF_list, iR_list = insert_list[i]
            ins = insert_seqs[i]
            
            logging.info('Pairing primers: ' + str(i))
        
            # If there is a hard constraint for primers that should be used
            # That is, if flanktag primers should be used
            if (len(args.primers) > 0):
                pass
            # Otherwise, no flanktags are specified
            else:
                # sF rN-oR
                logging.info("  Calculating: 'sF' 'r"+str(ins.genome_r)+ins.type+"-oR' paired_primers...")
                sF_oR_paired_primers = sorted(
                    args.selected_oligo.pair(sF_list, iR_list, amplicon_size=amplicon_size, tm_max_difference=tm_max_difference, intervening=ins.fus_dist, min_delta_g=min_delta_g, folder=temp_folder, time_limit=args.primer_pair_limit),
                    key=lambda x: x.get_joint_weight(),
                    reverse=True
                )
                for pp in sF_oR_paired_primers:
                    pp.forward_primer.name, pp.reverse_primer.name = 'sF', 'r'+str(ins.genome_r)+ins.type+'-oR'
                logging.info('  len(sF_oR_paired_primers) = {}'.format(len(sF_oR_paired_primers)))
            
                # rN-0F sR
                logging.info("  Calculating: 'r"+str(ins.genome_r)+ins.type+"-oF' 'sR' paired_primers...")
                oF_sR_paired_primers = sorted(
                    args.selected_oligo.pair(iF_list, sR_list, amplicon_size=amplicon_size, tm_max_difference=tm_max_difference, intervening=ins.fds_dist, min_delta_g=min_delta_g, folder=temp_folder, time_limit=args.primer_pair_limit),
                    key=lambda x: x.get_joint_weight(),
                    reverse=True
                )
                for pp in oF_sR_paired_primers:
                    pp.forward_primer.name, pp.reverse_primer.name = 'r'+str(ins.genome_r)+ins.type+'-oF', 'sR'
                logging.info('  len(oF_sR_paired_primers) = {}'.format(len(oF_sR_paired_primers)))
                
                # rN-iF rN-iR
                logging.info("  Calculating: 'r"+str(ins.genome_r)+ins.type+"-iF' 'r"+str(ins.genome_r)+ins.type+"-iR' paired_primers...")
                iF_iR_paired_primers = sorted(
                    args.selected_oligo.pair(iF_list, iR_list, amplicon_size=amplicon_size, tm_max_difference=tm_max_difference, intervening=0, same_template=True, min_delta_g=min_delta_g, folder=temp_folder, time_limit=args.primer_pair_limit),
                    key=lambda x: x.get_joint_weight(),
                    reverse=True
                )
                for pp in iF_iR_paired_primers:
                    pp.forward_primer.name, pp.reverse_primer.name = 'r'+str(ins.genome_r)+ins.type+'-iF', 'r'+str(ins.genome_r)+ins.type+'-iR'
                logging.info('  len(iF_iR_paired_primers) = {}'.format(len(iF_iR_paired_primers)))
                
                # Add calculated primer pairs to the list
                pair_list.append([sF_oR_paired_primers, oF_sR_paired_primers, iF_iR_paired_primers, ins])
                insert_pair_list.append(iF_iR_paired_primers)
        
        # If a primer in feature/insert of one round ALSO is identical
        # with a primer in another round,
        # Then it should only be weighted for a SINGLE entry
        starting_set, finished_set = self.calculate_them_best_set(args, sF_sR_paired_primers, pair_list)
        
        return finished_set, insert_pair_list
    
    def get_primer_location(self, filenames, genomes, primer_sequence):
        #[filename, contig, start, end, strand]
        locations = []
        for r in range(len(filenames)):
            locations.append([])
            for contig, seq in genomes[r].items():
                for m in regex.finditer(primer_sequence, seq, flags=regex.IGNORECASE):
                    locations[-1].append([filenames[r], contig, m.start(), m.end(), '+'])
                for m in regex.finditer(nucleotides.rc(primer_sequence), seq, flags=regex.IGNORECASE):
                    locations[-1].append([filenames[r], contig, m.start(), m.end(), '-'])
        return locations
    
    def calculate_them_best_set(self, args, sF_sR_paired_primers, pair_list, max_iterations=1000):
        """
        Takes input primer sets and calculates their weights
        
        This starts with the best from a sorted list of 'sF' 'sR' primer pairs
        Then it optimizes by swapping the worst component a finite number of times
        (or until the delta drops below a certain amount)
        """
        # Define helper functions specific to only this method
        def nsum(n):
            return n*(n+1)/2
        
        def rank_order(x, reverse=False, shift=0):
            return [y+shift for y in sorted(range(len(x)), key=x.__getitem__, reverse=reverse)]
        
        def rank(x, reverse=False, shift=0):
            return [z[0]+shift for z in sorted(enumerate(sorted(enumerate(x), key=lambda w: w[1], reverse=reverse)), key=lambda y: y[1][0])]
        
        def rank_dist(x):
            s = nsum(len(x))
            return [y/s for y in rank(x, reverse=False, shift=1)]
        
        def product(x):
            z = 1
            for y in x:
                z *= y
            return z
        
        def random_choices(population, weights, k=1):
            """
            Return a k sized list of population elements chosen with replacement.
            """
            weight_sum = sum(weights)
            choices = zip(population, weights)
            values = []
            for i in range(k):
                r = random.uniform(0, weight_sum)
                upto = 0
                for c, w in choices:
                    if upto + w >= r:
                        values.append(c)
                        break
                    upto += w
                else:
                    values.append(random.choice(population))
            return values
        
        def filter_primer_pairs(pairs, forward=None, reverse=None):
            ooo = []
            
            if (forward and reverse):
                for pp in pairs:
                    if ((pp.forward_primer.sequence == forward.sequence) and (pp.reverse_primer.sequence == reverse.sequence)):
                        ooo.append(pp)
            elif forward:
                for pp in pairs:
                    if (pp.forward_primer.sequence == forward.sequence):
                        ooo.append(pp)
            elif reverse:
                for pp in pairs:
                    if (pp.reverse_primer.sequence == reverse.sequence):
                        ooo.append(pp)
            else:
                ooo = pairs
            
            return ooo
        
        def random_primer_by_weight_old(pairs, forward=None, reverse=None):
            ooo = filter_primer_pairs(pairs, forward, reverse)
            
            if (len(ooo) == 0):
                return None
            elif ('choices' in random.__all__):
                return random.choices(ooo, [x.get_joint_weight() for x in ooo])[0]
            else:
                return random_choices(ooo, [x.get_joint_weight() for x in ooo])[0]
        
        def random_primer_by_weight(pairs):
            if (len(pairs) == 0):
                return None
            elif ('choices' in random.__all__):
                return random.choices(pairs, [x.get_joint_weight() for x in pairs])[0]
            else:
                return random_choices(pairs, [x.get_joint_weight() for x in pairs])[0]
        
        # This should be in a loop...
        #sF_oR_paired_primers, oF_sR_paired_primers, iF_iR_paired_primers, ins = pair_list[0]
        
        # Set initial iteration count to zero
        iteration_count = 0
        
        # Create lists to hold 5-set of primer pairs, and their joint weight
        starting_set = [] # initial set
        finished_set = [] # final set
        
        sF_sR_iterator = iter(sF_sR_paired_primers)
        
        if ((max_iterations != None) and (len(sF_sR_paired_primers) > max_iterations)):
            logging.info('sF_sR_paired_primers: skipping {}/{} calculated primer pairs'.format(max(0, len(sF_sR_paired_primers)-max_iterations), len(sF_sR_paired_primers)))
        
        while (iteration_count < max_iterations):
            # Increment the count
            iteration_count += 1
            
            # Create random 6-sets, with each pair's probability determined by its joint weight
            #uf_dr_pair = random_primer_by_weight(sF_sR_paired_primers) # comment this out to prevent random sampling of this one
            try:
                sF_sR_pair = next(sF_sR_iterator)
                #sF_sR_pair.forward_primer.name = 'sF'
                #sF_sR_pair.reverse_primer.name = 'sR'
            except StopIteration:
                break
            
            logging.info('loop {}:'.format(iteration_count))
            
            # Populate set with all primer pairs that have 'sF' and 'sR'
            pp_sources = []
            pp_setN = []
            p_setN = []
            for sF_oR_paired_primers, oF_sR_paired_primers, iF_iR_paired_primers, ins in pair_list:
                # Add the 'oR' from 'sF-oR' pair
                pp_sources.append(filter_primer_pairs(sF_oR_paired_primers, forward=sF_sR_pair.forward_primer))
                #for pp in pp_sources[-1]:
                #    if pp:
                #        pp.forward_primer.name, pp.reverse_primer.name = 'sF', 'r'+str(ins.genome_r)+ins.type+'-oR'
                pp = random_primer_by_weight(pp_sources[-1]) # Pick a random primer by weight
                pp_setN.append(pp) # Will be 'None' if 'pp_sources[-1]' is empty
                if pp:
                    p_setN.append(pp.reverse_primer)
                else:
                    p_setN.append(None)
                
                # Add the 'oF' from 'oF-sR' pair
                pp_sources.append(filter_primer_pairs(oF_sR_paired_primers, reverse=sF_sR_pair.reverse_primer))
                #for pp in pp_sources[-1]:
                #    if pp:
                #        pp.forward_primer.name, pp.reverse_primer.name = 'r'+str(ins.genome_r)+ins.type+'-oF', 'sR'
                pp = random_primer_by_weight(pp_sources[-1]) # Pick a random primer by weight
                pp_setN.append(pp) # Will be 'None' if 'pp_sources[-1]' is empty
                if pp:
                    p_setN.append(pp.forward_primer)
                else:
                    p_setN.append(None)
            
            logging.info('  length of sources: ' + str([len(x) for x in pp_sources]))
            
            # # Pick random primer by weight for each source
            # pp_setN = []
            # p_setN = []
            # for pp_list in pp_sources:
            #     pp = random_primer_by_weight(pp_list)
            #     pp_setN.append(pp) # Will be 'None' if 's' is empty
            # 
            # # Break down the 'pp_setN' which contains primer pairs into a set of primers
            # p_setN = []
            # for pp in pp_setN:
            #     if pp:
            #         p_setN.append(pp.reverse_primer)
            #         p_setN.append(pp.forward_primer)
            #     else:
            #         p_setN.append(None)
            
            group_weight = args.selected_oligo.group_weight([sF_sR_pair.forward_primer, sF_sR_pair.reverse_primer]+p_setN)
            joint_weight = group_weight * product(x.get_joint_weight() for x in [sF_sR_pair]+pp_setN if x)
            starting_set.append((joint_weight, [sF_sR_pair]+pp_setN))
            
            logging.info('  starting_set: ' + str(starting_set[-1]))
            
            # iteratively improve until a local maxima is found
            # By continually swapping out the least-weighted component primer
            min_weight_delta = 1e-20
            weight_delta = 1
            while(weight_delta > min_weight_delta):
                # Get list of indices from smallest weight to largest weight (excluding 'sF' and 'sR')
                wi_order = rank_order([x.weight if (x != None) else math.inf for x in p_setN])
                
                # If primer is 'None', then its weight will be 'math.inf', and the index
                # corresponding to it will be the last element in wi_order:
                #   rank_order([1, 2, 2.5, 0.001, math.inf]) # [3, 0, 1, 2, 4]
                
                # Start at the worst, and go to the next-worst, then next, then next
                for wi in wi_order:
                    # Copy the list of primer pairs
                    pp_setNN = pp_setN[:]
                    
                    # Swap the worst-performing primer with an alternative.
                    # If the alternative gives a better weight, then keep it
                    # and break out of the loop
                    for source_pp in pp_sources[wi]:
                        pp_setNN[wi] = source_pp
                        
                        p_setNN = []
                        for pp in pp_setNN:
                            if pp:
                                if (wi % 2 == 0):
                                    p_setNN.append(pp.reverse_primer)
                                else:
                                    p_setNN.append(pp.forward_primer)
                            else:
                                p_setNN.append(None)
                        
                        new_group_weight = args.selected_oligo.group_weight([sF_sR_pair.forward_primer, sF_sR_pair.reverse_primer]+p_setNN)
                        new_joint_weight = new_group_weight * product(x.get_joint_weight() for x in [sF_sR_pair]+pp_setNN if x)
                        
                        if (new_joint_weight > joint_weight):
                            logging.info('           set: ' + str((new_joint_weight, [sF_sR_pair]+pp_setNN)))
                            
                            weight_delta = new_joint_weight - joint_weight
                            
                            # replace values for next iteration
                            pp_setN = pp_setNN
                            joint_weight = new_joint_weight
                            break
                    else:
                        weight_delta = 0
                    
                    if (weight_delta > 0):
                        break
            
            # Does not re-calculate 'pp.weight', thus the amplicon_size will not be further penalized
            # However, the 'pp.get_amplicon_size()' will accurately reflect the updated 'pp.intervening' length
            #sF_insert_sR_pair = copy.deepcopy(sF_sR_pair)
            #sF_insert_sR_pair.intervening = len(q_hih_seq)
            #sF_feature_sR_pair = copy.deepcopy(sF_sR_pair)
            #sF_feature_sR_pair.intervening = len(s_ush_seq)+len(sseq_feature)+len(s_dsh_seq)
            
            # Duplicate 'sF_sR_pair' for each template, giving it the proper 'intervening' length (amplicon size)
            sF_sR_pair_list = []
            for sF_oR_paired_primers, oF_sR_paired_primers, iF_iR_paired_primers, ins in pair_list:
                sF_sR_pair_list.append(copy.deepcopy(sF_sR_pair))
                sF_sR_pair_list[-1].intervening = ins.fus_dist + len(ins.seq) + ins.fds_dist
            
            #finished_set.append((joint_weight, [sF_sR_pair]+pp_setN)) # <-- need to add calculations for feature vs insert
            finished_set.append((joint_weight, sF_sR_pair_list+pp_setN)) # two copies of 'sF_sR_pair'
        
        return starting_set, finished_set
    
    def _parser_general(self):
        '''general parser'''
        # Create the parent argument parser
        parser = argparse.ArgumentParser(
            description=__description__,
            epilog=__epilog__,
            formatter_class=CustomHelpFormatter
        )
        
        # Create the subparsers
        subparsers = parser.add_subparsers(metavar='action', help='choose an action to perform (required)')
        
        # Change the help text of the "-h" flag
        parser._actions[0].help='Show this help message and exit.'
        
        # Special version action optional argument
        parser.add_argument("-v", "--version", action='version',
            help="Show program's version number and exit.",
            version='{__program__} {__version__} (revision {__revision__})'.format(**globals()))
        
        return parser, subparsers
    
    def _parser_glossary(self, subparsers):
        ''' "glossary" parser '''
        __glossary_description__ = "description:\n  Show glossary of common CRISPR/Cas terms, then exit."
        __glossary_help__ = "Show glossary of common CRISPR/Cas terms, then exit."
        parser_glossary = subparsers.add_parser('glossary',
            description=__glossary_description__,
            #epilog=__glossary_epilog__,
            formatter_class=CustomHelpFormatter,
            help=__glossary_help__
        )
        parser_glossary.set_defaults(func=self._glossary)
        
        # Change the help text of the "-h" flag
        parser_glossary._actions[0].help='Show this help message and exit.'
        
        # Special version action optional argument
        parser_glossary.add_argument("-v", "--version", action='version',
            help="Show program's version number and exit.",
            version='{__program__} {__version__} (revision {__revision__})'.format(**globals()))
        
        # Add mandatory arguments
        # Add optional arguments
        return parser_glossary
    
    def _parser_motifs(self, subparsers):
        ''' "motifs" parser '''
        __motifs_description__ = "description:\n  Show list of common CRISPR/Cas SPACER>PAM arrangements, then exit."
        __motifs_help__ = "Show list of common CRISPR/Cas SPACER>PAM arrangements, then exit."
        parser_motifs = subparsers.add_parser('motifs',
            description=__motifs_description__,
            #epilog=__motif_epilog__,
            formatter_class=CustomHelpFormatter,
            help=__motifs_help__
        )
        parser_motifs.set_defaults(func=self._motifs)
        
        # Change the help text of the "-h" flag
        parser_motifs._actions[0].help='Show this help message and exit.'
        
        # Special version action optional argument
        parser_motifs.add_argument("-v", "--version", action='version',
            help="Show program's version number and exit.",
            version='{__program__} {__version__} (revision {__revision__})'.format(**globals()))
        
        # Add mandatory arguments
        # Add optional arguments
        return parser_motifs
    
    def _parser_algorithms(self, subparsers):
        ''' "algorithms" parser '''
        __algorithms_description__ = "description:\n  Show list of all implemented gRNA evaluation algorithms."
        __algorithms_help__ = "Show list of all implemented gRNA evaluation algorithms."
        parser_algorithms = subparsers.add_parser('algorithms',
            description=__algorithms_description__,
            #epilog=__algorithms_epilog__,
            formatter_class=CustomHelpFormatter,
            help=__algorithms_help__
        )
        parser_algorithms.set_defaults(func=self._algorithms)
        
        # Change the help text of the "-h" flag
        parser_algorithms._actions[0].help='Show this help message and exit.'
        
        # Special version action optional argument
        parser_algorithms.add_argument("-v", "--version", action='version',
            help="Show program's version number and exit.",
            version='{__program__} {__version__} (revision {__revision__})'.format(**globals()))
        
        # Add mandatory arguments
        # Add optional arguments
        
        return parser_algorithms
    
    def _parser_aligners(self, subparsers):
        ''' "aligners" parser '''
        __aligners_description__ = "description:\n  Show list of all supported alignment programs."
        __aligners_help__ = "Show list of all supported alignment programs."
        parser_aligners = subparsers.add_parser('aligners',
            description=__aligners_description__,
            #epilog=__aligners_epilog__,
            formatter_class=CustomHelpFormatter,
            help=__aligners_help__
        )
        parser_aligners.set_defaults(func=self._aligners)
        
        # Change the help text of the "-h" flag
        parser_aligners._actions[0].help='Show this help message and exit.'
        
        # Special version action optional argument
        parser_aligners.add_argument("-v", "--version", action='version',
            help="Show program's version number and exit.",
            version='{__program__} {__version__} (revision {__revision__})'.format(**globals()))
        
        # Add mandatory arguments
        # Add optional arguments
        
        return parser_aligners
    
    def _parser_oligos(self, subparsers):
        ''' "oligos" parser '''
        __oligos_description__ = "description:\n  Show list of all supported oligonucleotide thermodynamics property programs."
        __oligos_help__ = "Show list of all supported oligonucleotide thermodynamics property programs."
        parser_oligos = subparsers.add_parser('oligos',
            description=__oligos_description__,
            #epilog=__oligos_epilog__,
            formatter_class=CustomHelpFormatter,
            help=__oligos_help__
        )
        parser_oligos.set_defaults(func=self._oligos)
        
        # Change the help text of the "-h" flag
        parser_oligos._actions[0].help='Show this help message and exit.'
        
        # Special version action optional argument
        parser_oligos.add_argument("-v", "--version", action='version',
            help="Show program's version number and exit.",
            version='{__program__} {__version__} (revision {__revision__})'.format(**globals()))
        
        # Add mandatory arguments
        # Add optional arguments
        
        return parser_oligos
    
    def _parser_search(self, subparsers):
        ''' "search" parser '''
        # Searches input FASTA for DNA, and outputs a GFF file for use as input for 'generate'
        
        __search_description__ = """\
description:
  Search input FASTA for DNA sequence, and output a GFF file for use as input
  for the 'generate' subroutine.
"""
        __search_help__ = "Search for positions of DNA in FASTA to make GFF file."
        __search_epilog__ = """\
example:
  Running AddTag with the following arguments:
   $ python3 {__program__} search --fasta genome.fasta --query AGTCGCCAAC
     CCAATTAGGAG > search.gff
  
  Will produce a GFF3 output file 'search.gff' like the following:
     chr1A	addtag	search	928299	928308	.	-	.	ID=search_0_0
     chr1B	addtag	search	928338	928347	.	-	.	ID=search_0_1
     chr4A	addtag	search	3652	3661	.	-	.	ID=search_0_2
     chr4B	addtag	search	3652	3661	.	-	.	ID=search_0_3
     chr3A	addtag	search	819381	819391	.	+	.	ID=search_1_0
     chr3B	addtag	search	819367	819377	.	+	.	ID=search_1_1

  You can use the '--identifier' option, as follows:
   $ python3 {__program__} search --fasta genome.fasta --query AGTCGCCAAC
     CCAATTAGGAG --identifier feature1 feature2 > search.gff
  
  This assigns the identifier to each query:
     chr1A	addtag	search	928299	928308	.	-	.	ID=feature1_0
     chr1B	addtag	search	928338	928347	.	-	.	ID=feature1_1
     chr4A	addtag	search	3652	3661	.	-	.	ID=feature1_2
     chr4B	addtag	search	3652	3661	.	-	.	ID=feature1_3
     chr3A	addtag	search	819381	819391	.	+	.	ID=feature2_0
     chr3B	addtag	search	819367	819377	.	+	.	ID=feature2_1
""".format(**globals())
        parser_search = subparsers.add_parser('search',
            description=__search_description__,
            epilog=__search_epilog__,
            formatter_class=CustomHelpFormatter,
            help=__search_help__
        )
        parser_search.set_defaults(func=self._search)
        
        # Change the help text of the "-h" flag
        parser_search._actions[0].help='Show this help message and exit.'
        
        # Special version action optional argument
        parser_search.add_argument("-v", "--version", action='version',
            help="Show program's version number and exit.",
            version='{__program__} {__version__} (revision {__revision__})'.format(**globals()))
        
        # Add mandatory arguments
        required_group = parser_search.add_argument_group('required arguments')
        required_group.add_argument("--fasta", required=True, nargs="+", metavar="*.fasta", type=str,
            help="FASTA files with contigs to find unique gRNA sites. Input FASTA \
            must not be compressed. User can decide whether ambiguous bases can \
            be chosen for gRNA sites. All FASTA sequences should have unique \
            primary headers (everything between the '>' symbol and the first \
            whitespace should be unique). You should include FASTA of the genome \
            and any plasmids.")
        
        required_group.add_argument("--query", required=True, nargs="+", metavar="SEQUENCE",
            type=str, help="One or more sequences to search.")
        
        # Add optional arguments
        parser_search.add_argument("--identifier", metavar="FEATURE", type=str, nargs="+",
            help="Identifier for input query sequences.")
        
        parser_search.add_argument("--tag", metavar="TAG", type=str, default="ID",
            help="GFF3 attribute tag.")
        
        return parser_search
    
    def _parser_feature(self, subparsers):
        ''' "feature" parser '''
        __feature_description__ = "description:\n  Search GFF features for specific text."
        __feature_help__ = "Search GFF features for specific text."
        __feature_epilog__ = """\
example:
  In general, you can search all feature attributes for text as follows:
   $ python3 {__program__} feature --gff genome.gff --query HSP90 > features.gff
""".format(**globals())
#  If you want to limit your search to specific tags, then you could use
#  these parameters:
#   $ python3 {__program__} feature --gff genome.gff --query Gene=GAL4 > features.gff
#"""

        parser_feature = subparsers.add_parser('feature',
            description=__feature_description__,
            epilog=__feature_epilog__,
            formatter_class=CustomHelpFormatter,
            help=__feature_help__
        )
        parser_feature.set_defaults(func=self._feature)
        
        # Change the help text of the "-h" flag
        parser_feature._actions[0].help='Show this help message and exit.'
        
        # Special version action optional argument
        parser_feature.add_argument("-v", "--version", action='version',
            help="Show program's version number and exit.",
            version='{__program__} {__version__} (revision {__revision__})'.format(**globals()))
        
        # Add mandatory arguments
        required_group = parser_feature.add_argument_group('required arguments')
        required_group.add_argument("--gff", required=True, metavar="*.gff", type=str,
            help="GFF file specifying chromosomal features that will be searched.")
        
        required_group.add_argument("--query", required=True, nargs="+", metavar="TEXT",
            type=str, help="One or more words to search for within the GFF file.")
        
        # Add optional arguments
        parser_feature.add_argument("--linked_tags", metavar="TAG", nargs="*",
            type=str, default=['ID', 'Name', 'Alias', 'Parent', 'Gene'],
            help="If a feature is found that matches the query, include other \
            features that have similar values for these tags.")
        
        parser_feature.add_argument("--allow_errors", action="store_true",
            default=False, help="Include matches with minor differences from the query.")
        
        parser_feature.add_argument("--header", action="store_true",
            default=False, help="Begin output with a commented line containing field names.")
        
        return parser_feature
    
    def _parser_extract(self, subparsers):
        ''' "extract" parser '''
        
        __extract_description__ = """\
description:
  Extracts selected sequences from input FASTA by matching their primary
  sequence header (Everything between the '>' and the first whitespace), and
  outputs in FASTA format.
"""
        __extract_help__ = "Search FASTA headers for specific text."
        __extract_epilog__ = """\
example:
  Running AddTag with the following arguments:
   $ python3 {__program__} extract --fasta excision-spacers.fasta
   --query exTarget-33 exTarget-21 > extract.fasta
""".format(**globals())
        parser_extract = subparsers.add_parser('extract',
            description=__extract_description__,
            epilog=__extract_epilog__,
            formatter_class=CustomHelpFormatter,
            help=__extract_help__
        )
        parser_extract.set_defaults(func=self._extract)
        
        # Change the help text of the "-h" flag
        parser_extract._actions[0].help='Show this help message and exit.'
        
        # Special version action optional argument
        parser_extract.add_argument("-v", "--version", action='version',
            help="Show program's version number and exit.",
            version='{__program__} {__version__} (revision {__revision__})'.format(**globals()))
        
        # Add mandatory arguments
        required_group = parser_extract.add_argument_group('required arguments')
        required_group.add_argument("--fasta", required=True, nargs="+", metavar="*.fasta", type=str,
            help="FASTA files with contigs. All FASTA sequences should have unique \
            primary headers (everything between the '>' symbol and the first \
            whitespace should be unique).")
        
        required_group.add_argument("--query", required=True, nargs="+", metavar="TEXT",
            type=str, help="One or more words to search for within sequence headers.")
        
        # Add optional arguments
        parser_extract.add_argument("--allow_errors", action="store_true",
            default=False, help="Include matches with minor differences from the query.")
        
        return parser_extract
    
    def _parser_evaluate(self, subparsers):
        ''' "evaluate" parser '''
        __evaluate_description__ = 'description:\n  Evaluate pre-designed CRISPR/Cas oligonucleotide sequences.'
        __evaluate_help__ = "Evaluate pre-designed CRISPR/Cas oligonucleotide sequences."
        parser_evaluate = subparsers.add_parser('evaluate',
            description=__evaluate_description__,
            #epilog=__evaluate_epilog__,
            formatter_class=CustomHelpFormatter,
            help=__evaluate_help__
        )
        parser_evaluate.set_defaults(func=self._evaluate)
        
        # Change the help text of the "-h" flag
        parser_evaluate._actions[0].help='Show this help message and exit.'
        
        # Special version action optional argument
        parser_evaluate.add_argument("-v", "--version", action='version',
            help="Show program's version number and exit.",
            version='{__program__} {__version__} (revision {__revision__})'.format(**globals()))
        
        # Add required arguments
        required_group = parser_evaluate.add_argument_group('required arguments')
        required_group.add_argument("--fasta", required=True, nargs="+", metavar="*.fasta", type=str,
            help="FASTA files with contigs to find unique gRNA sites. Input FASTA \
            must not be compressed. User can decide whether ambiguous bases can \
            be chosen for gRNA sites. All FASTA sequences should have unique \
            primary headers (everything between the '>' symbol and the first \
            whitespace should be unique). You should include FASTA of the genome \
            and any plasmids.")
        required_group.add_argument("--gff", required=True, metavar="*.gff", type=str,
            help="GFF file specifying chromosomal features that should be \
                 multiplexed together.")
        required_group.add_argument("--folder", required=True, metavar="FOLDER",
            type=str, help="Path of folder to store generated files.")
        
        # Add required, mutually-exclusive group
        me_group = parser_evaluate.add_mutually_exclusive_group(required=True)
        me_group.add_argument("--spacers", metavar="*.fasta",
            help="Evaluate where spacers would target, and for each site, what \
                 its scores would be.")
        
        me_group.add_argument("--gRNAs", metavar="*.fasta",
            help="Evaluate where gRNAs would target, and for each site, what \
                 its scores would be.")
        
        me_group.add_argument("--dDNAs", nargs="+", metavar="*.fasta", type=str,
            help="Homology arms, presence of unique gRNA target site, and which \
                 features these dDNAs would disrupt, and whether-or-not the \
                 dDNA would introduce a mutation.")
        
        me_group.add_argument("--primers", metavar="*.fasta",
            help="Evaluate what amplicon sizes to expect for wt/ko/ki for \
                 input primers. Also evaluate the Tm/sensitivity/dimerization/etc \
                 for the input primers. Does not evaluate these in multiplex.")
        
        # Maybe: make --gRNAs, --dDNAs, and --primers all accept the same number of FASTA files
        #   Then, when it does calculations, it does all of them in serial steps.
        #  step 1: --gRNAs[0], --dDNAs[0], --primers[0]
        #  step 2: --gRNAs[1], --dDNAs[1], --primers[1]
        #  etc...
        
        # Add optional arguments
        
        # If evaluating knock-in gRNA, then include the ko-dDNA with the "--fasta" argument??
        #  addtag evaluate
        #   --ko-gRNA            
        #   --ko-dDNA            Evaluate knock-out dDNA for which DNA features it would knock-out, and evaluate whether it has a unique gRNA target
        #   --ki-gRNA            Evaulate if knock-in gRNA would target the input ko-dDNA
        #   --ki-dDNA            Evaluate if knock-in dDNA will restore wild type correctly, or where the mutations will be
        #   --primers FASTA      Evaluate what amplicon sizes to expect for wt/ko/ki for input primers.
        #                        Also evaluate the Tms/sensitivity/dimerization/etc for the input primers
        #                        (DOES NOT EVALUATE THESE IN MULTIPLEX)
        #                        >feature1 F
        #                        NNNNNNNNNNNNNNNNNNNN
        #                        >feature1 R
        #                        NNNNNNNNNNNNNNNNNNNN
        
        parser_evaluate.add_argument("--motifs", metavar="MOTIF", nargs="+", type=str,
            default=["N{20}>NGG"],
            help="Find only targets with these 'SPACER>PAM' motifs, written from \
            5' to 3'. '>' points toward PAM. IUPAC ambiguities accepted. '{a,b}' \
            are quantifiers. '/' is a sense strand cut, '\\' is an antisense strand \
            cut, and '|' is a double-strand cut. '.' is a base used for positional \
            information, but not enzymatic recognition. Be sure to enclose each \
            motif in quotes so your shell does not interpret STDIN/STDOUT redirection.")
        parser_evaluate.add_argument("--off_target_motifs", metavar="MOTIF", nargs="+", type=str,
            default=[],
            help="Defaults to the same as the on-target motif. Definition syntax is identical.")
        
        return parser_evaluate
    
    def _parser_generate(self, subparsers):
        ''' "generate" parser '''
        __generate_description__ = "description:\n  Design full sets of oligonucleotide sequences for CRISPR/Cas genome engineering experiment."
        __generate_help__ = "Design full sets of oligonucleotide sequences for CRISPR/Cas genome engineering experiment."
        parser_generate = subparsers.add_parser('generate',
            description=__generate_description__,
            #epilog=__generate_epilog__,
            formatter_class=CustomHelpFormatter,
            help=__generate_help__
        )
        parser_generate.set_defaults(func=self._generate)
        
        # Change the help text of the "-h" flag
        parser_generate._actions[0].help='Show this help message and exit.'
        
        # Special version action optional argument
        parser_generate.add_argument("-v", "--version", action='version',
            help="Show program's version number and exit.",
            version='{__program__} {__version__} (revision {__revision__})'.format(**globals()))
        
        # Add required arguments
        required_group = parser_generate.add_argument_group('required arguments')
        required_group.add_argument("--fasta", required=True, nargs="+", metavar="*.fasta", type=str,
            help="FASTA files with contigs to find unique gRNA sites. Input FASTA \
            must not be compressed. User can decide whether ambiguous bases can \
            be chosen for gRNA sites. All FASTA sequences should have unique \
            primary headers (everything between the '>' symbol and the first \
            whitespace should be unique). You should include FASTA of the genome \
            and any plasmids.")
        required_group.add_argument("--gff", required=True, nargs="+", metavar="*.gff", type=str,
            help="GFF files specifying chromosomal features that should be \
            targeted (in multiplex) or excluded.")
        required_group.add_argument("--folder", required=True, metavar="FOLDER",
            type=str, help="Path of folder to store generated files.")
        
        # Add optional arguments
        parser_generate.add_argument("--exhaustive", action="store_true", 
            help="Perform brute force search for optimal gRNA design. \
            This will significantly increase runtime.") # Default to on when not trimming???
        #parser.add_argument("--feature_homolog_regex", metavar="REGEX", type=str, default=None, help="regular expression with capturing group containing invariant feature. Example: '(.*)_[AB]' will treat features C2_10010C_A and C2_10010C_B as homologs")
        # okay idea, but needs more thought before implementation
        parser_generate.add_argument("--homologs", metavar="*.homologs", type=str, default=None,
            help="Path to text file containing homologous features on the same \
            line, separated by TAB characters")
        #parser.add_argument("--pams", metavar="SEQ", nargs="+", type=str,
        #    default=["NGG"], help="Constrain finding only targets with these PAM sites")
        #parser.add_argument("--target_lengths", nargs=2, metavar=('MIN', 'MAX'),
        #    type=int, default=[17, 20],
        #    help="The length range of the 'target'/'spacer'/gRNA site")
        # Replacement for --pams and --target_lengths with this:
        parser_generate.add_argument("--motifs", metavar="MOTIF", nargs="+", type=str,
            default=["N{17}|N{3}>NGG"],
            help="Find only targets with these 'SPACER>PAM' motifs, written from \
            5' to 3'. '>' points toward PAM. IUPAC ambiguities accepted. '{a,b}' \
            are quantifiers. '/' is a sense strand cut, '\\' is an antisense strand \
            cut, and '|' is a double-strand cut. '.' is a base used for positional \
            information, but not enzymatic recognition. Be sure to enclose each \
            motif in quotes so your shell does not interpret STDIN/STDOUT redirection.")
        parser_generate.add_argument("--off_target_motifs", metavar="MOTIF", nargs="+", type=str,
            default=[],
            help="Defaults to the same as the on-target motif. Definition syntax is identical.")
        # Need to decide if construct inputs should be TSV, or FASTA
        # And whether or not there should be an upstream parameter separate from
        # a downstream one. or if they are the same, then what?
        parser_generate.add_argument("--constructs", metavar="*.fasta", nargs="+", type=str,
            default=[], help="The first sequence will be prepended, and the second \
            sequence will be appended to the generated spacer sequences to form \
            the construct sequences. It is useful to put the gRNA promotor as the \
            first sequence, and the scaffold sequence and terminator as the \
            second. Specify one FASTA file for each motif.")
        parser_generate.add_argument("--tag", metavar='TAG', type=str, default='ID',
            help="GFF tag with feature names. Examples: 'ID', 'Name', 'Gene', 'Parent', or 'locus_tag'")
        parser_generate.add_argument("--selection", metavar='FEATURE', nargs="+", type=str, default=None, #default=argparse.SUPPRESS, # '==SUPPRESS=='
            help="Select only certain features rather than all features in input GFF file.")
        parser_generate.add_argument("--ambiguities", type=str,
            choices=["exclusive", "discard", "disambiguate", "keep"],
            default="discard",
            help="How generated gRNAs should treat ambiguous bases: \
            exclusive - gRNAs will only be created for ambiguous locations; \
            discard - no gRNAs will be created where the FASTA has an ambiguous base; \
            disambiguate - gRNAs containing ambiguous bases will be converted to a set of non-ambiguous gRNAs; \
            keep - gRNAs can have ambiguous bases.")
        parser_generate.add_argument("--case", type=str, default="ignore",
            choices=["ignore", "upper-only", "lower-only", "mixed-lower", "mixed-upper", "mixed-only"],
            help="Restrict generation of gRNAs based on case of nucleotides in input FASTA: \
            ignore - keep all spacers; \
            upper-only - discard any potential spacer with a lower-case character in the genome; \
            lower-only - discard all potential spacers with upper-case characters; \
            mixed-lower - discard spacers that are all upper-case; \
            mixed-upper - discard spacers that are all lower-case; \
            mixed-only - only use spacers that have both lower- and upper-case characters.")
        #parser.add_argument("--min_contig_edge_distance", metavar="N", type=int, default=500,
        #    help="Minimum distance from contig edge a site can be found")
        parser_generate.add_argument("--features", metavar="FEATURE", type=str, nargs="+", default=["gene"],
            help="Features to design gRNAs against. Must exist in GFF file. \
            Examples: 'CDS', 'gene', 'mRNA', 'exon', 'intron', 'tRNA', 'rRNA'.\
            The special 'all' feature type will include all listed features.")
        #parser_generate.add_argument("--warning_features", metavar='FEATURE', nargs="+", type=str, default=['all'],
        #    help="GFF tags that will trigger a warning if they overlap with the \
        #    target feature. Examples: 'CDS', 'gene', 'mRNA', 'exon', 'intron', 'tRNA', 'rRNA'")
        parser_generate.add_argument("--excluded_features", metavar='FEATURE', nargs='+', type=str, default=[],
            help="Prevent feature expansion if it would overlap one of these features present in the GFF. \
            Examples: 'CDS', 'gene', 'mRNA', 'exon', 'intron', 'tRNA', 'rRNA'.\
            The special 'all' feature type will exclude all listed features \
            except those specified by the '--features' option.")
        parser_generate.add_argument("--dDNA_gDNA_ratio", metavar="N", type=int, default=1000,
            help="Ratio of donor DNA to genomic DNA for calculating off-target scores")
        #parser_generate.add_argument("--target_gc", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[25, 75],
        #    help="Generated gRNA spacers must have %%GC content between these values (excludes PAM motif)") # Moved to prefilter
        
        #  center_feature: --------HHHH[...............FEATURE.........TARGET]HHHH------------------------
        #                  -----HHHH[..................FEATURE.........TARGET...]HHHH--------------------- pad=3
        #   center_target: -----------------------HHHH[FEATURE.........TARGET................]HHHH--------
        #                  --------------------HHHH[...FEATURE.........TARGET...................]HHHH----- pad=3
        #     center_both: -----------------------HHHH[FEATURE.........TARGET]HHHH------------------------
        #                  --------------------HHHH[...FEATURE.........TARGET...]HHHH--------------------- pad=3
        # justify_feature: -----------------------HHHH[FEATURE.........TARGET]HHHH------------------------
        #                  -----------------------HHHH[FEATURE.........TARGET...]HHHH--------------------- pad=3
        #  justify_target: -----------------------HHHH[FEATURE.........TARGET]HHHH------------------------
        #                  --------------------HHHH[...FEATURE.........TARGET]HHHH------------------------ pad=3
        parser_generate.add_argument("--feature_expansion_method", type=str, default=None,
            choices=['center_feature', 'center_target', 'center_both', 'justify_feature', 'justify_target'],
            help="If a feature needs to be expanded to contain a gRNA target, \
            expand the feature such that either the feature, the target, or \
            both the feature and the target are in the center.")
        parser_generate.add_argument("--feature_expansion_lengths", nargs=2, metavar=("MIN", "MAX"), type=int, default=[100, 4000],
            help="If a feature needs to be expanded to contain a gRNA target (SPACER>PAM site), \
            the expanded feature must be within this range, inclusive.")
        parser_generate.add_argument("--feature_expansion_pad", metavar="N", type=int, default=0,
            help="If a feature needs to be expanded to contain a gRNA target, \
            and the expanded feature is within the permitted lengths, then \
            the expanded feature will have N number of nucleotides padded.")
        
        parser_generate.add_argument("--excise_upstream_homology", nargs=2, metavar=("MIN", "MAX"), type=int, default=[50,50],
            help="Range of homology lengths acceptable for knock-out dDNAs, inclusive.")
        parser_generate.add_argument("--excise_downstream_homology", nargs=2, metavar=("MIN", "MAX"), type=int, default=[47,50],
            help="Range of homology lengths acceptable for knock-out dDNAs, inclusive.")
        #
        #
        parser_generate.add_argument("--excise_donor_lengths", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[100, 100],
            help="Range of lengths acceptable for knock-out dDNAs, inclusive.")
        parser_generate.add_argument("--excise_insert_lengths", nargs=2, metavar=("MIN", "MAX"), type=int, default=[0,3],
            help="Range for inserted DNA lengths, inclusive (mintag). \
            If MIN < 0, then regions of dDNA homology (outside the feature) will be removed.")
        parser_generate.add_argument("--excise_feature_edge_distance", metavar="N", type=int, default=0,
            help="If positive, gRNAs won't target any nucleotides within this distance \
            from the edge of the feature. If negative, gRNAs will target nucleotides \
            this distance outside the feature.")
        #
        # May need to remove '--excise_upstream_feature_trim' and '--excise_downstream_feature_trim'
        # And replace with '--expand_feature {upstream,downstream,both}'
        parser_generate.add_argument("--excise_upstream_feature_trim", nargs=2, metavar=('MIN', 'MAX'),
            type=int, default=[0, 0], help="Between MIN and MAX number of nucleotides \
            upstream of the feature will be considered for knock-out when designing \
            donor DNA.")
        parser_generate.add_argument("--excise_downstream_feature_trim", nargs=2, metavar=("MIN", "MAX"),
            type=int, default=[0, 0], help="Between MIN and MAX number of nucleotides \
            downstream of the feature will be considered for knock-out when designing \
            donor DNA.")
        #parser_generate.add_argument("--expand_feature", type=str, default="both",
        #    choices=["none", "upstream", "downstream", "both"],
        #    help="Sequence immediately up/down-stream of intended feature may \
        #    be knocked-out in order to generate good mintag sites.")
        #
        #
        #parser.add_argument("--min_donor_insertions", metavar="N", type=int, default=2,
        #    help="The uniqueness of final donor DNA compared to the rest of the genome")
        #parser.add_argument("--min_donor_deletions", metavar="N", type=int, default=2,
        #    help="The uniqueness of final donor DNA compared to the rest of the genome")
        #parser.add_argument("--min_donor_substitutions", metavar="N", type=int, default=2,
        #    help="The uniqueness of final donor DNA compared to the rest of the genome")
        #parser.add_argument("--min_donor_errors", metavar="N", type=int, default=3,
        #    help="The uniqueness of final donor DNA compared to the rest of the genome")
        parser_generate.add_argument("--revert_upstream_homology", nargs=2, metavar=("MIN", "MAX"), type=int, default=[200,400],
            help="Range of homology lengths acceptable for knock-in dDNAs, inclusive.")
        parser_generate.add_argument("--revert_downstream_homology", nargs=2, metavar=("MIN", "MAX"), type=int, default=[200,400],
            help="Range of homology lengths acceptable for knock-in dDNAs, inclusive.")
        #parser_generate.add_argument("--revert_donor_lengths", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[0, 100000],
        #    help="Range of lengths acceptable for knock-in dDNAs.")
    #    parser.add_argument("--min_donor_distance", metavar="N", type=int, default=36,
    #        help="The minimum distance in bp a difference can exist from the edge of donor DNA") # homology with genome
        #parser_generate.add_argument("--max_consecutive_ts", metavar="N", type=int, default=4, # Moved to prefilter algorithm
        #    help="The maximum number of Ts allowed in generated gRNA sequences.")
        parser_generate.add_argument("--max_number_sequences_reported", metavar="N", type=int, default=5,
            help="The maximum number of sequences to report for each step.")
        parser_generate.add_argument("--min_weight_reported", metavar="N", type=float, default=0.01,
            help="Only gRNA-dDNA pairs with at least this much weight will be reported.")
        # program currently will only search 'both' strands
        #parser.add_argument("--strands", type=str, choices=["+", "-", "both"], default="both",
        #    help="Strands to search for gRNAs")
        
        # Add command line arguments for the additional hard constraints:
        #  Only report potential targets that have no off targets with mismatches within 8, 12, N nt from 3' end
        parser_generate.add_argument("--processors", metavar="N", type=int, default=(os.cpu_count() or 1),
            help="Number of processors to use when performing pairwise sequence alignments.")
        
        aligner_choices = [x.name for x in aligners.aligners]
        parser_generate.add_argument("--aligner", type=str, choices=aligner_choices, default='bowtie2',
            help="Program to calculate pairwise alignments. Please note that the 'addtag' internal aligner is very slow.")
        
        oligo_choices = [x.name for x in oligos.oligos]
        parser_generate.add_argument("--oligo", type=str, choices=oligo_choices, default='UNAFold',
            help="Program to perform thermodynamic calculations.")
        
        prefilter_choices = [C.name for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.prefilter]
        parser_generate.add_argument("--prefilters", nargs='+', type=str,
            choices=prefilter_choices, default=['GC', 'PolyT'],
            help="Specific algorithms for determining gRNA goodness.")
        
        off_target_choices = [C.name for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.off_target]
        parser_generate.add_argument("--offtargetfilters", nargs='+', type=str,
            choices=off_target_choices, default=['CFD', 'Hsu-Zhang'],
            help="Specific algorithms for determining gRNA goodness.")
        
        on_target_choices = [C.name for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.on_target]
        parser_generate.add_argument("--ontargetfilters", nargs='+', type=str,
            choices=on_target_choices, default=['Azimuth'],
            help="Specific algorithms for determining gRNA goodness.")
        
        postfilter_choices = [C.name for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.postfilter]
        parser_generate.add_argument("--postfilters", nargs='+', type=str,
            choices=postfilter_choices, default=['Errors'],
            help="Specific algorithms for determining gRNA goodness.")
        
        parser_generate.add_argument("--primer_scan_limit", metavar="N", type=int, default=2*60,
            help="Amount of time (in seconds) to limit each primer scan.")
        
        parser_generate.add_argument("--primer_pair_limit", metavar="N", type=int, default=5*60,
            help="Amount of time (in seconds) to limit primer pairings.")
        
        # Temporary stopgap to make sure the calculations don't take too long
        #parser_generate.add_argument("--max_time", metavar="SECONDS", type=float, default=60,
        #    help="Maximum amount of time, in seconds, for each feature to spend calculating dDNAs.")
        
        parser_generate.add_argument("--bartag_number", metavar="N", type=int, default=1,
            help="Number of bartags per locus to generate. \
            For efficiency's sake, a maximum 'features × bartags = 200' is enforced.")
        
        parser_generate.add_argument("--bartag_motif", metavar="MOTIF", type=str, default='N{10}',
            help="Structure of nucleotides that should be generated. \
            Longer oligonucleotide lengths take longer to calculate.")
        
        parser_generate.add_argument("--bartag_distance", metavar="N", type=int, default=3,
            help="Minimum required edit distance between each bartag.")
        
        parser_generate.add_argument("--flanktags", nargs='+', action=ValidateFlanktags, type=str,
            metavar=('{*.fasta,uniform,specific}', '{single,pair}'), default=None,
            help="If 'uniform', then all features will share the same flanktags. \
            If 'specific', then each feature will have its own flanktags. \
            If '*.fasta', then flanktags will be taken from an input FASTA file with either a single or paired flanktags. \
            If 'single', then design a single oligo to serve as both the uptag and dntag. \
            If 'pair', then design uptag different to dntag.")
        
        #--flanktags myfile.fasta
        #--flanktags uniform single
        #--flanktags uniform pair
        #--flanktags feature single
        #--flanktags feature pair
        
        # By default, unlabeled flanktags won't be strand specific, but if they are labeled, then they will be strand-specific
        # The flanktag length will be 18-25 nt, preferring 20
        
        # Splitting up "generate" into specific sub-tasks
        parser_generate.add_argument("--ko-gRNA", action='store_true', default=False,
            help="Design gRNAs to target features in genome")
        
        #parser_generate.add_argument("--ko-dDNA", type=str, default=None,
        #    choices=['mintag', 'addtag', 'unitag', 'bartag'],
        #    help="'mintag' are unique us/i/ds junction targets specific to each feature. \
        #    'addtag' are unique targets for each feature. \
        #    'unitag' is a single, invariant target for ALL features. \
        #    'bartag' are unique barcodes for each feature (does not guarantee targets).")
        
        parser_generate.add_argument("--ko-dDNA", type=str, action=ValidateKodDNA,
            metavar='{*.fasta,mintag,addtag,unitag,bartag}', default=None,
            help="'*.fasta' is a FASTA file containing user-specified sequences. \
            'mintag' are unique us/i/ds junction targets specific to each feature. \
            'addtag' are unique targets for each feature. \
            'unitag' is a single, invariant target for ALL features. \
            'bartag' are unique barcodes for each feature (does not guarantee targets).")
        
        parser_generate.add_argument("--ki-gRNA", action='store_true', default=False,
            help="Design gRNAs to target the ko-dDNA. \
            Defaults to True if '--ko-dDNA mintag' is specified.")
        
        parser_generate.add_argument("--ki-dDNA", nargs='?', type=str, default=None,
            metavar='*.fasta', const=True, # If no command-line argument follows, the value of const will be assumed instead.
            help="If no file is specified, then design wild type dDNA. \
            If a FASTA file is specified, then design dDNA to replace features \
            with the sequences in this FASTA file. The primary sequence header \
            should either be the TAG identifier for the feature to replace \
            (specified by the '--tag' option) or the gene name (specified \
            within the '--homologs' file).")
        #    help="If no file is specified, then design wild type dDNA. \
        #    If a FASTA file is specified, then design dDNA to replace features \
        #    with the sequences in this FASTA file. The sequence header should \
        #    have a 'TAG=FEATURE' field where 'TAG' is the tag specified with \
        #    '--tag' option, and 'FEATURE' corresponds to the value of that \
        #    'TAG' within the input GFF.")
        
        # subparsers
        #  addtag generate
        #   --ko-gRNA                              Design gRNAs to target features in genome
        #   --ko-dDNA mintag|addtag|unitag|bartag  Design knock-out dDNAs for target features in genome
        #                                          (these will be constrained to dDNAs targetable by ki-gRNAs)
        #             mintag                          Designs unique us/i/ds junction target specific to each feature
        #                    addtag                   Designs unique target for each feature
        #                           unitag            Designs a single/uniform, unique target for ALL features
        #                                  bartag     Designs unique barcode for each feature
        #   --bartag_length INT                    Specify length of sigtag barcode in nt (Or should this be --sigtag_motif NNNNNN?)
        #   --flanktags generate                   Design identical uptag/dntag primer sequences to addtag/sigtag/unitag for each site/feature
        #   --flanktag_length INT INT              Specify length of generated uptag and dntag sequences in nt
        #   --flanktags input FASTA                Use file specifying uptag/dntag primers sequences for addtag/sigtag for all sites
        #                                          >seq1 ID=feature1 flanktag=uptag
        #                                          NNNNNNNNNNNNNNNNNNNN
        #                                          >seq2 ID=feature1 dntag
        #                                          NNNNNNNNNNNNNNNNNNNN
        #   
        #   --ki-gRNA     Design knock-in gRNAs that target the ko-dDNA
        #   --ki-dDNA wt        Design the dDNA for knocking the wild type feature back in
        #   --ki-dDNA FASTA     If a FASTA file is specified, then use mutants for these indicated features
        #                       >seq1 ID=feature1
        #                       NNNNNNNNNNNNNNNNNNN
        #                       The entire feature will be replaced by this sequence
        
        return parser_generate
    
    def _parser_confirm(self, subparsers):
        ''' "confirm" parser '''
        __confirm_description__ = """\
description:
  Design primers for confirming whether each step of genome engineering is
  successful or not. This does not design multiplexable primers.
"""
        __confirm_help__ = "Design primers for confirming whether each step of genome engineering is successful or not."
        parser_confirm = subparsers.add_parser('confirm',
            description=__confirm_description__,
            #epilog=__confirm_epilog__,
            formatter_class=CustomHelpFormatter,
            help=__confirm_help__
        )
        parser_confirm.set_defaults(func=self._confirm)
        
        # Change the help text of the "-h" flag
        parser_confirm._actions[0].help='Show this help message and exit.'
        
        # Special version action optional argument
        parser_confirm.add_argument("-v", "--version", action='version',
            help="Show program's version number and exit.",
            version='{__program__} {__version__} (revision {__revision__})'.format(**globals()))
        
        # Add mandatory arguments
        # Add required arguments
        required_group = parser_confirm.add_argument_group('required arguments')
        required_group.add_argument("--fasta", required=True, nargs="+", metavar="*.fasta", type=str,
            help="FASTA files with contigs to find unique gRNA sites. Input FASTA \
            must not be compressed. User can decide whether ambiguous bases can \
            be chosen for gRNA sites. All FASTA sequences should have unique \
            primary headers (everything between the '>' symbol and the first \
            whitespace should be unique). You should include FASTA of the genome \
            and any plasmids.")
        
        required_group.add_argument("--dDNAs", required=True, nargs="+", metavar="*.fasta", type=str,
            help="A dDNA FASTA file for each subsequent CRISPR/Cas transformation. \
            Typically, the first is KO, and the second is KI. However, any number \
            of serial genome engineering experiments can be specified.")
        
        required_group.add_argument("--folder", required=True, metavar="FOLDER",
            type=str, help="Path of folder to store generated files.")
        
        # Add optional arguments
        parser_confirm.add_argument("--processors", metavar="N", type=int, default=(os.cpu_count() or 1),
            help="Number of processors to use when performing pairwise sequence alignments.")
        
        aligner_choices = [x.name for x in aligners.aligners]
        parser_confirm.add_argument("--aligner", type=str, choices=aligner_choices, default='blastn',
            help="Program to calculate pairwise alignments. Please note that the 'addtag' internal aligner is very slow.")
        
        parser_confirm.add_argument("--number_pcr_conditions", metavar="N", type=int, default=None,
            help="Number of PCR conditions to develop primers for. All amplicons \
            within each condition will have similar size (nt) and melting temperatures. \
            If unspecified, will default to the number of target features.")
        
        parser_confirm.add_argument("--primer_scan_limit", metavar="N", type=int, default=2*60,
            help="Number of seconds to limit each primer scan.")
        
        parser_confirm.add_argument("--primer_pair_limit", metavar="N", type=int, default=5*60,
            help="Amount of time (in seconds) to limit primer pairings.")
        
        parser_confirm.add_argument("--primers", nargs="+", metavar="*.fasta", type=str, default=[],
            help="A FASTA file for each round containing primer sequences that \
            must be used for that round. Usually, these correspond to flanktags.")
        
        oligo_choices = [x.name for x in oligos.oligos]
        parser_confirm.add_argument("--oligo", type=str, choices=oligo_choices, default='UNAFold',
            help="Program to perform thermodynamic calculations.")
        
        parser_confirm.add_argument("--max_evalue", metavar="N", type=float, default=0.001,
            help="The maximum accepted alignment E-value to consider a recombination.")
        
        parser_confirm.add_argument("--min_length", metavar="N", type=int, default=35,
            help="The minimum accepted alignment length to consider a recombination.")
        
        parser_confirm.add_argument("--skip_round", metavar="N", nargs="+", type=int, default=[],
            help="Skip primer calculations for these rounds.")
        
        # Nucleotide matching stuff
        #  - number errors (for fuzzy regex)
        
        # PCR conditions:
        #  - primer_size (min, max)
        #  - monovalent_cation_concentration
        #  - divalent_cation_concentration
        #  - amplicon_size (min, max)
        #  - template_concentration
        #  - dNTP_concentration
        #  - temperature (25 usually)
        #  - max_tm_difference
        #  - melting_temperature (min, max)
        #  - max_3prime_homology_length
        #  - min_delta_g
        #  - gc (min, max)
        #  - gc_clamp_length (min, max)
        #  - max_run_length
        
        return parser_confirm
    
    def parse_arguments(self):
        '''
        Create parser object, populate it with options, then parse command line
        input.
        '''
        
        parser, subparsers = self._parser_general()
        parser_glossary = self._parser_glossary(subparsers)
        parser_motifs = self._parser_motifs(subparsers)
        parser_algorithms = self._parser_algorithms(subparsers)
        parser_aligners = self._parser_aligners(subparsers)
        parser_oligos = self._parser_oligos(subparsers)
        parser_search = self._parser_search(subparsers)
        parser_feature = self._parser_feature(subparsers)
        parser_extract = self._parser_extract(subparsers)
        parser_evaluate = self._parser_evaluate(subparsers)
        parser_generate = self._parser_generate(subparsers)
        parser_confirm = self._parser_confirm(subparsers)
        
        
        # Add special arguments for additional help messages --> These have been deprecated
        #parser.add_argument("-s", "--show", nargs=0, action=ValidateShowMotifs, default=argparse.SUPPRESS,
        #    help="Show list of common RGN motifs, then exit.")
        #parser.add_argument("-g", "--glossary", nargs=0, action=ValidateShowGlossary, default=argparse.SUPPRESS,
        #    help="Show glossary of common CRISPR/RGN terms, then exit.")
        
        # Add optional arguments --> Deprecated
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
        
        # Return the parsed arguments
        return args
    
    def process_general_arguments(self, args):
        """Perform ubiquitous argument parsing things"""
        if hasattr(args, 'folder'):
            # Create the project directory if it doesn't exist
            os.makedirs(args.folder, exist_ok=True)
        
            # Create the logger, and have it write to 'folder/log.txt'
            logging.basicConfig(filename=os.path.join(args.folder, 'log.txt'), level=logging.INFO, format='%(message)s') # format='%(levelname)s %(asctime)s: %(message)s'
        
        if hasattr(args, 'aligner'):
            # Add 'args.selected_aligner' to hold the actual aligner object
            for a in aligners.aligners:
                if (a.name == args.aligner):
                    args.selected_aligner = a
                    break
        
        if hasattr(args, 'oligo'):
            # Add 'args.selected_oligo' to hold the actual oligo object
            for o in oligos.oligos:
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
    
    def rolling_picker(self):
        """
        Rolling picker algorithm for multiplex gRNA/dDNA.
        
        Once all candidate target sequences are fully annotated and ranked, the
        sgRNA designer cycles through the list of candidates, attempting to pick
        sequences in order to achieve the best final set of sequences.
        
        We select all best-ranked gRNA, and if there is a conflict, then the
        worse one is replaced with a different candidate.
        
        After each round of picking, constraints are somewhat relaxed until
        a usable set of gRNA are found.
        
        """
        
        pass
    
    def merge_features(self, features):
        """Combine overlapping features?"""
        return features
    
    def get_exTarget_homologs(self, homologs):
        """Get ExcisionTarget objects for each homologous feature group"""
        ext_dict = {}
        if homologs:
            groups = set()
            for feature_name, f in Feature.features.items():
                #groups.add(tuple(sorted(homologs[feature_name])))
                groups.add(tuple(sorted(homologs[f.get_expand_parent().name])))
            for g in groups:
                ext_dict[g] = set()
                for name, obj in ExcisionTarget.indices.items():
                    obj_features = set(utils.flatten(x[0].split(',') for x in obj.locations))
                    if (len(obj_features.intersection(g)) == len(g)):
                        ext_dict[g].add(obj)
        return ext_dict
    
    def get_exTarget_allele_specific(self):
        """Gets allele-specific ExcisionTarget objects"""
        ext_dict = {}
        for feature_name, f in Feature.features.items():
        #for f in features:
            ext_dict[feature_name] = set()
            for name, obj in ExcisionTarget.indices.items():
                obj_features = set(utils.flatten(x[0].split(',') for x in obj.locations))
                if ((len(obj_features) == 1) and (feature_name in obj_features)):
                    ext_dict[feature_name].add(obj) # Store with the feature as key
        return ext_dict
    
    def get_reTarget_homologs(self, homologs):
        """Get ReversionTarget objects for each homologous feature group"""
        ret_dict2 = {}
        
        if homologs:
            groups = set() # A set of all the parent (non-derived) feature names
            for feature_name, f in Feature.features.items():
                #contig, start, end, strand = features[feature]
                # Convert to tuple
                #groups.add(tuple(sorted(homologs[feature_name])))
                groups.add(tuple(sorted(homologs[f.get_expand_parent().name])))
            # groups = {('F1_A', 'F1_B'), ('F2_A', 'F2_B')} # A set of tuples
            
            for g in groups:
                ret_dict2[g] = set()
                # Get all ReversionTargets that have all feature of this group in its location information
                
                for name, obj in ReversionTarget.indices.items():
                    #obj_features = set(x[0] for x in obj.locations) # these are comma-separated
                    obj_features = set(utils.flatten(x[0].split(',') for x in obj.locations))
                    
                    if (len(obj_features.intersection(g)) == len(g)):
                        ret_dict2[g].add(obj) # Store with the tuple form of the group as key
                
                #print('ret_dict2', len(g), g)
                #for obj in sorted(ret_dict2[g], key=lambda x: int(x.name.split('-')[1])):
                #    print(' ', obj)
        
        # key = ('F1_A', 'F1_B') # tuple of homologous features
        # value = [ReversionTarget(), ReversionTarget(), ...] # list of ReversionTarget objects
        return ret_dict2
    
    def get_reTarget_allele(self, features):
        """Gets all ReversionTarget objects for each feature"""
        ret_dict = {}
        for f in features:
            ret_dict[f] = set()
            for name, obj in ReversionTarget.indices.items():
                obj_features = set(utils.flatten(x[0].split(',') for x in obj.locations))
                if (f in obj_features):
                    ret_dict[f].add(obj) # Store with the feature as key
        return ret_dict
    
    def get_reTarget_allele_specific(self):
        """Gets allele-specific ReversionTarget objects"""
        ret_dict = {}
        for feature_name, f in Feature.features.items():
        #for f in features:
            ret_dict[feature_name] = set()
            for name, obj in ReversionTarget.indices.items():
                obj_features = set(utils.flatten(x[0].split(',') for x in obj.locations))
                if ((len(obj_features) == 1) and (feature_name in obj_features)):
                    ret_dict[feature_name].add(obj) # Store with the feature as key
        return ret_dict
    
    def rank_donors(self, donor_list):
        """Uses a heuristic to find the best Donor in the list"""
        # If the best ReversionTarget has multiple locations, then
        # choose the best ExcisionDonor location:
        #    1) minimize the distance between x and y: w..x:mAT:y..z
        #    2) minimize the length of mAT
        #    3) report all ties (don't break ties)
        
        rank_list = []
        for d in donor_list: # these should all be ExcisionDonor objects
            gaps = []
            for l in d.locations:
                gaps.append(len(l[4]) + l[5][0]-l[3][1])
            rank_list.append((min(gaps), d))
        
        return sorted(rank_list, key=lambda x: x[0]) # smallest insert size/gap length will be first
    
    def rank_targets(self, target_list):
        """Uses a heuristic to find the best Target in the list"""
        
        # Hsu-Zhang off-target score should be >95
        #   if Hsu-Zhange < 95, score should drop quickly
        # Azimuth on-target score should be >60
        #   if Azimuth < 60, score should drop quickly
        # CFD off-target score should be >50
        #   if CFD < 50, score should drop quickly
        
        #rank_azimuth = lambda x: 1/(1+1.17**(50-x))
        #rank_hsuzhang = lambda x: 1/(1+1.8**(90-x))
        #rank_cfd = lambda x: 1/(1+1.2**(40-x))
        #rank = lambda x, y, z: rank_azimuth(x)*rank_hsuzhang(y)*rank_cfd(z)
        
        # Returns list of targets, sorted such that the one with the highest
        # aggregate score is the 0th index
        rank_list = []
        for t in target_list:
            #rank_list.append((rank(t.score['Azimuth'], t.off_targets['Hsu-Zhang'], t.off_targets['CFD']), t))
            
            rank = 1.0
            for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms:
                if C.off_target:
                    rank *= C.weight(t.off_targets[C.name])
                else:
                    rank *= C.weight(t.score[C.name])
            rank_list.append((rank, t))
        
        #return sorted(target_list, key=lambda x: rank(x.score['Azimuth'], x.off_targets['Hsu-Zhang'], x.off_targets['CFD']), reverse=True)
        return sorted(rank_list, key=lambda x: x[0], reverse=True)
    
    def get_best_table(self, args, homologs, feature2gene):
        """
        Identify and print the best spacers and dDNAs for each feature, given
        each mAT insert size, and us/ds trim length.
        Thus, the user can see the best spacer for any given combination of these.
        """
        
        ########## Results table for the knock-out dDNAs and the efficiency of their gRNA for cutting them (ExcisionDonor+ReversionTarget or ko-dDNA+ki-gRNA) ##########
        # Add a column for "WARNING", that tells the user if the target feature overlaps with another feature
        #  it will just list the feature IDs that overlap:
        #   C4_03620C_A-T,C4_03620C_A-T-E1
        
        # Only include columns for algorithms whose weight does not equal 1 uniformly (alternatively, whose classes override the default 'weight()' method)
        #algorithms.weighted_algorithms
        
        #header = ['gene', 'features', 'insert', 'mAT', 'translations', '(us, ds) trim', 'weight', 'OT:Hsu-Zhang', 'OT:CFD', 'Azimuth', 'reTarget name', 'reTarget sequence', 'ExDonors']
        # Add columns for:
        #  contig:start..end
        #   contig_A:40221..40241,contig_B:40155..40175
        #  feature:start..end
        #   CR_02630C_A:221..241,CR_02630C_B:224..244
        #  feature,feature:contig:strand:start..end
        #   C1_06280C_A,C1_06280C_B:exDonor-222:-:43..66
        #  strand
        #   +     or     -
        #
        # 'hairpin ΔG', 'homodimer ΔG', 'heterodimer ΔG'
        # locations0s --> 'features:contig:strand:start..end'
        header = ['gene', 'features', 'contig:strand:start..end', 'us-trim:mAT:ds-trim', 'translations', 'weight', 'OT:Hsu-Zhang', 'OT:CFD', 'Azimuth', 'reTarget name', 'reTarget sequence', 'exDonors', 'warning']
        print('\t'.join(header))
        
        # Get best ReversionTargets by calculating their weights, and also getting
        # the ExcisionDonors that correspond to the ReversionTarget
        ret_dict2 = self.get_reTarget_homologs(homologs)
        for feature_homologs in sorted(ret_dict2):
            # Get these two things which are common to all the top hits
            gene = feature2gene[feature_homologs[0]] # Get the gene name
            csfeatures = ','.join(feature_homologs) # Add the features as comma-separated list
            
            features_pos = '' #features[gene]
            
            outputs = {}
            
            # Print the top N for each insert size and trim
            for weight, obj in self.rank_targets(ret_dict2[feature_homologs]):
                othz = round(obj.off_targets['Hsu-Zhang'], 2)
                otcfd = round(obj.off_targets['CFD'], 2)
                azimuth = round(obj.score['Azimuth'], 2)
                
                # Get the ExcisionDonor objects for this ki-spacer, and weigh them
                rds = self.rank_donors(obj.get_donors())
                # filter out all but the top-weighted ones
                rds = [x for x in rds if (x[0] == rds[0][0])]
                
                exdonors = ','.join(map(lambda x: x[1].name, rds))
                
                key0 = sorted(set(utils.flatten(map(lambda x: x[1].get_inserts_and_trims(), rds))))
                key0s = ','.join('{}:{}:{}'.format(x[1], x[0], x[2]) for x in key0)
                key1 = sorted(set([(len(x[0]), x[1], x[2]) for x in key0])) # replace mAT with length
                key1s = ','.join(map(str, key1))
                #locations0 = sorted(set(utils.flatten([xo[1].format_location(x) for x in xo[1].locations] for xo in rds)))
                #locations0s= ','.join(locations0)
                translations = None
                
                sline = [gene, csfeatures, features_pos, key0s, translations, weight, othz, otcfd, azimuth, obj.name, obj.spacer+'|'+obj.pam, exdonors]
                
                if (weight >= args.min_weight_reported):
                    if (len(outputs.get(key1s, [])) < args.max_number_sequences_reported):
                        outputs.setdefault(key1s, []).append(sline)
            
            for k in sorted(outputs): # sort by insert/trims
            #for k in sorted(outputs, key=lambda x: outputs[x][4], reverse=True): # sort by weight
                for sline in outputs[k]:
                    print('\t'.join(map(str, sline)))
            
            if (len(outputs) == 0):
                logging.info('No spacers for targeting knocked-out ' + csfeatures)
        
        # Print a table of allele-specific knock-in spacer and dDNA pairs
        ret_dict = self.get_reTarget_allele_specific()
        for feature in sorted(ret_dict):
            gene = feature2gene[feature] # Get the gene name
            outputs = {}
            
            for weight, obj in self.rank_targets(ret_dict[feature]):
                othz = round(obj.off_targets['Hsu-Zhang'], 2)
                otcfd = round(obj.off_targets['CFD'], 2)
                azimuth = round(obj.score['Azimuth'], 2)
                
                rds = self.rank_donors(obj.get_donors())
                rds = [x for x in rds if (x[0] == rds[0][0])]
                
                exdonors = ','.join(map(lambda x: x[1].name, rds))
                
                key0 = sorted(set(utils.flatten(map(lambda x: x[1].get_inserts_and_trims(), rds))))
                key0s = ','.join('{}:{}:{}'.format(x[1], x[0], x[2]) for x in key0)
                key1 = sorted(set([(len(x[0]), x[1], x[2]) for x in key0])) # replace mAT with length
                key1s = ','.join(map(str, key1))
                translations = None
                
                sline = [gene, feature, key0s, translations, weight, othz, otcfd, azimuth, obj.name, obj.spacer+'>'+obj.pam, exdonors] # should automatically determine the direction of the SPACER>PAM based on the motif
                
                if (weight >= args.min_weight_reported):
                    if (len(outputs.get(key1s, [])) < args.max_number_sequences_reported):
                        outputs.setdefault(key1s, []).append(sline)
            
            for k in sorted(outputs):
                for sline in outputs[k]:
                    print('\t'.join(map(str, sline)))
            
            # If there are no allele-specific records, then nothing is printed
            if (len(outputs) == 0):
                logging.info('No allele-specific spacers for targeting knocked-out ' + feature)
        
        ########## Results table for cutting the wild-type genome (ExcisionTarget or ko-gRNA) ##########
        
        # to add to header: contig:strand:start..end
        header = ['gene', 'features', 'weight', 'OT:Hsu-Zhang', 'OT:CFD', 'Azimuth', 'exTarget name', 'exTarget sequence', 'reDonors']
        print('\t'.join(header))
        
        # Print the best homozygous ExcisionTargets for each feature set
        ext_dict2 = self.get_exTarget_homologs(homologs)
        for feature_homologs in sorted(ext_dict2):
            gene = feature2gene[feature_homologs[0]] # Get the gene name
            csfeatures = ','.join(feature_homologs) # Add the features as comma-separated list
            
            red_list = set()
            for name, obj in ReversionDonor.indices.items():
                obj_features = set(obj.get_location_features())
                if (len(obj_features.intersection(feature_homologs)) > 0):
                    red_list.add(obj)
            red_list = sorted(red_list, key=lambda x: int(x.name.split('-')[1]))
            
            outputs = []
            for weight, obj in self.rank_targets(ext_dict2[feature_homologs]):
                othz = round(obj.off_targets['Hsu-Zhang'], 2)
                otcfd = round(obj.off_targets['CFD'], 2)
                azimuth = round(obj.score['Azimuth'], 2)
                
                redonors = ','.join(x.name for x in red_list)
                
                sline = [gene, csfeatures, weight, othz, otcfd, azimuth, obj.name, obj.spacer+'>'+obj.pam, redonors] # should automatically determine the direction of the SPACER>PAM based on the motif
                
                if (weight >= args.min_weight_reported):
                    if (len(outputs) < args.max_number_sequences_reported):
                        outputs.append(sline)
            
            for sline in outputs:
                print('\t'.join(map(str, sline)))
            
            if (len(outputs) == 0):
                logging.info('No spacers for targeting ' + csfeatures)
        
        # Print a table of allele-specific knock-out spacer and dDNA pairs
        ext_dict = self.get_exTarget_allele_specific()
        for feature in sorted(ext_dict):
            gene = feature2gene[feature] # Get the gene name
            
            red_list = set()
            for name, obj in ReversionDonor.indices.items():
                if feature in obj.get_location_features():
                    red_list.add(obj)
            red_list = sorted(red_list, key=lambda x: int(x.name.split('-')[1]))
            
            outputs = []
            
            for weight, obj in self.rank_targets(ext_dict[feature]):
                othz = round(obj.off_targets['Hsu-Zhang'], 2)
                otcfd = round(obj.off_targets['CFD'], 2)
                azimuth = round(obj.score['Azimuth'], 2)
                
                redonors = ','.join(x.name for x in red_list)
                
                sline = [gene, feature, weight, othz, otcfd, azimuth, obj.name, obj.spacer+'|'+obj.pam, redonors]
                
                if (weight >= args.min_weight_reported):
                    if (len(outputs) < args.max_number_sequences_reported):
                        outputs.append(sline)
            
            for sline in outputs:
                print('\t'.join(map(str, sline)))
            
            # If there are no allele-specific records, then nothing is printed
            if (len(outputs) == 0):
                logging.info('No allele-specific spacers for targeting knocked-out ' + feature)
        
    
    def log_results(self, args, homologs, n=None):
        """Function that prints to the log file the best spacers and dDNAs for each feature"""
        
        logging.info('Log of best results...')
        
        # Print best ReversionTargets calculated and their corresponding ExcisionDonors
        logging.info("Best 'ReversionTarget's calculated and their corresponding 'ExcisionDonor's...")
        ret_dict2 = self.get_reTarget_homologs(homologs)
        for k in ret_dict2:
            logging.info(str(k) + ' ' + str(len(ret_dict2[k])))
            
            # Get the value for N
            if (n == None):
                display_num = len(ret_dict2[k])
            else:
                display_num = n
            
            # Print the top N
            for rank, obj in self.rank_targets(ret_dict2[k])[:display_num]:
                logging.info(' ' + str(rank) + ' ' + str(obj))
                # Get the ExcisionDonor objects for this ki-spacer, and rank them
                rds = self.rank_donors(obj.get_donors())
                # filter out all but the top-ranked ones
                rds = [x for x in rds if (x[0] == rds[0][0])]
                for gap, exd_obj in rds:
                    logging.info('   ' + str(gap) + ' ' + str(exd_obj.get_trims()) + ' ' + str(exd_obj))
        
        # Print best ExcisionTargets (not necessarily homozygous) for each feature
        logging.info("Best 'ExcisionTarget's for each feature...")
        # and the ReversionDonor
        for feature in sorted(Feature.features):
            logging.info(feature)
            et_list = []
            for name, obj in ExcisionTarget.indices.items():
                if feature in obj.get_location_features():
                    et_list.append(obj)
            
            # Get the value for N
            if (n == None):
                display_num = len(et_list[k])
            else:
                display_num = n
            
            # Print the top N
            for rank, obj in self.rank_targets(et_list)[:display_num]:
                logging.info('  ' + str(rank) + ' ' + str(obj))
            
            red_list = []
            for name, obj in ReversionDonor.indices.items():
                if feature in obj.get_location_features():
                    red_list.append(obj)
            for obj in red_list:
                logging.info('  ' + str(obj))
    
    def old_get_best(self, args, features, contigs):
        #exd_dict = {}
        red_dict = {}
        ext_dict = {}
        ret_dict = {}
        
        # Non-exclusively separate all instances into feature groups
        for feature in features:
            red_dict[feature] = []
            for name, obj in ReversionDonor.indices.items():
                if feature in obj.get_location_features():
                    red_dict[feature].append(obj)
            
            ext_dict[feature] = []
            for name, obj in ExcisionTarget.indices.items():
                if feature in obj.get_location_features():
                    ext_dict[feature].append(obj)
            
            ret_dict[feature] = set()
            for name, obj in ReversionTarget.indices.items():
                for c in obj.get_contigs():
                    if feature in ExcisionDonor.indices[c].get_location_features():
                        ret_dict[feature].add(obj)
        
        # Find the best instances for ko/ki
        for feature in sorted(features):
            # Find the best ReversionTarget
            on_target_sorted = sorted(ret_dict[feature], key=lambda x: x.score['Azimuth'], reverse=True)
            ret_best = None
            for i in range(len(on_target_sorted)):
                if (on_target_sorted[i].off_targets['Hsu-Zhang'] >= 90):
                    ret_best = on_target_sorted[i]
                    break
            #if not ret_best:
            #    print(len(on_target_sorted))
            #    for tmp in on_target_sorted:
            #        print(' ', tmp)
            #if not ret_best: # this could fail if there are no sites
            #    off_target_sorted = sorted(ret_dict[feature], key=lambda x: x.score['Hsu-Zhang'], reverse=True)
            #    ret_best = off_target_sorted[0]
            
            # Find the best ExcisionTarget
            on_target_sorted = sorted(ext_dict[feature], key=lambda x: x.score['Azimuth'], reverse=True)
            ext_best = None
            for i in range(len(on_target_sorted)):
                if (on_target_sorted[i].off_targets['Hsu-Zhang'] >= 95):
                    ext_best = on_target_sorted[i]
                    break
            
            # Find the ExcisionDonors that correspond with the best ReversionTarget
            exd_best = []
            if ret_best:
                exd_best = [ExcisionDonor.indices[x] for x in ret_best.get_contigs()]
            
            # Find the ReversionDonor
            red_best = red_dict[feature]
            
            logging.info("")
            logging.info("  feature = " + feature)
            logging.info("ko spacer = " + str(ext_best))
            for d in exd_best:
                logging.info("  ko dDNA = " + str(d))
            logging.info("ki spacer = " + str(ret_best))
            for d in red_best:
                logging.info("  ki dDNA = " + str(d))

def target_filter(sequence, target, pam, upstream, downstream, args):
    '''
    Filters the candidate gRNA sequence based on the following criteria:
     1) case: ignore, upper-only, lower-only, mixed-lower, mixed-upper, mixed-only
     2) ambiguous character expansion: exclusive, discard, keep, disambiguate
     3) SPACER>PAM check using regex (following disambiguation expansion)
     4) Prefilters: maximum consecutive Ts, %GC
    
    Returns list of validated sequences
    '''
    
    seq = sequence
    
    # Check the case of the potential gRNA sequence
    if (args.case == "upper-only"):
        if regex.search('[a-z]', seq):
            return [] # Reject this sequence because it has lower-case characters
    elif (args.case == "lower-only"):
        if regex.search('[A-Z]', seq):
            return [] # Reject this sequence because it has upper-case characters
    elif (args.case == "mixed-lower"):
        if not regex.search('[a-z]', seq):
            return []
    elif (args.case == "mixed-upper"):
        if not regex.search('[A-Z]', seq):
            return []
    elif (args.case == "mixed-only"):
        if not (regex.search('[a-z]', seq) and regex.search('[A-Z]', seq)):
            return [] # Reject this sequence because it does not have both lower-case and upper-case characters
    #elif (args.case == "ignore") # then do nothing
    #    pass
    
    # Molecule  5'-sequence-3'                                          Description     ignore  upper-only  lower-only  mixed-lower  mixed-upper  mixed-only
    # ========  ======================================================  ==============  ======  ==========  ==========  ===========  ===========  ==========
    # genome    ACCATAGGAATCCAGCGGCGATCTTAAaggaggatctaggtcgatagcggaata  -               -       -           -           -            -            -
    # spacer           GAATCCAGCGGCGATCTTAA                             All upper-case  keep    keep        discard     discard      keep         discard
    # spacer                     GCGATCTTAAaggaggatct                   Mixed case      keep    discard     discard     keep         keep         keep
    # spacer                               aggaggatctaggtcgatag         All lower-case  keep    discard     keep        keep         discard      discard
    
    # Discard potential gRNAs that have mismatches with their target site
    #if (args.case == "invariant-lower"):
    #    pass
    #elif (args.case == "invariant-upper"):
    #    pass
    
    # Convert input sequence to upper-case so it can be evaluated by the scoring algorithms
    seq = seq.upper()
    
    
    if (args.ambiguities == 'discard'):
        # If target sequence has any ambiguities, then discard it
        if regex.search('[^ATCGatcg]', seq):
            return []
        # Otherwise, proceed
        seqs1 = [seq]
    elif (args.ambiguities == 'disambiguate'):
        # Disambiguate sequences if necessary
        seqs1 = nucleotides.disambiguate_iupac(seq)
    elif (args.ambiguities == 'exclusive'):
        # If no ambiguous characters are found, hten discard it
        if not regex.search('[^ATCGatcg]', seq):
            return [] # Reject this sequence
        # Otherwise, disambiguate
        seqs1 = nucleotides.disambiguate_iupac(seq)
    else:
        # Do nothing if just 'keep'
        seqs1 = [seq]
    
    
    # After disambiguation, some spacers may no longer confine to their motifs.
    # Thus, we must remove targets that do not confine to at least one of
    # the defined SPACER>PAM motifs.
    #
    # For instance, if motif is AN{4}>NGG
    #     genome ACCTACATCWAGCTAGGCTCTAA
    #     spacer          WAGCT
    #        pam               AGG
    # expansion1          TAGCT           <-- Discard
    # expansion2          AAGCT           <-- Keep
    # 
    # Also, separate the SPACER and PAM motifs
    seqs2 = []
    targets2 = []
    pams2 = []
    for nt in seqs1:
        # for i in range(len(args.parsed_motifs)):
        for motif in OnTargetMotif.motifs:
        #for spacers, pams, side in args.parsed_motifs:
            #m = nucleotides.motif_conformation(nt, spacers, pams, side)
            #spacers, pams, side = args.parsed_motifs[i]
            spacers, pams, side = motif.parsed_list
            #compiled_regex = args.compiled_motifs[i]
            m = nucleotides.motif_conformation2(nt, side, motif.compiled_regex)
            if m:
                seqs2.append(nt)
                targets2.append(m[0])
                pams2.append(m[1])
                #break # consider breaking here
    #seqs = temp_seqs
    
    # Remove targets with T{5,} (old code)
    #seqs = [ nt for nt in seqs if ('T'*(args.max_consecutive_ts+1) not in nt) ]
    
    # Apply all prefilters to each sequence
    seqs3 = []
    for i in range(len(seqs2)):
        prefilter_passes = []
        prefilter_choices = [C for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.prefilter] # Really, should only be SingleSequenceAlgorithms
        this = (seqs2[i], targets2[i], pams2[i], upstream, downstream) # For GC and PolyT, only 'target' is used
        for C in prefilter_choices:
            c_score = C. calculate(this)
            if (C.minimum <= c_score <= C.maximum):
                prefilter_passes.append(True)
            else:
                prefilter_passes.append(False)
        # Code to handle PairedSequenceAlgorithms
        # parent = (sequence, target, pam, upstream, downstream)
        #for C in algorithms.paired_algorithms:
        #    c_score = C.calculate(parent, this)
        #    if C.prefilter:
        #        if (C.minimum <= c_score <= C.maximum):
        #            prefilter_passes.append(True)
        #        else:
        #            prefilter_passes.append(False)
        
        # If alignment meets all prefilter criteria, then set it as True
        if all(prefilter_passes): # FYI: all([]) returns True
            #seqs3.append(seqs2[i])
            seqs3.append((seqs2[i], targets2[i], pams2[i]))
            
    
    
    
    # # Remove targets whose %GC is outside the chosen bounds
    # temp_seqs3 = []
    # temp_targets3 = []
    # temp_pams3 = []
    # for i in range(len(temp_seqs2)):
    #     if (args.target_gc[0] <= scores.gc_score(temp_targets2[i]) <= args.target_gc[1]):
    #         temp_seqs3.append(temp_seqs2[i])
    #         temp_targets3.append(temp_targets2[i])
    #         temp_pams3.append(temp_pams2[i])
    # #seqs = temp_seqs2
    
    # [(seq, spacer, pam), (seq, spacer, pam)]
    #rets = []
    #for i in range(len(temp_seqs2)):
    #    rets.append((temp_seqs2[i], temp_targets2[i], temp_pams2[i]))
    #return rets
    return seqs3

def save_object(obj, name, path='/dev/shm/addtag'):
    if isinstance(name, str):
        if (len(name) > 0):
            os.makedirs(path, exist_ok=True)
            with open(os.path.join(path, name), 'wb') as flo:
                pickle.dump(obj, flo, pickle.HIGHEST_PROTOCOL)
        else:
            raise Exception("object name is not a string")

def load_object(name, path='/dev/shm/addtag'):
    with open(os.path.join(path, name), 'rb') as flo:
        data = pickle.load(flo)
    return data

def make_labeled_primer_alignments(label_list, sequence_list, contig_name, primer_pair_label_list, primer_pair_list, shifts=None):
    
    lengths = [len(x) for x in sequence_list]
    output = make_alignment_labels(label_list, lengths)
    output_names = ['' for x in output]
    
    output_names.append(contig_name)
    output.append(''.join(sequence_list))
    
    for i, (pp_label, pp) in enumerate(zip(primer_pair_label_list, primer_pair_list)):
        output_names.append(pp_label)
        if (shifts == None):
            output.append(' '*pp.forward_primer.position + pp.get_formatted())
        else: # Need to improve this so it uses shifts for real
            if pp:
                if (i in [0, 1]):
                    output.append(' '*pp.forward_primer.position + pp.get_formatted())
                elif (i == 2):
                    output.append(' '*(lengths[0]+lengths[1]) + ' '*pp.forward_primer.position + pp.get_formatted()) # NOT elegant AT ALL
            else:
                output.append('None')
    
    n = max(len(x) for x in output_names)
    lines = []
    for i in range(len(output)):
        lines.append(output_names[i].ljust(n, ' ') + ' ' + output[i])
    
    return lines

def make_alignment_labels(labels, lengths, edges_open=False):
    """
    Prints multi-line set of exclusive sequences and their labels
    
    Example:
      print_alignment_labels(['us region', 'label', 'us homology', 'insert',
          'ds homology', 'ds region'], [8, 5, 15, 3, 15, 10])
       ┌us region
       │       ┌label              ┌insert           ┌ds region
      ┌┴─────┐┌┴──┐┌─us homology─┐┌┴┐┌─ds homology─┐┌┴───────┐
    """
    lines = ['']
    for i, (label, slen) in enumerate(zip(labels, lengths)):
        
        # Modify the labels
        if (len(label)+2 > slen):
            if (slen == 0):
                if (len(label) > 0):
                    corner = '┌'
                else:
                    corner = ''
            elif (slen == 1):
                if (len(label) > 0):
                    corner = '┌'
                else:
                    corner = ''
            else:
                corner = ' ┌'
            #for j in range(1, len(lines)):
            for j in range(len(lines)-1, 0, -1):
                if (len(lines[j]) <= len(lines[0])):
                    lines[j] += ' '*(len(lines[0])-len(lines[j])) + corner+label
                    break
            else:
                if (len(lines) > 1):
                    new_line = ''
                    for k in range(0, len(lines[0])):
                        if (lines[-1][k] in ['┌', '│', '¦']):
                            if lines[0][k] in ['│', '┤', '┴']:
                                new_line += '│'
                            else:
                                new_line += '¦'
                        else:
                            new_line += ' '
                    if (lengths[i-1] == 0):
                        if (slen == 0):
                            corner = '┌'
                        else:
                            corner = '¦┌'
                    lines.append(new_line + corner+label)
                else:
                    lines.append(' '*len(lines[0]) + corner+label)
        
        # Modify lines[0]
        if (slen == 0):
            pass
        elif (slen == 1):
            if (len(label) > 0):
                lines[0] += '│'
            else:
                lines[0] += '╥'
        elif (slen == 2):
            if (len(label) > 0):
                lines[0] += '┌┤'
            else:
                lines[0] += '┌┐'
        elif (len(label)+2 > slen):
            lines[0] += '┌┴'+'─'*(slen-3)+'┐'
        else:
            line = '┌'+'─'*(slen-2)+'┐'
            pos = slen//2-len(label)//2
            lines[0] += line[:pos] + label + line[pos+len(label):]
    
    #for line in lines[1:]:
    #    print(line)
    #print(lines[0])
    return lines[1:] + [lines[0]]
