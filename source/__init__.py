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

# Import non-standard packages
import regex

# Import included AddTag-specific modules
from . import utils
from . import nucleotides
from . import scores
from . import algorithms
from . import aligners
from . import oligos

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
  folder/excision-query.sam         SAM file for alignment of candidate spacers
  folder/excision-query.err         STDOUT/STDERR from alignment
  folder/excision-spacers.fasta     Spacer sequences with scores
  folder/excision-constructs.fasta  Construct sequences to be synthesized
  folder/excision-dDNAs.fasta       dDNA sequences to be synthesized for
                                    knock-out
  folder/reversion-query.fasta      FASTA file of candidate knock-in spacers
  folder/reversion-query.sam        SAM file for alignment of candidate spacers
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
    """Class primarily representing a SAM alignment"""
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
            for sequence, obj in cls.sequences.items():
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
    def exhaustive_site_search(cls, args, f, contig_sequence, orientation):
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
                                        # If it is not one of the two failure states, then add it
                                        if not ((t[2] < len(upstream)) or (t[1] >= len(upstream) + len(mAT))):
                                            cls(f.name, f.contig, orientation, dDNA, (start1, end1), mAT, (start2, end2), spacer=t)
    
    @classmethod
    def generate_donors(cls, args, contigs):
        """
        Creates the DNA oligo with the structure:
        [upstream homology][unique gRNA][downstream homology]
        that excises the target feature
        """
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
            #cls.exhaustive_site_search(args, feature, contig, start, end, strand, contig_sequence, orientation)
            cls.exhaustive_site_search(args, f, contig_sequence, orientation)
    
    @classmethod
    def generate_alignments(cls):
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
    def load_sam(cls, filename, args, contigs, sep=':'):
        """
        Read in SAM file.
        sep is the separator for the header. Positions are converted to 0-index
        Creates a list of Sequence objects
        """
        
        # Sequence Alignment/Map (SAM) format is TAB-delimited. Apart from the
        # header lines, which are started with the '@' symbol, each alignment line
        # consists of:
        #   Col  Field  Description
        #   1    QNAME  Query template/pair NAME
        #   2    FLAG   bitwise FLAG
        #   3    RNAME  Reference sequence NAME
        #   4    POS    1-based leftmost POSition/coordinate of clipped sequence
        #   5    MAPQ   MAPping Quality (Phred-scaled)
        #   6    CIAGR  extended CIGAR string
        #   7    MRNM   Mate Reference sequence NaMe ('=' if same as RNAME)
        #   8    MPOS   1-based Mate POSistion
        #   9    TLEN   inferred Template LENgth (insert size)
        #   10   SEQ    query SEQuence on the same strand as the reference
        #   11   QUAL   query QUALity (ASCII-33 gives the Phred base quality)
        #   12+  OPT    variable OPTional fields in the format TAG:VTYPE:VALUE
        
        # Code to decompress a *.bam file should go here
        
        with open(filename, 'r') as flo:
            for line in flo:
                if not line.startswith('@'):
                    sline = line.rstrip().split("\t")
                    if ((len(sline) > 5) and (sline[2] != '*')):
                        target = cls.indices[sline[0]] # key=exTarget-519
                        
                        # Get orientation
                        alignment_orientation = utils.sam_orientation(int(sline[1]))
                        
                        # Get alignment position
                        alignment_contig = sline[2]
                        alignment_start = int(sline[3])-1
                        alignment_end = int(sline[3])-1+utils.cigar_length(sline[5])
                        
                        # Reverse-complement if needed
                        alignment_contig_sequence = contigs[alignment_contig]
                        alignment_sequence = alignment_contig_sequence[alignment_start:alignment_end]
                        alignment_upstream = alignment_contig_sequence[alignment_start-10:alignment_start]
                        alignment_downstream = alignment_contig_sequence[alignment_end:alignment_end+10]
                        actual_sequence = sline[9]
                        if (alignment_orientation == '-'):
                            alignment_sequence = nucleotides.rc(alignment_sequence)
                            alignment_upstream, alignment_downstream = nucleotides.rc(alignment_downstream), nucleotides.rc(alignment_upstream)
                            actual_sequence = nucleotides.rc(actual_sequence)
                        
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
        
        logging.info(cls.__name__ + ' SAM file parsed: {!r}'.format(filename))
    
    @classmethod
    def generate_query_fasta(cls, filename, sep=':'):
        """
        Creates a FASTA file (That becomes header for SAM file)
        >exTarget-id feature:contig:orientation:start..end ...
        """
        with open(filename, 'w') as flo:
            for sequence, obj in cls.sequences.items():
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
                ]), file=flo)
                print(sequence, file=flo)
        logging.info(cls.__name__ + ' spacers FASTA generated: {!r}'.format(filename))
        return filename
    
    def format_location(self, location, sep=':'):
        feature, contig, orientation, start, end, upstream, downstream = location
        return sep.join([feature, contig, orientation, '..'.join([str(start), str(end)])])
    
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
        targets = set()
        #for seq_i, sequence in enumerate(dDNAs):
        for orientation in ['+', '-']:
            if (orientation == '-'):
                sequence = nucleotides.rc(sequence)
            
            #for i in range(len(args.parsed_motifs)):
            for mymotif in OnTargetMotif.motifs:
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
        r_spacer = None
        r_pam = None
        r_motif = None
        #for i in range(len(args.parsed_motifs)):
        for motif in OnTargetMotif.motifs:
            #spacers, pams, side = args.parsed_motifs[i]
            spacers, pams, side = motif.parsed_list
            #compiled_regex = args.compiled_motifs[i]
            m = nucleotides.motif_conformation2(sequence, side, motif.compiled_regex)
            
            # If the motif matches, then return it
            if m:
                r_spacer = m[0]
                r_pam = m[1]
                #r_motif = args.motifs[i]
                r_motif = motif
                break
            # If the motif does not match, then fudge it for this motif
            # This will return the last-fudged motif, and not necessarily the best-fudged motif
            # Will need to improve this code later
            else:
                if (side == '>'):
                    l = max(map(len, pams))
                    r_spacer = sequence[:-l]
                    r_pam = sequence[-l:]
                    #r_motif = args.motifs[i]
                    r_motif = motif
                elif (side == '<'):
                    l = max(map(len, pams))
                    r_spacer = seq[:l]
                    r_pam = seq[l:]
                    #r_motif = args.motifs[i]
                    r_motif = motif
        return r_spacer, r_pam, r_motif
    
    def get_features(self):
        """Return (sorted) list of all feature names this Target maps to"""
        return sorted(set(x[0] for x in self.locations))
    
    def get_contigs(self):
        return sorted(set(x[1] for x in self.locations))
    
    def get_location_features(self):
        """Return (sorted) list of all feature names this Target maps to"""
        return sorted(set(x[0] for x in self.locations))
    
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
            self.spacer + '|' + self.pam,
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
            ef = f.expand_feature(args, contig_sequence)
            
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
        #  --motifs, where '.' is any non-spacer nucleotide
        #    'CCN<N{17,20}.{13,18}N{17,20}>NGG'     # PAM-out
        #    'N{17,20}>NGG.{13,18}CCN<N{17,20}'     # PAM-in
        
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

class Feature(object):
    features = {}
    # Origin flags
    NONE=0
    INPUT=1
    DERIVED=2
    def __init__(self, contig, start, end, strand, name=None, attributes=None, source=None, feature_type=None, score=None, frame=None, origin=0, sep=';'):
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
                start = int(sline[3])-1
                end = int(sline[4])
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
    def load_gff_file(cls, filename, feature_types, selected_features, tag):
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
                    if obj.feature_type in feature_types:
                        obj.name = obj.attributes[tag]
                        
                        # Reduce the total number of features to just the ones indicated in the selection
                        if selected_features:
                            if (obj.name in selected_features):
                                cls.features[obj.name] = obj
                        else:
                            cls.features[obj.name] = obj
        
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
    
    def expand_feature(self, args, contig_sequence, expansion_size=100, minimum_targets_per_feature=10):
        feature_start = self.start
        feature_end = self.end
        
        contig_length = len(contig_sequence)
        
        if (feature_end == None):
            feature_end = len(contig_sequence)
        
        # Gradually expand feature size until the minimum number of targets is found
        feature_sequence = contig_sequence[feature_start:feature_end]
        targets = Target.get_targets(args, feature_sequence) # Does both orientations (+/-)
        count = 0
        while ((len(targets) < minimum_targets_per_feature) and ((feature_start, feature_end) != (0, contig_length))):
            feature_start = max(0, feature_start-expansion_size)
            feature_end = min(contig_length, feature_end + expansion_size)
            
            feature_sequence = contig_sequence[feature_start:feature_end]
            targets = Target.get_targets(args, feature_sequence) # Does both orientations (+/-)
            count += 1
        
        # Save the new feature as a Feature object
        if (count > 0):
            new_name = f.name + '_derived'
            new_attributes = f.attributes.copy()
            new_attributes[args.tag] = new_name
            new_feature = Feature(f.contig, feature_start, feature_end, f.strand, source=f.source, feature_type=f.feature_type, score=f.score, frame=f.frame, attributes=new_attributes, origin=Feature.DERIVED)
            new_feature.name = new_name
            Feature.features[new_feature.name] = new_feature
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
        nmin = 1
        nmax = 2
        if not (nmin <= len(values) <= nmax):
            raise argparse.ArgumentTypeError('argument --{f}: expected between {nmin} and {nmax} arguments'.format(f=self.dest, nmin=nmin, nmax=nmax))
        
        # If only 1 argument, then it is a filename
        if (len(values) == 1):
            pass
        elif (len(values) == 2):
            # If 2 arguments, then the first is either "uniform" or "specific"
            s0 = set(values[0])
            valid_s0_values = {'uniform', 'specific'} # set
            d0 = s0.difference(valid_s0_values)
            if (len(d0) > 0):
                raise ValueError('invalid design TYPE: %s (choose from {uniform,specific})' % d.pop())
            
            s1 = set(values[1]) # in case there are duplicates, reduce them to a single occurrence with set()
            valid_s1_values = {'single', 'paired'} # set
            d1 = s1.difference(valid_s1_values)
            if (len(d1) > 0):
                raise ValueError('invalid design TYPE: %s (choose from {single,paired})' % d.pop())
        
        setattr(args, self.dest, list(s))

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
        
        # Get timestamp for analysis beginning
        start_time = time.time()
        
        if hasattr(args, 'folder'):
            # Create the project directory if it doesn't exist
            os.makedirs(args.folder, exist_ok=True)
        
            # Create the logger
            logging.basicConfig(filename=os.path.join(args.folder, 'log.txt'), level=logging.INFO, format='%(message)s') # format='%(levelname)s %(asctime)s: %(message)s'
        
        # Echo the command line parameters to STDOUT and the log
        print(args, flush=True)
        
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
        # End
    
    def _feature(self, args):
        print("Search GFF here.")
    
    def _evaluate(self, args):
        '''UNDER DEVELOPMENT'''
        print("Perform the evaluation here.")
        
        contig_sequences = utils.load_multiple_fasta_files(args.fasta)
        
        dDNA_file = ExcisionDonor.generate_fasta(os.path.join(args.folder, 'excision-dDNAs.fasta'))
        genome_fasta_file = utils.write_merged_fasta(contig_sequences, os.path.join(args.folder, 'genome.fasta'))
        
        # search spacer FASTA against genome+dDNA FASTA
        genome_index_file = args.selected_aligner.index(genome_fasta_file, os.path.basename(genome_fasta_file), args.folder, args.processors)
        dDNA_index_file = args.selected_aligner.index(dDNA_file, os.path.basename(dDNA_file), args.folder, args.processors)
        
        q2gDNA_sam_file = args.selected_aligner.align(ex_query_file, genome_index_file, 'excision-query-2-gDNA.sam', args.folder, args.processors)
        q2exdDNA_sam_file = args.selected_aligner.align(ex_query_file, ex_dDNA_index_file, 'excision-query-2-excision-dDNA.sam', args.folder, args.processors)
        
        ExcisionTarget.load_sam(q2gDNA_sam_file, args, contig_sequences)
        ExcisionTarget.load_sam(q2exdDNA_sam_file, args, ExcisionDonor.get_contig_dict())
        
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
        
        
        
    
    def _generate(self, args):
        """Perform complete CRISPR/Cas analysis for input"""
        
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
        
        # Parse the on-target motifs
        for motif in args.motifs:
            OnTargetMotif(motif)
        
        # Parse the off-target motifs
        for motif in set(args.off_target_motifs).difference(args.motifs): # These motifs are exclusive to off-target list
            OffTargetMotif(motif)
        
        # Motif.motifs           list of ALL motifs, both on-target and off-target?  <-- not implemented yet
        # OnTargetMotif.motifs   list of on-target motifs
        # OffTargetMotif.motifs  list of off-target motifs
        
        # Add 'args.selected_aligner' to hold the actual aligner object
        for a in aligners.aligners:
            if (a.name == args.aligner):
                args.selected_aligner = a
                break
        
        #contigs = utils.load_fasta_file(args.fasta[0]) # To do --> Do this for all FASTA files, then merge them into the same dictionary
        
        # Load the FASTA file specified on the command line
        # Merge all sequence information into the same dictionary
        #fasta_index, contig_index, contig_sequences = utils.load_indexed_fasta_files(args.fasta)
        contig_sequences = utils.load_multiple_fasta_files(args.fasta)
        
        # Open and parse the GFF file specified on the command line
        #features = utils.load_gff_file(args.gff, args.features, args.tag)
        # Filter features by what is selected
        #features = self.filter_features(features, args.selection)
        Feature.load_gff_file(args.gff, args.features, args.selection, args.tag)
        Feature.assert_features(args.selection, contig_sequences)
        
        # Make index of homologs
        if args.homologs:
            homologs, feature2gene = utils.load_homologs(args.homologs)
        else:
            homologs, feature2gene = None, None
        
        # Merge features?
        #features = merge_features(features)
        
        
        #### Some checks that need to be added ####
        # The feature being targeted MUST not go up to the edge of the contigs
        # Otherwise there will be no junction. i.e. the upstream and downstream
        # regions won't exist.
        # Actually, these regions must be at minimum, 50 nt
        ###########################################
        
        
        if (args.ko_gRNA): # (True/False)
            # Search for good targets within specified features
            pass
            
        if (args.ko_dDNA): # (mintag/addtag/unitag/bartag)
            if (args.ko_dDNA == 'mintag'):
                pass
            elif (args.ko_dDNA == 'addtag'):
                pass
            elif (args.ko_dDNA == 'unitag'):
                pass
            elif (args.ko_dDNA == 'bartag'):
                pass
        
        if (args.ki_gRNA): # (True/False)
            pass
        
        if (args.ki_dDNA != None): # (None/True/'*.fasta')
            #if (len(args.ki_dDNA) == 1):
            if isinstance(args.ki_dDNA, str):
                # If a file is specified, then it has the knock-in DNA
                # that should be stitched to flanking homology arms
                pass
            else:
                # Otherwise, generate KI dDNA that are wild type
                pass
        
        
        # Search features within contigs for targets that match the motifs
        # Old code (without feature expansion) ExcisionTarget.get_targets(args, contig_sequences, features)
        ExcisionTarget.search_all_features(args, contig_sequences)
        
        # Write the query list to FASTA
        ex_query_file = ExcisionTarget.generate_query_fasta(os.path.join(args.folder, 'excision-query.fasta'))
        
        # Generate excision dDNAs and their associated reversion gRNA spacers
        ExcisionDonor.generate_donors(args, contig_sequences)
        ReversionTarget.get_targets()
        ex_dDNA_file = ExcisionDonor.generate_fasta(os.path.join(args.folder, 'excision-dDNAs.fasta')) # Program will fail with error if this file is empty...
        re_query_file = ReversionTarget.generate_query_fasta(os.path.join(args.folder, 'reversion-query.fasta'))
        
        # Generate reversion dDNAs and write them to FASTA
        ReversionDonor.generate_donors(args, contig_sequences)
        re_dDNA_file = ReversionDonor.generate_fasta(os.path.join(args.folder, 'reversion-dDNAs.fasta'))
        
        # Merge input FASTA files into a single one
        genome_fasta_file = utils.write_merged_fasta(contig_sequences, os.path.join(args.folder, 'genome.fasta'))
        
        # Index args.fasta for alignment
        #index_file = index_reference(args)
        genome_index_file = args.selected_aligner.index(genome_fasta_file, os.path.basename(genome_fasta_file), args.folder, args.processors)
        ex_dDNA_index_file = args.selected_aligner.index(ex_dDNA_file, os.path.basename(ex_dDNA_file), args.folder, args.processors)
        re_dDNA_index_file = args.selected_aligner.index(re_dDNA_file, os.path.basename(re_dDNA_file), args.folder, args.processors)
        
        # Use selected alignment program to find all matches in the genome and dDNAs
        #ex_genome_sam_file = align(ex_query_file, genome_index_file, args)
        exq2gDNA_sam_file = args.selected_aligner.align(ex_query_file, genome_index_file, 'excision-query-2-gDNA.sam', args.folder, args.processors)
        exq2exdDNA_sam_file = args.selected_aligner.align(ex_query_file, ex_dDNA_index_file, 'excision-query-2-excision-dDNA.sam', args.folder, args.processors)
        
        #print("ExcisionTarget before SAM parsing")
        #for et_seq, et_obj in ExcisionTarget.sequences.items():
        #    print(et_obj)
        
        # Load the SAM files and add Alignments to ExcisionTarget sequences
        ExcisionTarget.load_sam(exq2gDNA_sam_file, args, contig_sequences)
        ExcisionTarget.load_sam(exq2exdDNA_sam_file, args, ExcisionDonor.get_contig_dict())
        
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
        #re_sam_file = align(re_query_file, genome_index_file, args)
        req2gDNA_sam_file = args.selected_aligner.align(re_query_file, genome_index_file, 'reversion-query-2-gDNA.sam', args.folder, args.processors)
        req2exdDNA_sam_file = args.selected_aligner.align(re_query_file, ex_dDNA_index_file, 'reversion-query-2-excision-dDNA.sam', args.folder, args.processors)
        req2redDNA_sam_file = args.selected_aligner.align(re_query_file, re_dDNA_index_file, 'reversion-query-2-reversion-dDNA.sam', args.folder, args.processors)
        
        # Load the SAM files and add Alignments to ReversionTarget sequences
        ReversionTarget.load_sam(req2gDNA_sam_file, args, contig_sequences)
        ReversionTarget.load_sam(req2exdDNA_sam_file, args, ExcisionDonor.get_contig_dict())
        ReversionTarget.load_sam(req2redDNA_sam_file, args, ReversionDonor.get_contig_dict())
        
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
        self.get_best(args, homologs)
        self.get_best_table(args, homologs, feature2gene)
        
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
    
    def _confirm(self, args):
        print("Design conformation primers here.")
        
        # Add 'args.selected_oligo' to hold the actual aligner object
        for o in oligos.oligos:
            if (o.name == args.oligo):
                args.selected_oligo = o
                break
        
        genome_contigs = utils.load_multiple_fasta_files(args.fasta)
        
        rounds = utils.load_fasta_files_into_list(args.dDNAs)
        
        #args.number_pcr_conditions
        
        
        # use local alignment to find where left- and right- homology arms of each dDNA are in the genome
        # genome     ─genome┐┌─us_homology─┐┌─feature─┐┌─ds_homology─┐┌genome─
        # dDNA               ┌─us_homology─┐┌──*tag*──┐┌─ds_homology─┐
        
        # dDNA    ATCA.CG.ATAC
        # dDNA   (ATCA).*(ATAC)
        # genome AGTCAGCATCAGACGCGACTCAGCGCGGAGCATCTATCAGCCGCGAAATGATATACGCGCTCTGTGTGAATTAACACATATAGAGAAAAGCGCCTGATTATATATCTCTCGGTGTGCGCGATGGGGACTAG
        #               ATCA-----------------------------------------ATAC
        #                                           ATCA-------------ATAC
        
        results = []
        for r, dDNA_contigs in enumerate(rounds):
            print(r)
            for dDNA_name, dDNA_seq in dDNA_contigs.items():
                print(dDNA_name, flush=True)
                C = nucleotides.build_dDNA_regex(dDNA_seq)
                print('compiled', flush=True)
                for contig_name, contig_seq in genome_contigs.items():
                    print(contig_name, flush=True)
                    
                    # Finds most (but not all) of the possible places the dDNA could integrate
                    for m in C.finditer(contig_seq, overlapped=True):
                        print(m)
                        
                        
                        # Measuring wild-type alleles
                        feature_sequence = m.group(1) # The capturing group holds the non-flanking DNA (i.e. the feature that is disrupted)
                        upstream_sequence = contig_seq[m.start()-500:m.start()]
                        downstream_sequence = contig_seq[m.end():m.end()+500]
                        
                        
                        #primer_pairs = args.selected_oligo.scan_sequence(seq) # <---- replace with correct stuff
                        #----> Need a sanity check to test for how unique each generated primer is across the genome <-----
                        fseq = None
                        rseq = None
                        fpos = None
                        rpos = None
                        size = None
                        tm = None
                        
                        results.append([r, args.dDNAs[r], dDNA_name, contig_name, fseq, fpos, rseq, rpos, size, tm])
        
        # Output
        # Round   Filename       Contig     Mapping   Fseq           Fpos    Rseq         Rpos    size    tm
        # 0       genome.fasta   -          chr1      ACAAAGGCTAGG   12000   GTGATCGAAG   14000   2000    61.2
        # 1       dDNA1.fasta    feature1   chr1      GCAAATCAAAG    11000   CGCTGCATACC  13000   2000    60.8
        # 1       dDNA1.fasta    feature2   chr1      CCATAGCCCAGC   12345   CACAGGTGC    12300   123     58.6
        headers = ['Round', 'Filename', 'Contig', 'Mapping', 'Fseq', 'Fpos', 'Rseq', 'Rpos', 'size', 'tm']
        print('\t'.join(headers))
        for sline in results:
            print('\t'.join(map(str, sline)))
    
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
        __glossary_description__ = "Show glossary of common CRISPR/Cas terms, then exit."
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
        __motifs_description__ = "Show list of common CRISPR/Cas SPACER>PAM arrangements, then exit."
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
        __algorithms_description__ = "Show list of all implemented gRNA evaluation algorithms."
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
        __aligners_description__ = "Show list of all supported alignment programs."
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
        __oligos_description__ = "Show list of all supported oligonucleotide thermodynamics property programs."
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
        
        __search_description__ = "Search input FASTA for DNA sequence, and output a GFF file for use as input to 'generate'"
        __search_help__ = "Search for positions of DNA in FASTA to make GFF file."
        parser_search = subparsers.add_parser('search',
            description=__search_description__,
            #epilog=__feature_epilog__,
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
        
        required_group.add_argument("--query", required=True, nargs="+", metavar="SEQUENCE", type=str,
            help="Sequence to search")
        
        # Add optional arguments
        parser_search.add_argument("--identifier", metavar="FEATURE", type=str, nargs="+",
            help="Identifier for input query sequences")
        
        parser_search.add_argument("--tag", metavar="TAG", type=str, default="ID",
            help="GFF3 attribute tag")
        
        return parser_search
    
    def _parser_feature(self, subparsers):
        ''' "feature" parser '''
        __feature_description__ = "Search features for specific gene name or id."
        __feature_help__ = "Search features for specific gene name or id."
        parser_feature = subparsers.add_parser('feature',
            description=__feature_description__,
            #epilog=__feature_epilog__,
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
        # Add optional arguments
        return parser_feature
    
    def _parser_evaluate(self, subparsers):
        ''' "evaluate" parser '''
        __evaluate_description__ = 'Evaluate pre-designed CRISPR/Cas oligonucleotide sequences.'
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
        __generate_description__ = "Design full sets of oligonucleotide sequences for CRISPR/Cas genome engineering experiment."
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
        required_group.add_argument("--gff", required=True, metavar="*.gff", type=str,
            help="GFF file specifying chromosomal features that should be \
            multiplexed together.")
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
            default=["N{20}>NGG"],
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
            Examples: 'CDS', 'gene', 'mRNA', 'exon', 'intron', 'tRNA', 'rRNA'")
        parser_generate.add_argument("--warning_features", metavar='FEATURE', nargs="+", type=str, default=['all'],
            help="GFF tags that will trigger a warning if they overlap with the \
            target feature. Examples: 'CDS', 'gene', 'mRNA', 'exon', 'intron', 'tRNA', 'rRNA'")
        parser_generate.add_argument("--dDNA_gDNA_ratio", metavar="N", type=int, default=1000,
            help="Ratio of donor DNA to genomic DNA for calculating off-target scores")
        #parser_generate.add_argument("--target_gc", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[25, 75],
        #    help="Generated gRNA spacers must have %%GC content between these values (excludes PAM motif)") # Moved to prefilter
        #
        #
        parser_generate.add_argument("--excise_upstream_homology", nargs=2, metavar=("MIN", "MAX"), type=int, default=[50,50],
            help="Range of homology lengths acceptable for knock-out dDNAs, inclusive.")
        parser_generate.add_argument("--excise_downstream_homology", nargs=2, metavar=("MIN", "MAX"), type=int, default=[47,50],
            help="Range of homology lengths acceptable for knock-out dDNAs, inclusive.")
        #
        #
        parser_generate.add_argument("--excise_donor_lengths", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[100, 100],
            help="Range of lengths acceptable for knock-out dDNAs, inclusive.")
        parser_generate.add_argument("--excise_insert_lengths", nargs=2, metavar=("MIN", "MAX"), type=int, default=[0,3],
            help="Range for inserted DNA lengths, inclusive (mini-AddTag, mAT). \
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
        parser_generate.add_argument("--revert_upstream_homology", nargs=2, metavar=("MIN", "MAX"), type=int, default=[300,300],
            help="Range of homology lengths acceptable for knock-in dDNAs, inclusive.")
        parser_generate.add_argument("--revert_downstream_homology", nargs=2, metavar=("MIN", "MAX"), type=int, default=[300,300],
            help="Range of homology lengths acceptable for knock-in dDNAs, inclusive.")
        parser_generate.add_argument("--revert_donor_lengths", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[0, 100000],
            help="Range of lengths acceptable for knock-in dDNAs.")
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
        
        
        # Temporary stopgap to make sure the calculations don't take too long
        parser_generate.add_argument("--max_time", metavar="SECONDS", type=float, default=60,
            help="Maximum amount of time, in seconds, for each feature to spend calculating dDNAs.")
        
        
        
        parser_generate.add_argument("--bartag_motif", type=str, default='N{20}',
            help="Structure of nucleotides that should be generated")
        
        parser_generate.add_argument("--bartag_distance", type=int, default=3,
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
        
        parser_generate.add_argument("--ko-dDNA", type=str, default=None,
            choices=['mintag', 'addtag', 'unitag', 'bartag'],
            help="'mintag' are unique us/i/ds junction targets specific to each feature. \
            'addtag' are unique targets for each feature. \
            'unitag' is a single and uniform, unique target for ALL features. \
            'bartag' are unique barcodes for each feature (does not guarantee targets).")
        
        parser_generate.add_argument("--ki-gRNA", action='store_true', default=False,
            help="Design gRNAs to target either the dDNA. \
            Defaults to True if '--ko-dDNA mintag' is specified.")
        
        parser_generate.add_argument("--ki-dDNA", nargs='?', type=str, default=None,
            metavar='*.fasta', const=True,
            help="If only the option with no file is specified, then design wild type dDNA. \
            If FASTA file is specified, then design dDNA to replace features with foreign sequences. \
            The sequence header should have a 'TAG=*' field where 'TAG' is the tag specified with '--tag' option.")
        
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
        __confirm_description__ = "Design primers for confirming whether each step of genome engineering is successful or not. This does not design multiplexable primers."
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
        
        # Add optional arguments
        parser_confirm.add_argument("--number_pcr_conditions", type=int, default=None,
            help="Number of PCR conditions to develop primers for. All amplicons \
            within each condition will have similar size (nt) and melting temperatures. \
            If unspecified, will default to the number of target features.")
        
        oligo_choices = [x.name for x in oligos.oligos]
        parser_confirm.add_argument("--oligo", type=str, choices=oligo_choices, default='Primer3',
            help="Program to perform thermodynamic calculations.")
        
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
            #for f in features:
                groups.add(tuple(sorted(homologs[feature_name])))
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
            groups = set()
            for feature_name, f in Feature.features.items():
            #for f in features:
                #contig, start, end, strand = features[feature]
                # Convert to tuple
                groups.add(tuple(sorted(homologs[feature_name])))
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
                
                sline = [gene, feature, key0s, translations, weight, othz, otcfd, azimuth, obj.name, obj.spacer+'|'+obj.pam, exdonors]
                
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
                
                sline = [gene, csfeatures, weight, othz, otcfd, azimuth, obj.name, obj.spacer+'|'+obj.pam, redonors]
                
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
        
    
    def get_best(self, args, homologs):
        """Function that returns the best spacers and dDNAs for each feature"""
        
        display_num = 5
        
        # Print best ReversionTargets calculated and their corresponding ExcisionDonors
        ret_dict2 = self.get_reTarget_homologs(homologs)
        for k in ret_dict2:
            logging.info(str(k) + ' ' + str(len(ret_dict2[k])))
            # Print the top 5
            for rank, obj in self.rank_targets(ret_dict2[k])[:display_num]:
                logging.info(' ' + str(rank) + ' ' + str(obj))
                # Get the ExcisionDonor objects for this ki-spacer, and rank them
                rds = self.rank_donors(obj.get_donors())
                # filter out all but the top-ranked ones
                rds = [x for x in rds if (x[0] == rds[0][0])]
                for gap, exd_obj in rds:
                    logging.info('   ' + str(gap) + ' ' + str(exd_obj.get_trims()) + ' ' + str(exd_obj))
        
        # Print best ExcisionTargets (not necessarily homozygous) for each feature
        # and the ReversionDonor
        for feature in sorted(Feature.features):
            logging.info(feature)
            et_list = []
            for name, obj in ExcisionTarget.indices.items():
                if feature in obj.get_location_features():
                    et_list.append(obj)
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
