#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/__init__.py

# Import standard packages
import sys
import os
import argparse
import textwrap
import time

# Import non-standard packages
import regex

# Import included AddTag-specific modules
from . import utils
from . import nucleotides
from . import scores
from . import algorithms

from .aligners import bowtie2

# Define meta variables
__author__ = "Thaddeus D. Seher (@tdseher) & Aaron Hernday"
__date__ = utils.load_git_date()
__fullversion__ = utils.load_git_version()
__version__ = __fullversion__[:7]
__commits__ = utils.load_git_commits()
__program__ = os.path.basename(sys.argv[0])
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
                            cut┐    3'╤╤╤╤╗            ╔╤╤╗   ││
           ┌──╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥ ╥╥╥┐      ╚╦╦╦╦╦╦╦╦╦╦╦╦╝  ╢   ││
           │5'╩╩╩╩╩╩╩╩╩╩╩╩╩╩╩╩╩═╩╩╩╪╧╧╧╧╧╧╧╩╩╩╩╩╩╩╩╩╩╩╩╧╧╧╝   ││
           │  └────────────────────│──────crRNA───────┘└──────┘│
           │  └──────spacer───────┘│└──────────scaffold────────┘
    3'╥╥╥╥╥┘                cut┐   └╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥╥5'
    5'╨╨╨╨╨───┴┴┴┴┴┴┴┴┴┴┴┴┴┴┴┴┴ ┴┴┴─╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨╨3' genome
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
  
  Copyright (c) 2017 {__author__}.
  All rights reserved.

version:
  short   {__version__}
  full    {__fullversion__}
  commits {__commits__}
  date    {__date__}

protein:
  The Cas9 or Cpf1 protein you use should be engineered specifically for your
  organism. It should be codon-optomized, and if using eukarya, contain an
  appropriate nuclear localization sequence. Cas9 from different species bind
  to different PAM sequences, which is useful when no suitable PAM sequence is
  present within your gene of interest. Additionally, the different Cas9 gene
  sequences can have huge length differences. Remember that each Cas9 is only
  compatible with the tracrRNA and crRNA (or synthetic gRNA) derived from the
  same species.

glossary:
  dDNA        Donor DNA. The DNA that is knocked-in using endogenous cellular
              homologous recombination machinery.
  gDNA        Genomic DNA.
  HDR         Homology-Directed Repair, a DNA repair mechanism that uses a
              template to repair nicks or double-stranded breaks.
  PAM         The Protospacer Adjacent Motif is a short nucleotide sequence that
              serves as a binding signal for Cas9, and its presence is a strict
              requirement for Cas9-mediated DNA cleavage. Must be at 3' end for
              Cas9, and 5' end for Cpf1. The PAM sequence is on the genomic DNA,
              not the gRNA.
  spacer      The ~20 nt element that is homologous to the region of your
              feature of interest, and will direct Cas9 nuclease activity.
              The portion of the crRNA (or sgRNA) that is complementary to the
              genomic DNA target sequence ~20 nt.
  target      The ~20 nt genomic sequence that precedes the PAM sequence.
  protospacer Short genomic DNA sequences ~20 nt of foreign DNA separated by a
              short palindromic repeat and kept like a record against future
              encounters.
  pre-crRNA   CRISPR array of protospacers is transcribed into pre-crRNA.
  crRNA       pre-crRNA is processed (cut up) to produce a set of crRNAs.
              CRISPR-targeting RNA that contains both the ~20 base spacer
              element and additional nucleotides which are complementary to
              the tracrRNA. crRNA is variable.
  tracrRNA    Hybridizes to the crRNA and binds to the CAS9 protein activating
              the complex to creating double-stranded breaks at specific sites
              within genomic sequence.
              Trans-activating crRNA (which serves as the Cas9 nuclease-
              recruiting sequence?) that has sequence complementary to the
              palindromic repeat. When the tracrRNA hybridizes to the short
              palindromic repeat, it triggers processing by the bacterial
              double-stranded RNA-specific ribonuclease, RNase III. Any crRNA
              and the tracrRNA can then both bind to the Cas9 nuclease, which
              then becomes activated and specific to the DNA sequence
              complimentary to the crRNA. tracrRNA is invariable, and is
              specific to each Cas9 protein.
  gRNA        guide RNA, spacer + scaffold
              A synthetic fusion of the endogenous bacterial crRNA and tracrRNA
              sequences. Provides both targeting specificity and
              scaffolding/binding ability for Cas9 nuclease. Does not exist in
              nature. Also referred to as sgRNA.
  sgRNA       Synthetic guide RNA, or single guide RNA (synonymous with 'gRNA').
              Combines the tracrRNA and crRNA, which are separate molecules in
              the native CRISPR/Cas9 system in S. pyogenes, into a single RNA
              construct, simplifying the components needed to use CRISPR/Cas9
              for genome editing (for plasmid or IVT expression).
              A linker loop sequence is included between the two.
  scaffold    The sequence within the gRNA that is responsible for Cas9 binding.
              Does not include the 20 nt spacer/targeting sequence that is used
              to guide Cas9 to target DNA.
  Cas9        Cas9 family nucleases
  dCas9       Catalytically "dead" Cas9 protein
  FokI-dCas9  dCas9 fused with the dimerization-dependent FokI nuclease domain:
              creates a dimeric RNA-guided FokI-dCas9 nuclease (RFN)
              architecture requiring recognition of extended double-length
              target sites for efficient cleavage. Amino-terminal fusions of
              FokI to dCas9 can recognize two 20-nucleotide 'half-sites' in a
              'PAM-out' orientation separated by a 13-18 bp spacer and can
              efficiently cleave in this intervening region.
  SpCas9      Cas9 from Streptococcus pyogenes
  eCas9       Any engineered Cas9 variant
  eSpCas9     SpCas9 variant bearing alanine substitutions at three positions
              predicted to interact with the non-target DNA strand
  NmCas9      Cas9 from Neisseria meningitidis
  SaCas9      Cas9 from Staphylococcus aureus
  SpCas9-HF1  Alanine substitutions introduced at four residues in SpCas9,
              identified from previously published crystal structures, to
              disrupt non-specific contacts with the phosphate backbone of the
              target DNA strand (which interacts with the gRNA) to create
              SpCas9-HF1 (high-fidelity variant 1).
  Cas9n       Engineered variants of Cas9 in which one of the two nuclease
              domains has been catalytically inactivated, which results in the
              nicking of only one DNA strand and leaving the other strand
              intact. Another strategy proposed to reduce off-target effects
              is to use paired Cas9 nickases (Cas9n), mutated versions of Cas9
              in which one of the two nuclease domains (RuvC or HNH) has been
              catalytically inactivated (for example, by introduction of a
              D10A or H840A mutation). Paired nickases can be directed by two
              gRNAs targeted to neighbouring sites to create offset nicks that
              can induce indel mutations.
  Cpf1        Cpf1 family nucleases

motifs:
  Below are common SPACER>PAM arrangements.
       5'-Motif-3'        Protein  System                          Citation
          N{{20}}>NGG       SpCas9   Streptococcus pyogenes (Sp)     ?
       N{{17,20}}>NGG       SpCas9   Streptococcus pyogenes          ?
          N{{20}}>NGA       SpCas9   Streptococcus pyogenes VQR      Kleinstiver, et al., 2015
          N{{20}}>NGNG      SpCas9   Streptococcus pyogenes EQR      Kleinstiver, et al., 2015
          N{{20}}>NGAG      SpCas9   Streptococcus pyogenes EQR      Kleinstiver, et al., 2015
          N{{20}}>NGCG      SpCas9   Streptococcus pyogenes VRER     Kleinstiver, et al., 2015
          N{{20}}>NAAG      SpCas9   Streptococcus pyogenes QQR1     Anders, et al., 2016
          N{{20}}>NAG       SpCas9   Streptococcus pyogenes          ?
          N{{20}}>NRG       SpCas9   Streptococcus pyogenes          ?
  G{{,2}}N{{19,20}}>NGG       ??Cas9   ?                               ?
         RN{{19}}>NGG       Cas9p    Plants                          Ma, et al., 2015
        RYN{{19}}>NGG       Cas9p    Plants                           + Ma & Liu, 2016
         N{{20?}}>NNAGAAW   StCas9   Streptococcus thermophilus (St) Cong et al., 2013
       N{{20,23}}>NNAGAA    StCas9   Streptococcus thermophilus      Kleinstiver, et al., 2015
          N{{20}}>NGGNG     StCas9   Streptococcus thermophilus      ?
          N{{21}}>NNGRRT    SaCas9   Staphylococcus aureus (Sa)      Ran et al., 2015
       N{{21,23}}>NNGRRT    SaCas9   Staphylococcus aureus           Kleinstiver, et al., 2015
         N{{20?}}>NGRRT     SaCas9   Staphylococcus aureus           ?
         N{{20?}}>NGRRN     SaCas9   Staphylococcus aureus           ?
          N{{20}}>NNNNGMTT  NmCas9   Neisseria meningitidis (Nm)     Hou et al., 2013
          N{{20}}>NNNNACA   CjCas9   Campylobacter jejuni (Cj)       ?
         N{{20?}}>NAAAAC    TdCas9   Treponema denticola (Td)        ?
           TTTN<N{{20,23}}  Cpf1     Acidaminococcus/Lachnospiraceae ?
            TTN<N{{20,23}}  Cpf1     Francisella novicida (putative) ?

outputs:
  STDOUT                            Abbreviated program status
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
  folder/protection-primers.fasta   Primers for amplifying protection dDNAs
  folder/primers.fasta              Primers to check for KO/KI (contains
                                    expected amplicon sizes in header)
""".format(**locals())
__epilog__ = """\
example:
 $ python3 {__program__} --fasta genome.fasta --gff genome.gff --folder addtag-output
""".format(**locals())

class Alignment(object):
    """Class primarily representing a SAM alignment"""
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
        self.prefilter = None
        
        # Variable to hold the scores
        self.score = {}
        
        # Haven't yet added:
        #  lidentities, ridentities, r_score, bae, chari, oof, proxgc, want, xu
        #self.ridentities = nucleotides.ridentities(self.contig_target, aligned_target)
        #self.r_scores = {}
        #for i in [4, 8, 12, 16]:
        #    self.r_scores[i] = scores.r_score(self.contig_target, aligned_target, i)
    
    def calculate_scores(self, parent):
        # parent = (sequence, target, pam, upstream, downstream)
        this = (self.sequence, self.target, self.pam, self.upstream, self.downstream)
        prefilter = []
        for C in algorithms.single_algorithms:
            c_score = C.calculate(this)
            self.score[C.name] = c_score
            if C.prefilter:
                if (C.minimum <= c_score <= C.maximum):
                    prefilter.append(True)
                else:
                    prefilter.append(False)
        for C in algorithms.paired_algorithms:
            c_score = C.calculate(parent, this)
            self.score[C.name] = c_score
            if C.prefilter:
                if (C.minimum <= c_score <= C.maximum):
                    prefilter.append(True)
                else:
                    prefilter.append(False)
        
        # If alignment meets all prefilter criteria, then set it as True
        self.prefilter = all(prefilter) # all([]) returns True
    
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
            'motif=' + self.motif,
            'prefilter=' + str(self.prefilter)] +
            [x + '=' + str(round(self.score[x], 2)) for x in self.score]
            ) + ')'

class Donor(object):
    prefix = 'Donor'
    sequences = {}
    indices = {}
    
    @classmethod
    def generate_fasta(cls, filename, sep=':'):
        """
        Creates a FASTA file
        >id feature:contig:orientation:start1..end1:mAT:start2..end2 ...
        """
        with open(filename, 'w') as flo:
            for sequence, obj in cls.sequences.items():
                #don_entry = tuple(["dDNA-"+str(i), feature, contig, '+', start1, start2, end1, end2, dDNA])
                print(' '.join(['>'+obj.name, 'spacers='+str(len(obj.spacers))] + [obj.format_location(x, sep) for x in obj.locations]), file=flo)
                print(sequence, file=flo)
        print(cls.__name__ + ' dDNA FASTA generated: {!r}'.format(filename))
        return filename
    
    def format_location(self, location, sep=':'):
        feature, contig, orientation, *segments = location
        output_list = [feature, contig, orientation]
        for x in segments:
            if isinstance(x, str):
                x = (x,)
            output_list.append('..'.join(map(str, x)))
        return sep.join(output_list)
    
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
    
    def __repr__(self):
        return self.__class__.__name__ + '(' + ' '.join([
            self.name,
            self.contig + ':' + self.orientation + ':' +
            str(self.segment1[0]) + '..' + str(self.segment1[1]) + ':' +
            self.insert + ':' + str(self.segment2[0]) + '..' + str(self.segment2[1]),
            'features=' + ','.join(self.features) or 'None',
            'spacers=' + str(len(self.spacers))
            ]) + ')'

class ExcisionDonor(Donor):
    prefix = 'exDonor'
    sequences = {}
    indices = {}
    
    @classmethod
    def get_targets(cls, args, sequence):
        targets = set()
        #for seq_i, sequence in enumerate(dDNAs):
        for orientation in ['+', '-']:
            if (orientation == '-'):
                sequence = nucleotides.rc(sequence)
            
            for i in range(len(args.parsed_motifs)):
                spacers, pams, side = args.parsed_motifs[i]
                compiled_regex = args.compiled_motifs[i]
                #matches = nucleotides.motif_search(sequence, spacers, pams, side)
                matches = nucleotides.motif_search2(sequence, side, compiled_regex)
                for seq, start, end, spacer, pam in matches:
                    if (orientation == '-'):
                        start, end = len(sequence) - end, len(sequence) - start
                    filtered_targets = target_filter(seq, args)
                    for filt_seq, filt_spacer, filt_pam in filtered_targets:
                        t_upstream = sequence[start-10:start]
                        t_downstream = sequence[end:end+10]
                        
                        targets.add((orientation, start, end, t_upstream, t_downstream, filt_seq, side, filt_spacer, filt_pam, args.motifs[i]))
        return sorted(targets) # becomes a list
    
    @classmethod
    def generate_donors(cls, args, features, contigs):
        """
        Creates the DNA oligo with the structure:
        [upstream homology][unique gRNA][downstream homology]
        that excises the target feature
        """
        # Generate the full set of potential dDNAs
        for feature in features:
            contig, start, end, strand = features[feature]
            # start & end are 0-based indices, inclusive/exclusive
            
            # assumes start < end
            # DNA 5' of feature is upstream
            # DNA 3' of feature is downstream
            
            my_contig = contigs[contig]
            orientation = '+'
            
            # For each potential dDNA, evaluate how good it is
            for us_trim in range(args.excise_upstream_feature_trim[0], args.excise_upstream_feature_trim[1]+1):
                for ds_trim in range(args.excise_downstream_feature_trim[0], args.excise_downstream_feature_trim[1]+1):
                    for us_hom in range(args.excise_upstream_homology[0], args.excise_upstream_homology[1]+1):
                        for ds_hom in range(args.excise_downstream_homology[0], args.excise_downstream_homology[1]+1):
                            for insert_length in range(args.excise_insert_lengths[0], args.excise_insert_lengths[1]+1):
                                if (args.excise_donor_lengths[0] <= us_hom+insert_length+ds_hom <= args.excise_donor_lengths[1]):
                                    start1, end1 = start-us_hom-us_trim, start-us_trim
                                    start2, end2 = end+ds_trim, end+ds_hom+ds_trim
                                    upstream = my_contig[start1:end1]
                                    downstream = my_contig[start2:end2]
                                    #upstream = contigs[contig][start - args.excise_donor_homology[1]:start]
                                    #downstream = contigs[contig][end:end + args.excise_donor_homology[1]]
                                    
                                    # when insert_length = 0, then the kmers are [''] (single element, empty string)
                                    for mAT in nucleotides.kmers(insert_length):
                                        # Add this candidate dDNA to the list of all candidate dDNAs
                                        dDNA = upstream + mAT + downstream
                                        
                                        #if (args.excise_donor_lengths[0] <= len(dDNA) <= args.excise_donor_lengths[1]):
                                        #dDNAs.append(dDNA)
                                        new_targets = cls.get_targets(args, dDNA) # [(orientation, start, end, filt_seq, side, filt_spacer, filt_pam), ...]
                                        
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
                                                cls(feature, contig, orientation, dDNA, (start1, end1), mAT, (start2, end2), spacer=t)
    
    @classmethod
    def generate_alignments(cls):
        for ind, obj in cls.indices.items():
            print(obj.name)
            loc = next(iter(obj.locations)) # Pull an arbitrary location record
            for segment in loc[3:]:
                if isinstance(segment, str):
                    length = len(segment)
                else:
                    length = segment[1] - segment[0]
                print(cls.make_label(length, ''), end='')
            print('')
            print(obj.sequence)
            for s in obj.spacers:
                orientation = s[0]
                start = s[1]
                end = s[2]
                spacer = s[7]
                pam = s[8]
                if (orientation == '+'):
                    print(' '*start + spacer + pam)
                else:
                    print(' '*start + nucleotides.rc(spacer+pam))
    
    @staticmethod
    def make_label(length, label):
        if (length == 0):
            return ''
        if (length == 1):
            return '|'
        if (length > 1):
            out = ['-'] * length
            out[0] = '['
            out[-1] = ']'
            if (length >= len(label) + 4):
                start = int(length/2 - len(label)/2)
                for i, c in enumerate(label):
                    out[start+i] = c
            return ''.join(out)

class ReversionDonor(Donor):
    prefix = 'reDonor'
    sequences = {}
    indices = {}
    
    @classmethod
    def generate_donors(cls, args, features, contigs):
        """
        Creates the DNA oligo with the structure:
        [upstream homology][original feature][downstream homology]
        """
        # There should be one ReversionDonor for each feature
        for feature in features:
            # First, get the homology blocks up- and down-stream of the feature
            contig, start, end, strand = features[feature]
            
            feature_length = end - start
            
            my_contig = contigs[contig]
            orientation = '+'
            
            for us_hom in range(args.revert_upstream_homology[0], args.revert_upstream_homology[1]+1):
                for ds_hom in range(args.revert_downstream_homology[0], args.revert_downstream_homology[1]+1):
                    if (args.revert_donor_lengths[0] <= us_hom+feature_length+ds_hom <= args.revert_donor_lengths[1]):
                        # to implement:
                        #   make revert_upstream_homology and revert_downstream_homology exclude the trim sequences
                        if (strand == '+'):
                            start1, end1 = start-us_hom, start
                            start2, end2 = end, end+ds_hom
                        else:
                            start1, end1 = start-ds_hom, start
                            start2, end2 = end, end+us_hom
                        upstream = my_contig[start1:end1]
                        downstream = my_contig[start2:end2]
                        feature_sequence = my_contig[start:end]
                        
                        dDNA = upstream + feature_sequence + downstream
                        #for re_seq, re_target in ReversionTarget.indices.items():
                        #    if feature in [x[0] for x in re_target.locations]:
                        #        cls(feature, contig, orientation, (start1, end1), '..'.join(map(str, [start, end])), (start2, end2), dDNA, re_target)
                        #cls(feature, contig, orientation, (start1, end1), '..'.join(map(str, [start, end])), (start2, end2), dDNA, None)
                        cls(feature, contig, orientation, dDNA, (start1, end2), spacer=None)

class Target(object):
    """Data structure defining a gRNA Target"""
    prefix = 'Target'
    sequences = {} # key = nucleotide sequence, value = ExcisionTarget object
    indices = {} # key = exTarget-102, value = ExcisionTarget object
    
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
        
        print(cls.__name__ + ' SAM file parsed: {!r}'.format(filename))
    
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
                print(' '.join(['>'+obj.name] + [obj.format_location(x, sep) for x in obj.locations]), file=flo)
                print(sequence, file=flo)
                
        print(cls.__name__ + ' query FASTA generated: {!r}'.format(filename))
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
                    'alignments=' + str(len([a for a in obj.alignments if a.prefilter])) + '/' + str(len(obj.alignments)),
                    'on-target=' + str(round(obj.score['Azimuth'], 2)),
                    'off-target=' + str(round(obj.off_targets['Hsu-Zhang'], 2)),
                ]), file=flo)
                print(sequence, file=flo)
        print(cls.__name__ + ' spacers FASTA generated: {!r}'.format(filename))
        return filename
    
    def format_location(self, location, sep=':'):
        feature, contig, orientation, start, end, upstream, downstream = location
        return sep.join([feature, contig, orientation, '..'.join([str(start), str(end)])])
    
    def __init__(self, feature, contig, orientation, start, end, upstream, downstream, sequence, side, spacer, pam, motif):
        """Create a structure for holding individual sequence information"""
        location = (feature, contig, orientation, start, end, upstream, downstream)
        self.locations = set()
        
        self.sequence = sequence
        self.side = side
        self.spacer = spacer
        self.pam = pam
        self.motif = motif
        
        # List to store alignments
        self.alignments = []
        
        # Get the index number
        self.index = len(self.sequences)
        self.name = self.prefix + '-' + str(self.index)
        
        the_target = self.sequences.setdefault(self.sequence, self)
        the_target.locations.add(location)
        self.indices.setdefault(the_target.name, the_target)
        
        # Scores for this sequence only (not PairedSequenceAlgorithm)
        self.score = {}
        self.off_targets = {}
        if (len(self.locations) > 0):
            self.calculate_default_scores()
    
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
    
    def overlap_coverage(self, start1, end1, start2, end2, index_base=0):
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
    
    def calculate_default_scores(self):
        """Populate this scores for this Sequence"""
        loc = next(iter(self.locations)) # Pull an arbitrary location record
        #parent = (self.sequence, self.target, self.pam, self.upstream, self.downstream)
        parent = (self.sequence, self.spacer, self.pam, loc[5], loc[6])
        for C in algorithms.single_algorithms:
            if (C.default != None):
                self.score[C.name] = C.default
            else:
                self.score[C.name] = C.calculate(parent)
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
        for i in range(len(args.parsed_motifs)):
            spacers, pams, side = args.parsed_motifs[i]
            compiled_regex = args.compiled_motifs[i]
            m = nucleotides.motif_conformation2(sequence, side, compiled_regex)
            if m:
                r_spacer = m[0]
                r_pam = m[1]
                r_motif = args.motifs[i]
                break
            else:
                if (side == '>'):
                    l = max(map(len, pams))
                    r_spacer = sequence[:-l]
                    r_pam = sequence[-l:]
                    r_motif = args.motifs[i]
                elif (side == '<'):
                    l = max(map(len, pams))
                    r_spacer = seq[:l]
                    r_pam = seq[l:]
                    r_motif = args.motifs[i]
        return r_spacer, r_pam, r_motif
    
    def get_features(self, features, contig, start, end):
        the_features = []
        for feature in features:
            f_contig, f_start, f_end, f_strand = features[feature]
            if (contig == f_contig):
                if (self.overlap_coverage(start, end, f_start, f_end) > 0):
                    the_features.append(feature)
        
        return the_features
    
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
            batch_scores = C.calculate(queries)
            
            # Assign the score to the appropriate Target
            for i, t_index in enumerate(t_sorted):
                t_obj = cls.indices[t_index]
                t_obj.score[C.name] = batch_scores[i]
    
    def score_off_targets(self, args, homologs, features):
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
        on_targets = {}
        off_targets = {}
        for C in calculators:
            on_targets[C.name] = []
            off_targets[C.name] = []
        
        # Get list of on-target features
        on_target_features = set([x[0] for x in self.locations])
        if homologs:
            for f in list(on_target_features):
                on_target_features.update(homologs.get(f, set()))
        
        # Check each alignment
        for a in self.alignments:
            a_features = self.get_features(features, a.contig, a.start, a.end) # contig, start, end
            
            if a.prefilter:
                for i, C in enumerate(calculators):
                    c_score = a.score[C.name]
                    if (C.minimum <= c_score <= C.maximum):
                        # if alignment is NOT a target:
                        if ((a_features == None) or (len(on_target_features.intersection(a_features)) == 0)):
                            off_targets[C.name].append(c_score)
                        # if alignment is a target
                        else:
                            on_targets[C.name].append(c_score)
        
        # Perform off-target calculations
        for C in calculators:
            try:
                self.off_targets[C.name] = scores.off_target_score(off_targets[C.name], on_targets[C.name])
            except ZeroDivisionError:
                self.off_targets[C.name] = 0.0 # This happens if the reversion-gRNA targets the excision-dDNA
    
    def __repr__(self):
        """Return a string containing a printable representation of the Target object."""
        return self.__class__.__name__ + '(' + ' '.join([
            self.name,
            self.spacer + '|' + self.pam,
            'motif=' + self.motif,
            'locations=' + str(len(self.locations)),
            'alignments=' + str(len([a for a in self.alignments if a.prefilter])) + '/' + str(len(self.alignments))] +
            [x + '=' + str(round(self.score[x], 2)) for x in self.score] + 
            ['OT:' + x + '=' + str(round(self.off_targets[x], 2)) for x in self.off_targets]
            ) + ')'

class ExcisionTarget(Target):
    prefix = 'exTarget'
    # Non-redundant dicts
    sequences = {} # key = nucleotide sequence, value = ExcisionTarget object
    indices = {} # key = exTarget-102, value = ExcisionTarget object
    
    @classmethod
    def get_targets(cls, args, contigs, features):
        """
        Searches within the annotated features on all contigs for targets
        that match the args.motifs criteria, then filters them.
        
        Returns list of valid gRNA sites with the following format
          [(feature, contig, start, end, seq, side, spacer, pam), ...]
        """
        # Find unique gRNA sites within each feature
        # Use a sliding window to make a list of queries
        for feature in features:
            
            feature_contig, feature_start, feature_end, feature_strand = features[feature]
            if (feature_end == None):
                feature_end = len(contigs[feature_contig])
            
            # Make sure the contig the feature is on is present in the FASTA
            if feature_contig in contigs:
                # Find a site within this feature that will serve as a unique gRNA
                
                # for each orientation:
                # if (args.strands in ['+', 'both']):
                #  targets.extend...
                # if (args.strands in ['-', 'both']):
                #  targets.extend...
                
                # Search both the '+' and '-' strands
                for orientation in ['+', '-']:
                    if (orientation == '+'):
                        sequence = contigs[feature_contig][feature_start:feature_end]
                    else:
                        sequence = nucleotides.rc(contigs[feature_contig][feature_start:feature_end])
                    
                    for i in range(len(args.parsed_motifs)):
                        spacers, pams, side = args.parsed_motifs[i]
                        compiled_regex = args.compiled_motifs[i]
                        #matches = nucleotides.motif_search(sequence, spacers, pams, side)
                        matches = nucleotides.motif_search2(sequence, side, compiled_regex)
                        for seq, start, end, spacer, pam in matches:
                            if (orientation == '+'):
                                real_start = feature_start + start
                                real_end = feature_start + end
                            else:
                                real_start = len(sequence) - end + feature_start
                                real_end = len(sequence) - start + feature_start
                            filtered_targets = target_filter(seq, args)
                            for filt_seq, filt_spacer, filt_pam in filtered_targets:
                                t_upstream = sequence[start-10:start]
                                t_downstream = sequence[end:end+10]
                                cls(feature, feature_contig, orientation, real_start, real_end, t_upstream, t_downstream, filt_seq, side, filt_spacer, filt_pam, args.motifs[i])
                                # Maybe add to ReversionDonor here

class ReversionTarget(Target):
    prefix = 'reTarget'
    # Non-redundant dicts
    sequences = {} # key = nucleotide sequence, value = ReversionTarget object
    indices = {} # key = reTarget-234, value = ReversionTarget object
    
    @classmethod
    def get_targets(cls):
        for dDNA, obj in ExcisionDonor.sequences.items():
            obj_features = ','.join([x[0] for x in list(obj.locations)])
            
            # Populate the ReversionTarget sequences indices dicts
            for t in obj.spacers: # [(orientation, start, end, filt_seq, side, filt_spacer, filt_pam), ...]
                #final_targets.append(tuple([obj.name, obj_feature, obj_contig] + list(t)))
                # t = (orientation, start, end, t_upstream, t_downstream, filt_seq, side, filt_spacer, filt_pam, args.motifs[i])
                ReversionTarget(obj_features, obj.name, t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9])

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

def parse_arguments():
    # Create the argument parser
    parser = argparse.ArgumentParser(
        description=__description__,
        epilog=__epilog__,
        formatter_class=CustomHelpFormatter
    )
    
    # Add required arguments
    required_group = parser.add_argument_group('required arguments')
    required_group.add_argument("--fasta", required=True, metavar="*.fasta", type=str,
        help="FASTA file with contigs to find unique gRNA sites. Input FASTA \
            must not be compressed. Ambiguous bases within the FASTA will not \
            be chosen for gRNA sites. All FASTA sequences should have unique \
            primary headers (everything between the '>' symbol and the first \
            whitespace should be unique).")
    required_group.add_argument("--gff", required=True, metavar="*.gff", type=str,
        help="GFF file specifying chromosomal features")
    required_group.add_argument("--folder", required=True, metavar="FOLDER",
        type=str, help="Path of folder to store generated files")
    
    # Special version action optional argument
    parser.add_argument("-v", "--version", action='version', version='{__program__} {__version__}'.format(**globals()))
    
    # Add optional arguments
    #parser.add_argument("--feature_homolog_regex", metavar="REGEX", type=str, default=None, help="regular expression with capturing group containing invariant feature. Example: '(.*)_[AB]' will treat features C2_10010C_A and C2_10010C_B as homologs")
    # okay idea, but needs more thought before implementation
    parser.add_argument("--feature_homologs", metavar="*.homologs", type=str, default=None,
        help="Path to text file containing homologous features on the same \
            line, separated by TAB characters")
    #parser.add_argument("--pams", metavar="SEQ", nargs="+", type=str,
    #    default=["NGG"], help="Constrain finding only targets with these PAM sites")
    #parser.add_argument("--target_lengths", nargs=2, metavar=('MIN', 'MAX'),
    #    type=int, default=[17, 20],
    #    help="The length range of the 'target'/'spacer'/gRNA site")
    # Replacement for --pams and --target_lengths with this:
    parser.add_argument("--motifs", metavar="MOTIF", nargs="+", type=str,
        default=["N{17,20}>NGG"],
        help="Find only targets with these 'SPACER>PAM' motifs, written from \
        5' to 3'. '>' points toward PAM. IUPAC ambiguities accepted. '{a,b}' \
        are quantifiers. Be sure to enclose motif parameters in quotes so your \
        shell does not interpret STDIN/STDOUT redirection.")
    # Need to decide if construct inputs should be TSV, or FASTA
    # And whether or not there should be an upstream parameter separate from
    # a downstream one. or if they are the same, then what?
    parser.add_argument("--constructs", metavar="*.fasta", nargs="+", type=str,
        default=[], help="The first sequence will be prepended, and the second \
        sequence will be appended to the generated spacer sequences to form \
        the construct sequences. It is useful to put the gRNA promotor as the \
        first sequence, and the scaffold sequence and terminator as the \
        second. Specify one FASTA file for each motif.")
    parser.add_argument("--tag", metavar='TAG', type=str, default='ID',
        help="GFF tag with feature names. Examples: 'ID', 'Name', 'Parent', or 'locus_tag'")
    parser.add_argument("--ambiguities", type=str, choices=["discard", "disambiguate", "keep"], default="discard",
        help="How generated gRNAs should treat ambiguous bases: \
        discard - no gRNAs will be created where the FASTA has an ambiguous base; \
        disambiguate - gRNAs containing ambiguous bases will be converted to a set of non-ambiguous gRNAs; \
        keep - gRNAs can have ambiguous bases")
    parser.add_argument("--case", type=str, default="ignore",
        choices=["ignore", "discard-lower", "discard-upper", "invariant-lower", "invariant-upper"],
        help="Restrict generation of gRNAs based on case of nucleotides in input FASTA")
    #parser.add_argument("--min_contig_edge_distance", metavar="N", type=int, default=500,
    #    help="Minimum distance from contig edge a site can be found")
    parser.add_argument("--features", metavar="FEATURE", type=str, nargs="+", default=["gene"],
        help="Features to design gRNAs against. Must exist in GFF file. Examples: 'CDS', 'gene', 'mRNA', 'exon'")
    parser.add_argument("--target_gc", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[25, 75],
        help="Generated gRNA spacers must have %%GC content between these values (excludes PAM motif)")
    parser.add_argument("--excise_upstream_homology", nargs=2, metavar=("MIN", "MAX"), type=int, default=[50,50],
        help="Range of homology lengths acceptable for knock-out dDNAs, inclusive.")
    parser.add_argument("--excise_downstream_homology", nargs=2, metavar=("MIN", "MAX"), type=int, default=[47,50],
        help="Range of homology lengths acceptable for knock-out dDNAs, inclusive.")
    parser.add_argument("--excise_donor_lengths", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[100, 100],
        help="Range of lengths acceptable for knock-out dDNAs, inclusive.")
    parser.add_argument("--excise_insert_lengths", nargs=2, metavar=("MIN", "MAX"), type=int, default=[0,3],
        help="Range for inserted DNA lengths, inclusive (mini-AddTag, mAT). If MIN < 0, then regions of dDNA homology (outside the feature) will be removed.")
    parser.add_argument("--excise_feature_edge_distance", metavar="N", type=int, default=0,
        help="If positive, gRNAs won't target any nucleotides within this distance \
             from the edge of the feature. If negative, gRNAs will target nucleotides \
             this distance outside the feature.")
    parser.add_argument("--excise_upstream_feature_trim", nargs=2, metavar=('MIN', 'MAX'),
        type=int, default=[0, 0], help="Between MIN and MAX number of nucleotides \
        upstream of the feature will be considered for knock-out when designing \
        donor DNA.")
    parser.add_argument("--excise_downstream_feature_trim", nargs=2, metavar=("MIN", "MAX"),
        type=int, default=[0, 0], help="Between MIN and MAX number of nucleotides \
        downstream of the feature will be considered for knock-out when designing \
        donor DNA.")
    #parser.add_argument("--min_donor_insertions", metavar="N", type=int, default=2,
    #    help="The uniqueness of final donor DNA compared to the rest of the genome")
    #parser.add_argument("--min_donor_deletions", metavar="N", type=int, default=2,
    #    help="The uniqueness of final donor DNA compared to the rest of the genome")
    #parser.add_argument("--min_donor_substitutions", metavar="N", type=int, default=2,
    #    help="The uniqueness of final donor DNA compared to the rest of the genome")
    #parser.add_argument("--min_donor_errors", metavar="N", type=int, default=3,
    #    help="The uniqueness of final donor DNA compared to the rest of the genome")
    parser.add_argument("--revert_upstream_homology", nargs=2, metavar=("MIN", "MAX"), type=int, default=[300,300],
        help="Range of homology lengths acceptable for knock-in dDNAs, inclusive.")
    parser.add_argument("--revert_downstream_homology", nargs=2, metavar=("MIN", "MAX"), type=int, default=[300,300],
        help="Range of homology lengths acceptable for knock-in dDNAs, inclusive.")
    parser.add_argument("--revert_donor_lengths", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[0, 100000],
        help="Range of lengths acceptable for knock-in dDNAs.")
    parser.add_argument("--min_donor_distance", metavar="N", type=int, default=36,
        help="The minimum distance in bp a difference can exist from the edge of donor DNA") # homology with genome
    parser.add_argument("--max_consecutive_ts", metavar="N", type=int, default=4,
        help="The maximum number of Ts allowed in generated gRNA sequences")
    # program currently will only search 'both' strands
    #parser.add_argument("--strands", type=str, choices=["+", "-", "both"], default="both",
    #    help="Strands to search for gRNAs")
    
    # Add command line arguments for the additional hard constraints:
    #  Only report potential targets that have no off targets with mismatches within 8, 12, N nt from 3' end
    parser.add_argument("--processors", metavar="N", type=int, default=(os.cpu_count() or 1),
        help="Number of processors to use when performing pairwise sequence alignments")
    parser.add_argument("--aligner", type=str, choices=['addtag', 'blast+',
        'blat', 'bowtie', 'bowtie2', 'bwa', 'cas-offinder'], default='bowtie2',
        help="Program to calculate pairwise alignments. Please note that the 'addtag' internal aligner is very slow.")
    # Other aligners to consider: 'rmap', 'maq', 'shrimp2', 'soap2', 'star', 'rhat', 'mrsfast', 'stampy'
    parser.add_argument("--python2_path", type=str, default="python",
        help="Path to the Python 2.7+ program")
    parser.add_argument("--bowtie_path", type=str, default="bowtie",
        help="Path to the 'bowtie' executable")
    parser.add_argument("--bowtie-build_path", type=str, default="bowtie-build",
        help="Path to the 'bowtie-build' executable")
    parser.add_argument("--bowtie2_path", type=str, default="bowtie2",
        help="Path to the 'bowtie2' executable")
    parser.add_argument("--bowtie2-build_path", type=str, default="bowtie2-build",
        help="Path to the 'bowtie2-build' executable")
    parser.add_argument("--bwa_path", type=str, default="bwa",
        help="Path to the 'bwa' executable")
    parser.add_argument("--blastn_path", type=str, default="blastn",
        help="Path to the 'blastn' executable")
    parser.add_argument("--blat_path", type=str, default="blat",
        help="Path to the 'blat' executable")
    parser.add_argument("--test", action="store_true", help="Perform tests only")
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Parse the motifs
    # Compile regex for motifs
    args.parsed_motifs = []
    args.compiled_motifs = []
    for motif in args.motifs:
        spacers, pams, side = parse_motif(motif) # Parse the motif
        args.parsed_motifs.append((spacers, pams, side)) # Add to args
        args.compiled_motifs.append(nucleotides.compile_motif_regex(spacers, pams, side, anchored=False)) # Add to args
    
    # Return the parsed arguments
    return args

def _parse_motif_helper(submotif):
    """
    Helper function that parses either the SPACER or PAM motif.
    Decodes quantifiers and returns a list
    """
    # Keep track of expanded sequences
    sequences = ['']
    
    # Keep track if a quantifier is being parsed
    quantifier = None
    
    # Iterate through the characters
    for c in submotif:
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
            for i, s in enumerate(sequences):
                for length in range(min_length, max_length+1):
                    new_sequences.append(s + last_chars[i]*length)
            sequences = new_sequences
            quantifier = None
        elif (quantifier != None): # add current character to quantifier if it is open
            quantifier += c
        else: # add the current character to the expanded sequences
            for i in range(len(sequences)):
                sequences[i] = sequences[i] + c
    
    return sequences

def parse_motif(motif):
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
    
    return _parse_motif_helper(spacer_motif), _parse_motif_helper(pam_motif), side

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

def merge_features(features):
    """Combine overlapping features?"""
    return features

def target_filter(seq, args):
#def new_target_filter(seq, side, spacer, pam, args):
    '''
    Filters the candidate gRNA sequence based on the following criteria:
     1) case: ignore, discard-lower, discard-upper (does not process invariant-lower/invariant-upper)
     2) ambiguous character expansion: discard, keep, disambiguate
     3) maximum consecutive Ts
     4) SPACER>PAM check using regex
     5) %GC
    
    Returns list of validated sequences
    '''
    
    # Check the case of the potential gRNA sequence
    if (args.case == "discard-lower"):
        if regex.search('[a-z]', seq):
            return [] # Reject this sequence because it has lower-case characters
    elif (args.case == "discard-upper"):
        if regex.search('[A-Z]', seq):
            return [] # Reject this sequence because it has upper-case characters
    # if (args.case == "ignore"), then do nothing
    
    # Convert sequences to upper-case so it can be evaluated by the scoring algorithms
    seq = seq.upper()
    
    # Check if target sequence has any ambiguities
    if (args.ambiguities == 'discard'):
        if regex.search('[^ATCGatcg]', seq):
            # Reject this sequence
            return []
        seqs = [seq]
    # Disambiguate sequences if necessary
    elif (args.ambiguities == 'disambiguate'):
        seqs = nucleotides.disambiguate_iupac(seq)
    # Do nothing if just 'keep'
    else:
        seqs = [seq]
    
    # Remove targets with T{5,}
    seqs = [ nt for nt in seqs if ('T'*(args.max_consecutive_ts+1) not in nt) ]
    
    # Remove targets that do not confine to at least one of the defined
    # SPACER>PAM motifs.
    # Also, separate the SPACER and PAM motifs
    temp_seqs = []
    temp_targets = []
    temp_pams = []
    for nt in seqs:
        for i in range(len(args.parsed_motifs)):
        #for spacers, pams, side in args.parsed_motifs:
            #m = nucleotides.motif_conformation(nt, spacers, pams, side)
            spacers, pams, side = args.parsed_motifs[i]
            compiled_regex = args.compiled_motifs[i]
            m = nucleotides.motif_conformation2(nt, side, compiled_regex)
            if m:
                temp_seqs.append(nt)
                temp_targets.append(m[0])
                temp_pams.append(m[1])
                #break # consider breaking here
    #seqs = temp_seqs
    
    # Remove targets whose %GC is outside the chosen bounds
    temp_seqs2 = []
    temp_targets2 = []
    temp_pams2 = []
    for i in range(len(temp_seqs)):
        if (args.target_gc[0] <= scores.gc_score(temp_targets[i]) <= args.target_gc[1]):
            temp_seqs2.append(temp_seqs[i])
            temp_targets2.append(temp_targets[i])
            temp_pams2.append(temp_pams[i])
    #seqs = temp_seqs2
    
    # [(seq, spacer, pam), (seq, spacer, pam)]
    rets = []
    for i in range(len(temp_seqs2)):
        rets.append((temp_seqs2[i], temp_targets2[i], temp_pams2[i]))
    return rets

def index_reference(args):
    if (args.aligner == 'addtag'):
        index_file = fasta
    elif (args.aligner == 'blast+'):
        pass
    elif (args.aligner == 'blat'):
        pass
    elif (args.aligner == 'bowtie'):
        pass
    elif (args.aligner == 'bowtie2'):
        index_file = bowtie2.index_reference(args.fasta, tempdir=args.folder, threads=args.processors)
    elif (args.aligner == 'bwa'):
        pass
    elif (args.aligner == 'cas-offinder'):
        pass
    return index_file

def align(query_file, index_file, args):
    if (args.aligner == 'addtag'):
        sam_file = None
    elif (args.aligner == 'blast+'):
        pass
    elif (args.aligner == 'blat'):
        pass
    elif (args.aligner == 'bowtie'):
        pass
    elif (args.aligner == 'bowtie2'):
        sam_file = bowtie2.align(query_file, index_file, folder=args.folder, threads=args.processors)
    elif (args.aligner == 'bwa'):
        pass
    elif (args.aligner == 'cas-offinder'):
        pass
    return sam_file

def main():
    """Function to run complete AddTag analysis"""
    
    # Obtain command line arguments and parse them
    args = parse_arguments()
    
    if args.test:
        # Perform test code
        test(args)
    else:
        # Get timestamp
        start = time.time()
        
        # Echo the command line parameters to STDOUT
        print(args)
        
        # Create the project directory if it doesn't exist
        os.makedirs(args.folder, exist_ok=True)
        
        # Load the FASTA file specified on the command line
        contigs = utils.load_fasta_file(args.fasta)
        
        # Open and parse the GFF file specified on the command line
        features = utils.load_gff_file(args.gff, args.features, args.tag)
        
        # Make index of homologs
        if args.feature_homologs:
            homologs = utils.load_homologs(args.feature_homologs)
        else:
            homologs = None
        
        # Merge features?
        #features = merge_features(features)
        
        # Search features within contigs for targets that match the motifs
        ExcisionTarget.get_targets(args, contigs, features)
        
        # Index the reference FASTA
        index_file = index_reference(args)
        
        # Write the query list to FASTA
        ex_query_file = ExcisionTarget.generate_query_fasta(os.path.join(args.folder, 'excision-query.fasta'))
        
        # Use selected alignment program to find all matches in the genome
        ex_sam_file = align(ex_query_file, index_file, args)
        
        print("ExcisionTarget before SAM parsing")
        for et_seq, et_obj in ExcisionTarget.sequences.items():
            print(et_obj)
        
        # Load the SAM file and add Alignments to ExcisionTarget sequences
        ExcisionTarget.load_sam(os.path.join(args.folder, 'excision-query.sam'), args, contigs)
        
        # Calculate off-target/guide scores for each algorithm
        print("ExcisionTarget after SAM parsing and off-target scoring")
        for et_seq, et_obj in ExcisionTarget.sequences.items():
            et_obj.score_off_targets(args, homologs, features)
            print(et_obj)
            for a in et_obj.alignments:
                print('  ', a)
        
        # Batch calculate with new ExcisionTarget class
        ExcisionTarget.score_batch()
        
        print("ExcisionTarget after Azimuth calculation")
        for et_seq, et_obj in ExcisionTarget.sequences.items():
            print(et_obj)
        
        # Generate the FASTA with the final scores
        excision_spacers_file = ExcisionTarget.generate_spacers_fasta(os.path.join(args.folder, 'excision-spacers.fasta'))
        
        # Discard potential gRNAs that have mismatches with their target site
        #if (args.case == "invariant-lower"):
        #    pass
        #elif (args.case == "invariant-upper"):
        #    pass
        
        # Generate excision dDNAs and their associated reversion gRNA spacers
        ExcisionDonor.generate_donors(args, features, contigs)
        ReversionTarget.get_targets()
        ex_dDNA_file = ExcisionDonor.generate_fasta(os.path.join(args.folder, 'excision-dDNAs.fasta'))
        re_query_file = ReversionTarget.generate_query_fasta(os.path.join(args.folder, 'reversion-query.fasta'))
        
        # Use selected alignment program to find all matches in the genome
        re_sam_file = align(re_query_file, index_file, args)
        
        # Load the SAM file and add Alignments to ReversionTarget sequences
        ReversionTarget.load_sam(os.path.join(args.folder, 'reversion-query.sam'), args, contigs)
        
        # Calculate off-target/guide scores for each algorithm
        print("ReversionTarget after SAM parsing and off-target scoring")
        for re_seq, re_obj in ReversionTarget.sequences.items():
            re_obj.score_off_targets(args, homologs, features)
            print(re_obj)
            for a in re_obj.alignments:
                print('  ', a)
        
        # Batch calculate with new ReversionTarget class
        ReversionTarget.score_batch()
        
        print("ReversionTarget after Azimuth calculation")
        for rt_seq, rt_obj in ReversionTarget.sequences.items():
            print(rt_obj)
        
        # Generate the FASTA with the final scores
        reversion_spacers_file = ReversionTarget.generate_spacers_fasta(os.path.join(args.folder, 'reversion-spacers.fasta'))
        
        # Generate reversion dDNAs and write them to FASTA
        ReversionDonor.generate_donors(args, features, contigs)
        re_dDNA_file = ReversionDonor.generate_fasta(os.path.join(args.folder, 'reversion-dDNAs.fasta'))
        
        # Test code to generate alignments
        ExcisionDonor.generate_alignments()
        
        # Print time taken for program to complete
        print('Runtime: {}s'.format(time.time()-start))

def test(args):
    """Code to test the classes and functions in 'source/__init__.py'"""
    # Echo the command line parameters
    print(args)
    
    sys.exit(10)
    
    # Get timestamp
    start = time.time()
    
    # Load the FASTA file specified on the command line
    print("=== FASTA ===")
    #contigs = utils.load_fasta_file(args.fasta)
    #print(list(contigs.keys())[:5])
    
    # Test SAM file parsing
    print("=== SAM ===")
    
    # Open and parse the GFF file specified on the command line
    # returns a dictionary:
    #  features[ID] = (contig, start(bp), end(bp), frame)
    print("=== GFF ===")
    #features = utils.load_gff_file(args.gff, args.features, args.tag)
    #for k in list(features.keys())[:10]:
    #    print(k, features[k])
    
    # Test code to find all similar oligonucleotides in the FASTA
    print("=== Align ===")
    #target = 'TCCGGTACAKTGAKTTGTAC'
    #regex = nucleotides.build_regex(target, max_errors=2)
    #matches = nucleotides.find_target_matches(regex, contigs, overlap=True)
    #for m in matches:
    #    print(m)
    #    for seq in nucleotides.disambiguate_iupac(m[4]):
    #        print(seq, len(seq), hsuzhang.hsuzhang_score(target, seq))
    
    # Test score calculations
    print("=== Score calculations ===")
    
    # Print time taken for program to complete
    print('Runtime: {}s'.format(time.time()-start))
