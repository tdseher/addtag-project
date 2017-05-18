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
  Program for identifying unique endogenous gRNA sites 
  and creating unique synthetic gRNA sites.
  
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
       5'-Motif-3'        Protein  Origin system                   Citation
       N{{17,20}}>NGG       SpCas9   Streptococcus pyogenes (Sp)     ?
          N{{20}}>NGG       SpCas9   Streptococcus pyogenes          ?
          N{{20}}>NGA       SpCas9   Streptococcus pyogenes VQR      ?
          N{{20}}>NAG       SpCas9   Streptococcus pyogenes          ?
          N{{20}}>NRG       SpCas9   Streptococcus pyogenes          ?
          N{{20}}>NGCG      SpCas9   Streptococcus pyogenes VRER     ?
  G{{,2}}N{{19,20}}>NGG       ??Cas9   ?                               ?
         N{{20?}}>NNAGAAW   StCas9   Streptococcus thermophilus (St) Cong et al., 2013
          N{{20}}>NNAGAA    StCas9   Streptococcus thermophilus      ?
          N{{20}}>NGGNG     StCas9   Streptococcus thermophilus      ?
          N{{21}}>NNGRRT    SaCas9   Staphylococcus aureus (Sa)      Ran et al., 2015
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

class Sequence(object):
    """Data structure defining a sequence"""
    
    #substitution_threshold = 4
    substitution_threshold = 15 # There is a bug in the function that counts these errors
    insertion_threshold = 2
    deletion_threshold = 2
    error_threshold = 5
    
    def __init__(self, feature, contig_sequence, args, contig=None, contig_orientation='+', contig_start=None, contig_end=None, feature_orientation=None, contig_upstream='', contig_downstream=''):
        """Create a structure for holding individual sequence information"""
        self.feature = feature
        self.feature_orientation = feature_orientation
        
        self.contig_sequence = contig_sequence
        self.contig_target, self.contig_pam, self.contig_motif = self.split_spacer_pam(self.contig_sequence, args)
        self.contig_upstream = contig_upstream
        self.contig_downstream = contig_downstream
        
        #self.contig_target, self.contig_pam = nucleotides.split_target_sequence(self.contig_sequence, pams, force=True)
        #self.disambiguated_sequences = disambiguate_iupac(self.contig_sequence, kind="dna") # need to apply pre-filters
        self.contig = contig
        self.contig_orientation = contig_orientation
        self.contig_start = contig_start
        self.contig_end = contig_end
        
        # List to store alignments
        self.alignments = []
        
        # Scores for this sequence only (not PairedSequenceAlgorithm)
        self.score = {}
        self.calculate_default_scores()
        
        self.off_targets = {}
    
    def calculate_default_scores(self):
        """Populate this scores for this Sequence"""
        #parent = (self.sequence, self.target, self.pam, self.upstream, self.downstream)
        parent = (self.contig_sequence, self.contig_target, self.contig_pam, self.contig_upstream, self.contig_downstream)
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
    
    def calculate_alignment_scores(self):
        parent = (self.contig_sequence, self.contig_target, self.contig_pam, self.contig_upstream, self.contig_downstream)
        for a in self.alignments:
            a.calculate_scores(parent)
    
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
    
    def add_alignment(self, aligned_sequence, args, aligned_contig, aligned_start, aligned_end, aligned_orientation, aligned_upstream, aligned_downstream):
        """Add a genomic position to the list of alignments"""
        #parent = (self.sequence, self.target, self.pam, self.upstream, self.downstream)
        parent = (self.contig_sequence, self.contig_target, self.contig_pam, self.contig_upstream, self.contig_downstream)
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
        
    def old_add_alignment(self, aligned_sequence, args, aligned_contig, aligned_start, aligned_end, aligned_orientation, aligned_upstream, aligned_downstream):
        """Add a genomic position to the list of alignments"""
        aligned_target, aligned_pam, aligned_motif = self.split_spacer_pam(aligned_sequence, args)
        if ((aligned_target != None) and (aligned_pam != None)):
            #aligned_target, aligned_pam = nucleotides.split_target_sequence(aligned_sequence, pams, force=True)
            substitutions, insertions, deletions = nucleotides.count_errors(self.contig_sequence, aligned_sequence)
            seq = (
                aligned_sequence,
                aligned_target,
                aligned_pam,
                aligned_contig,
                aligned_start,
                aligned_end,
                aligned_orientation,
                substitutions,
                insertions,
                deletions,
                doench.on_target_score_2014(aligned_target, aligned_pam),
                doench.on_target_score_2016(self.contig_target, aligned_target, aligned_pam),
                hsuzhang.hsuzhang_score(self.contig_target, aligned_target, iupac=False),
                scores.linear_score(self.contig_target, aligned_target),
                housden.housden_score(aligned_target),
                morenomateos.morenomateos_score(aligned_target, aligned_pam),
                0, #azimuth.azimuth_score(aligned_target, aligned_pam, upstream=aligned_upstream, downstream=aligned_downstream),
                nucleotides.ridentities(self.contig_target, aligned_target),
                scores.r_score(self.contig_target, aligned_target, 4),
                scores.r_score(self.contig_target, aligned_target, 8),
                scores.r_score(self.contig_target, aligned_target, 12),
                scores.r_score(self.contig_target, aligned_target, 16),
            )
            self.alignments.append(seq)
        else:
            print('Cannot add alignment:', aligned_sequence, aligned_contig, aligned_start, aligned_end, aligned_orientation, file=sys.stderr)
    
    def get_features(self, features, contig, start, end):
        the_features = []
        for feature in features:
            f_contig, f_start, f_end, f_strand = features[feature]
            if (contig == f_contig):
                if (self.overlap_coverage(start, end, f_start, f_end) > 0):
                    the_features.append(feature)
        
        return the_features
    
    def overlap_distance(self, start1, end1, start2, end2):
        coverage = self.overlap_coverage(start1, end1, start2, end2)
        if (coverage > 0):
            return -coverage
        else:
            #return min(abs(start2-end1), abs(start1-end2))
            return max(start2-end1, start1-end2)
            # 1 ..........
            # 2              XXXXXX
            
            # 1          ........
            # 2 XXXXXX
    
    def overlap_coverage(self, start1, end1, start2, end2):
        coverage = 0
        
        if (start2 <= start1 <= end1 <= end2):
            # 1      .............
            # 2   XXXXXXXXXXXXXXXXXX
            coverage = abs(end1 - start1) + 1
        elif (start1 <= start2 <= end2 <= end1):
            # 1      .....................
            # 2            XXXXXXXXXX
            coverage = abs(end2 - start2) + 1
        elif (start1 <= start2 <= end1 <= end2):
            # 1      ........
            # 2        XXXXXXXX
            coverage = abs(end1 - start2) + 1
        elif (start2 <= start1 <= end2 <= end1):
            # 1      .............
            # 2  XXXXXXXXXX
            coverage = abs(end2 - start1) + 1
            
        return coverage
    
    def overlap_percent_coverage(self, start1, end1, start2, end2):
        coverage = self.overlap_coverage(start1, end1, start2, end2)
        return float(coverage) / max(abs(end1 - start1), abs(end2 - start2))
    
    def percent_full_length_coverage(self, start1, end1, len1, start2, end2, len2):
        coverage = self.overlap_coverage(start1, end1, start2, end2)
        return float(coverage) / max(len1, len2)
    
    def score_off_targets(self, args, homologs, features):
        """
        Calculate Guide Score (off-target score) for all single/paired
        algorithms with the 'off_target=True' attribute.
        """
        calculators = []
        
        for C in algorithms.single_algorithms:
            if C.off_target:
                calculators.append(C)
        for C in algorithms.paired_algorithms:
            if C.off_target:
                calculators.append(C)
        #for C in algorithms.batched_single_algorithms:
        #    if C.off_target:
        #        calculators.append(C)
        
        on_targets = {}
        off_targets = {}
        for C in calculators:
            on_targets[C.name] = []
            off_targets[C.name] = []
        
        on_target_features = set([self.feature])
        if homologs:
            on_target_features.update(homologs.get(self.feature, set()))
        
        for a in self.alignments:
            #a_features = self.get_features(features, a[3], a[4], a[5]) # contig, start, end
            a_features = self.get_features(features, a.contig, a.start, a.end) # contig, start, end
            
            # If the divergence is reasonable
            if ((a.score['substitutions'] <= self.substitution_threshold) and
                (a.score['insertions'] <= self.insertion_threshold) and
                (a.score['deletions'] <= self.deletion_threshold) and 
                (a.score['errors'] <= self.error_threshold)
            ):
                for i, C in enumerate(calculators):
                    #c_score = a[10+i]
                    c_score = a.score[C.name]
                    if (C.minimum <= c_score <= C.maximum):
                        # if alignment is NOT a target:
                        if ((a_features == None) or (len(on_target_features.intersection(a_features)) == 0)):
                            off_targets[C.name].append(c_score)
                        # if alignment is a target
                        else:
                            on_targets[C.name].append(c_score)
        
        for C in calculators:
            try:
                self.off_targets[C.name] = scores.off_target_score(off_targets[C.name], on_targets[C.name])
            except ZeroDivisionError:
                self.off_targets[C.name] = 0.0 # This happens if the reversion-gRNA targets the excision-dDNA
        
    def old_score(self, args, homologs, features):
        """Calculate Guide scores for each algorithm"""
        calculations = [
            'doench2014',
            'doench2016',
            'hsuzhang',
            'linear',
            'housden',
            'morenomateos',
            #'azimuth'
        ]
        on_targets = {}
        off_targets = {}
        for calculation in calculations:
            on_targets[calculation] = []
            off_targets[calculation] = []
        # on_targets = {
        #     'doench2014': [self.doench2014],
        #     'doench2016': [self.doench2016],
        #     'hsuzhang': [self.hsuzhang],
        #     'linear': [self.linear],
        #     'housden': [self.housden],
        #     'morenomateos': [self.morenomateos],
        #     'azimuth': [self.azimuth],
        # }
        
        on_target_features = set([self.feature])
        if homologs:
            on_target_features.update(homologs.get(self.feature, set()))
        #print(self.feature, on_target_features, file=sys.stderr)
        
        for a in self.alignments:
            a_features = self.get_features(features, a[3], a[4], a[5]) # contig, start, end
            
            # If the divergence is reasonable
            if ((a[7] <= self.substitution_threshold) and
                (a[8] <= self.insertion_threshold) and
                (a[9] <= self.deletion_threshold) and 
                (sum(a[7:10]) <= self.error_threshold)
            ):
                for i, calculation in enumerate(calculations):
                    c_score = a[10+i]
                    if (c_score >= self.score_thresholds[calculation]):
                        # if alignment is NOT a target:
                        if ((a_features == None) or (len(on_target_features.intersection(a_features)) == 0)):
                            off_targets[calculation].append(c_score)
                        # if alignment is a target
                        else:
                            on_targets[calculation].append(c_score)
                    #else:
                    #    on_targets[calculation].append(0)
            #if (a[3:6] != (self.contig, self.contig_start, self.contig_end)):
                # lambda a, b: all([a[0]<=b[0], a[1]<=b[1], a[2]<b[2], sum(a)<=b[3]])
                    # if (a[10] >= self.doench2014_threshold):
                    #     doench2014_list.append(a[10])
                    # if (a[11] >= self.doench2016_threshold):
                    #     doench2016_list.append(a[11])
                    # if (a[12] >= self.hsuzhang_threshold):
                    #     hsuzhang_list.append(a[12])
                    # if (a[13] >= self.linear_threshold):
                    #     linear_list.append(a[13])
                    # if (a[15] >= self.morenomateos_threshold):
                    #     morenomateos_list.append(a[15])
                    # if (a[16] >= self.azimuth_threshold):
                    #     azimuth_list.append(a[16])
        # self.off_target_doench2014 = scores.off_target_score(doench2014_list, (self.doench2014,))
        # self.off_target_doench2016 = scores.off_target_score(doench2016_list, (self.doench2016,))
        # self.off_target_hsuzhang = scores.off_target_score(hsuzhang_list, (self.hsuzhang,))
        # self.off_target_linear = scores.off_target_score(linear_list, (self.linear,))
        # self.off_target_morenomateos = scores.off_target_score(morenomateos_list, (self.morenomateos,))
        # self.off_target_azimuth = 0 #scores.off_target_score(azimuth_list, (self.azimuth))
        
        for calculation in calculations:
            try:
                self.off_targets[calculation] = scores.off_target_score(off_targets[calculation], on_targets[calculation])
            except ZeroDivisionError:
                self.off_targets[calculation] = 0.0
        # try:
        #     self.off_target_doench2014 = scores.off_target_score(off_targets['doench2014'], on_targets['doench2014'])
        #     self.off_target_doench2016 = scores.off_target_score(off_targets['doench2016'], on_targets['doench2016'])
        #     self.off_target_hsuzhang = scores.off_target_score(off_targets['hsuzhang'], on_targets['hsuzhang'])
        #     self.off_target_linear = scores.off_target_score(off_targets['linear'], on_targets['linear'])
        #     self.off_target_housden = scores.off_target_score(off_targets['housden'], on_targets['housden'])
        #     self.off_target_morenomateos = scores.off_target_score(off_targets['morenomateos'], on_targets['morenomateos'])
        #     self.off_target_azimuth = 0 #scores.off_target_score(off_targets['azimuth'], on_targets['azimuth'])
        # except ZeroDivisionError:
        #     print(self, file=sys.stderr)
        #     for i in calculations:
        #         print(i, off_targets[i], on_targets[i], file=sys.stderr)
    
    def __repr__(self):
        """Return a string containing a printable representation of the Sequence object."""
        return 'Sequence(feature=' + self.feature + \
            ', ' + self.contig + \
            ':' + str(self.contig_start) + '..' + str(self.contig_end) + \
            ', ' + self.contig_target + '|' + self.contig_pam + \
            ', motif=' + self.contig_motif + \
            ', alignments=' + str(len(self.alignments)) + \
            ', ' + ', '.join([x + '=' + str(round(self.score[x], 2)) for x in self.score]) + \
            ', ' + ', '.join(['OT:' + x + '=' + str(round(self.off_targets[x], 2)) for x in self.off_targets]) + \
            ')'
            #', azimuth=' + str(round(self.azimuth, 2)) + \
            #', off-target=' + str(round(self.off_targets['hsuzhang'], 2)) + \
            

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
        for C in algorithms.single_algorithms:
            self.score[C.name] = C.calculate(this)
        for C in algorithms.paired_algorithms:
            self.score[C.name] = C.calculate(parent, this)
        
        # Since the distance algorithms aren't working correctly, temporarily
        # hard code these to 0
        self.score['substitutions'] = 0
        self.score['insertions'] = 0
        self.score['deletions'] = 0
        self.score['errors'] = 0
    
    def get_upstream_sequence(self, length, contigs):
        """Returns upstream sequence"""
        return ''
    
    def get_downstream_sequence(self, length, contigs):
        """Returns downstream sequence"""
        return ''
    
    def __repr__(self):
        return 'Alignment(' + \
            self.contig + ':' + self.orientation + ':' + str(self.start) + '..' + str(self.end) + \
            ', ' + self.target + '|' + self.pam + \
            ', motif=' + self.motif + \
            ', ' + ', '.join([x + '=' + str(round(self.score[x], 2)) for x in self.score]) + \
            ')'

class excision_dDNA(object):
    def __init__(self, contig, segment1, segment2, insert):
        self.contig = contig,
        self.segment1 = segment1 # (start1, end1)
        self.segment2 = segment2 # (start2, end2)
        self.insert = insert
    # >dDNA-0 C2_10210C_B:Ca22chr2B_C_albicans_SC5314:+:2096471..2096521:CCA:2097397:2097444
    

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

def load_sam_file_test(filename, args, contigs, sep=':'):
    """Read in SAM file.
    sep is the separator for the header. Positions are converted to 0-index
    Creates a list of Sequence objects
    """
    
    alignments = {}
    with open(filename, 'r') as flo:
        for line in flo:
            if not line.startswith('@'):
                sline = line.rstrip().split("\t")
                if ((len(sline) > 5) and (sline[2] != '*')):
                    # >feature:contig:orientation:start..end
                    feature, source_contig, source_orientation, source_start_end = sline[0].split(sep)
                    source_start, source_end = source_start_end.split('..')
                    source_start = int(source_start)
                    source_end = int(source_end)
                    source = (feature, source_contig, source_orientation, source_start, source_end)
                    source_upstream = contigs[source_contig][source_start-10:source_start]
                    source_downstream = contigs[source_contig][source_end:source_end+10]
                    if (source_orientation == '-'):
                        source_upstream, source_downstream = nucleotides.rc(source_downstream), nucleotides.rc(source_upstream)
                    
                    # Get orientation
                    alignment_orientation = utils.sam_orientation(int(sline[1]))
                    
                    # Get alignment position
                    alignment_contig = sline[2]
                    alignment_start = int(sline[3])-1
                    alignment_end = int(sline[3])-1+utils.cigar_length(sline[5])
                    
                    # Reverse-complement if needed
                    alignment_sequence = contigs[alignment_contig][alignment_start:alignment_end]
                    alignment_upstream = contigs[alignment_contig][alignment_start-10:alignment_start]
                    alignment_downstream = contigs[alignment_contig][alignment_end:alignment_end+10]
                    actual_sequence = sline[9]
                    if (alignment_orientation == '-'):
                        alignment_sequence = nucleotides.rc(alignment_sequence)
                        alignment_upstream, alignment_downstream = nucleotides.rc(alignment_downstream), nucleotides.rc(alignment_upstream)
                        actual_sequence = nucleotides.rc(actual_sequence)
                    
                    # if source not in alignments:
                    #    alignments[source] = s
                    # alignments[sournce].add_alignment(...)
                    
                    # Assuming creating an instance of Sequence() is cheaper
                    # than traversing alignments dict()
                    s = Sequence(
                        feature,
                        actual_sequence,
                        # contigs[source_contig][int(source_start):int(source_end)], # contig_sequence
                        args,
                        contig=source_contig,
                        contig_orientation=source_orientation,
                        contig_start=int(source_start),
                        contig_end=int(source_end),
                        feature_orientation=None,
                        contig_upstream=source_upstream,
                        contig_downstream=source_downstream,
                    )
                    
                    #dbg1=0
                    #try:
                    #    dbg1 = len(alignments[source].alignments)
                    #except KeyError:
                    #    pass
                    alignments.setdefault(source, s).add_alignment(
                        alignment_sequence, # aligned_sequence (as when matched with reference, thus may be revcomp of initial query)
                        args,
                        alignment_contig, # aligned_contig
                        alignment_start, # aligned_start
                        alignment_end, # aligned_end
                        alignment_orientation, # aligned_orientation (+/-)
                        alignment_upstream,
                        alignment_downstream,
                    )
                    #dbg2 = len(alignments[source].alignments)
                    #if (dbg1 == dbg2):
                    #    print("problem: ", source, file=sys.stderr)
    
    print('SAM file parsed: {!r}'.format(filename))
    #return list(alignments.values()) # unsorted list
    return list(map(lambda x: alignments[x], sorted(alignments))) # sorted

def load_reversion_sam_file_test(filename, args, contigs, donors, sep=':'):
    """Read in SAM file.
    sep is the separator for the header. Positions are converted to 0-index
    Creates a list of Sequence objects
    """
    
    dDNA_dict = {}
    for line in donors:
        #don_entry = tuple(["dDNA-"+str(i), feature, contig, '+', start1, start2, end1, end2, dDNA])
        dDNA_id, feature, contig, orientation, start1, end1, mAT, start2, end2, seq = line
        dDNA_dict[dDNA_id] = line
    
    alignments = {}
    with open(filename, 'r') as flo:
        for line in flo:
            if not line.startswith('@'):
                sline = line.rstrip().split("\t")
                if ((len(sline) > 5) and (sline[2] != '*')):
                    # >id:orientation:start..end:feature:contig
                    dDNA_id, dDNA_orientation, dDNA_start_end, feature, source_contig = sline[0].split(sep)
                    dDNA_start, dDNA_end = dDNA_start_end.split('..')
                    dDNA_start = int(dDNA_start)
                    dDNA_end = int(dDNA_end)
                    dDNA = (feature, dDNA_id, dDNA_orientation, dDNA_start, dDNA_end)
                    
                    dDNA_list = dDNA_dict[dDNA_id]
                    source_contig = dDNA_list[2]
                    source_orientation = dDNA_orientation
                    source_start = dDNA_list[4]
                    source_end = dDNA_list[8]
                    
                    # >feature:contig:orientation:start..end
                    #feature, source_contig, source_orientation, source_start_end = sline[0].split(sep)
                    #source_start, source_end = source_start_end.split('..')
                    #source_start = int(source_start)
                    #source_end = int(source_end)
                    #source = (feature, source_contig, source_orientation, source_start, source_end)
                    source_upstream = contigs[source_contig][source_start-10:source_start]
                    source_downstream = contigs[source_contig][source_end:source_end+10]
                    dDNA_upstream = source_upstream + dDNA_list[9][:dDNA_start]
                    dDNA_downstream = dDNA_list[9][dDNA_end:] + source_downstream
                    if (source_orientation == '-'):
                        source_upstream, source_downstream = nucleotides.rc(source_downstream), nucleotides.rc(source_upstream)
                        dDNA_upstream, dDNA_downstream = nucleotides.rc(dDNA_downstream), nucleotides.rc(dDNA_upstream)
                    
                    # Get orientation
                    alignment_orientation = utils.sam_orientation(int(sline[1]))
                    
                    # Get alignment position
                    alignment_contig = sline[2]
                    alignment_start = int(sline[3])-1
                    alignment_end = int(sline[3])-1+utils.cigar_length(sline[5])
                    
                    # Reverse-complement if needed
                    alignment_sequence = contigs[alignment_contig][alignment_start:alignment_end]
                    alignment_upstream = contigs[alignment_contig][alignment_start-10:alignment_start]
                    alignment_downstream = contigs[alignment_contig][alignment_end:alignment_end+10]
                    actual_sequence = sline[9]
                    if (alignment_orientation == '-'):
                        alignment_sequence = nucleotides.rc(alignment_sequence)
                        alignment_upstream, alignment_downstream = nucleotides.rc(alignment_downstream), nucleotides.rc(alignment_upstream)
                        actual_sequence = nucleotides.rc(actual_sequence)
                    
                    # if source not in alignments:
                    #    alignments[source] = s
                    # alignments[sournce].add_alignment(...)
                    
                    # Assuming creating an instance of Sequence() is cheaper
                    # than traversing alignments dict()
                    s = Sequence(
                        feature,
                        actual_sequence,
                        # contigs[source_contig][int(source_start):int(source_end)], # contig_sequence
                        args,
                        contig=dDNA_id, #source_contig,
                        contig_orientation=dDNA_orientation, #source_orientation,
                        contig_start=dDNA_start, #int(source_start),
                        contig_end=dDNA_end, #int(source_end),
                        feature_orientation=None,
                        contig_upstream=dDNA_upstream, #source_upstream,
                        contig_downstream=dDNA_downstream, #source_downstream,
                    )
                    #s.dDNA_id = dDNA_id
                    
                    alignments.setdefault(dDNA, s).add_alignment(
                        alignment_sequence, # aligned_sequence (as when matched with reference, thus may be revcomp of initial query)
                        args,
                        alignment_contig, # aligned_contig
                        alignment_start, # aligned_start
                        alignment_end, # aligned_end
                        alignment_orientation, # aligned_orientation (+/-)
                        alignment_upstream,
                        alignment_downstream,
                    )
    
    print('SAM file parsed: {!r}'.format(filename))
    return list(map(lambda x: alignments[x], sorted(alignments))) # sorted

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
    parser.add_argument("--revert_donor_lengths", nargs=2, metavar=('MIN', 'MAX'), type=int, default=[300, 600],
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

# Template 1
def generate_excise_target(args, feature):
    """Finds the gRNA sequence to cut within the specified places
    on the feature. May target either the + or - strand.
    """
    pass

def generate_revert_target(args, feature):
    """Creates the gRNA sequence that targets the excise donor DNA oligo.
    May target either the + or - strand.
    """
    pass

def generate_excise_donor(args, features, contigs):
    """
    Creates the DNA oligo with the structure:
    [upstream homology][unique gRNA][downstream homology]
    that excises the target feature
    """
    final_targets = []
    final_dDNAs = []
    
    # Generate the full set of potential dDNAs
    for feature in features:
        # First, get the homology blocks up- and down-stream of the feature
        contig, start, end, strand = features[feature]
        # start & end are 0-based indices, inclusive
        
        # assumes start < end
        # DNA 5' of feature is upstream
        #upstream = contigs[contig][start-args.excise_donor_homology[1]:start] # This is the max homology length
        
        # DNA 3' of feature is downstream
        #downstream = contigs[contig][end:end+args.excise_donor_homology[1]] # Max homology length
        
        # mini_addtags = itertools.something()
        # dDNAs = [upstream + x + downstream for x in itertools.something()]
        
        #dDNAs = []
        all_targets = {}
        #targets = set()
        my_contig = contigs[contig]
        
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
                                    new_targets = get_targets_temp(args, dDNA)
                                    
                                    for t in new_targets:
                                        if t in all_targets:
                                            if (len(dDNA) < len(all_targets[t][0])):
                                                all_targets[t] = (dDNA, start1, end1, mAT, start2, end2)
                                        else:
                                            all_targets[t] = (dDNA, start1, end1, mAT, start2, end2)
        
        for i, t in enumerate(all_targets):
            # Each entry should be:
            # (feature, contig, orientation, start, end, seq, side, spacer, pam)
            dDNA, start1, end1, mAT, start2, end2 = all_targets[t]
            rev_entry = tuple(["dDNA-"+str(i), feature, contig] + list(t))
            don_entry = tuple(["dDNA-"+str(i), feature, contig, '+', start1, end1, mAT, start2, end2, dDNA])
            
            # This non-redundancy code is REALLY REALLY REALLY slow
            #if (entry not in final_targets):
            #    final_targets.append(rev_entry)
            #    final_dDNAs.append(don_entry)
            final_targets.append(rev_entry)
            final_dDNAs.append(don_entry)
    
    return final_targets, final_dDNAs
    #    # For each potential dDNA, evaluate how good it is
    #    for insert_length in range(args.excise_insert_lengths[0], args.excise_insert_lengths[1]+1):
    #        # when insert_length = 0, then the kmers are [''] (single element, empty string)
    #        print(' ', insert_length, file=sys.stderr, end='', flush=True)
    #        for mAT in nucleotides.kmers(insert_length):
    #            # Add this candidate dDNA to the list of all candidate dDNAs
    #            dDNAs.append(upstream + mAT + downstream)
    #            
    #            # Use code to generate targets for this revised region
    #            # between [maximum 5' distance] upstream mAT downstream [maximum 3' distance]
    #            #targets.extend(get_targets(args, pff_contigs, pff_features))
    #            
    #            targets.extend(get_targets_temp(args, dDNAs))
    #    print('', file=sys.stderr)
    #    for t in targets:
    #        print(t)
    
    #return [(feature, seq_i, orientation, start, end, filt_seq, side, filt_spacer, filt_pam) for seq_i, orientation, start, end, filt_seq, side, filt_spacer, filt_pam in targets ]
    
#    # Write the targets to a FASTA file
#    query_file = utils.generate_query(os.path.join(args.folder, 'reversion-query.fasta'), targets)
#    
#    # Use selected alignment program to find all potential off-targets in the genome
#    sam_file = align(query_file, index_file, args.folder, args.processors)
    
    # Calculate scores for each target
    

def generate_revert_donor(args, feature):
    """Use template DNA sequence to create oligo that will be used for fixing
    the DSB
    """
    # Optionally provide a list of restriction enzyme targeting sites
    # on either side for easier cloning
    pass

def format_output():
    # output should take the following format
    # contig    gRNA-start-pos   gRNA-end-pos    gRNA-strand   features    gRNA-sequence    PAM     on-target-score    off-target-score
    pass

def merge_features(features):
    """Combine overlapping features?
    """
    return features

def process(args):
    # The Cas9 cuts 3-4bp upstream of the PAM sequence:
    #  ======= gRNA ==== === PAM
    #  CGATGCATCGACTTTAC CGA AGG
    #                   ^ Cut
    
    features = []
    for f in features:
        # identify_pam_positions()
        et = generate_excise_target(args, f)
        rt = generate_revert_target(args, f)
        
        ed = generate_excise_donor(args, f, rt)
        rd = generate_revert_donor(args, f)
    
    # Cut site can be anywhere within target feature
    # genome          ACGGATTAGAGAGAGGCCTCCTCCCGGAT-GAATGGAAGACTAAACGGTAGATATAAG
    # feature         -------|ORF->                                <-ORF|-------
    # PAM                                              PAM         PAM
    # excise target                         CCCGGAT^GAATGG
    # excise genome   ACGGATTAGAGAGAGGCCTCCTCCCGGAT^GAATGGAAGACTAAACGGTAGATATAAG
    # excise donor    ACGGATT-----------------------------------AAACGGTAGATATAAG
    # revert target      GATT-----------------------------------AAACGG
    # revert donor     CGGATTAGAGAGAGGCCTCCTCCCGGAT-GAATGGAAGACTAAACGGTAGATATA
    
    pass

# Template 2
def find_candidate_targets(args):
    pass

def align_targets():
    pass

def find_targets():
    pass

def make_primers():
    pass

def make_donors():
    pass

def excise():
    candidate_targets = find_targets()
    
    excise_primers = make_primers()
    
    excise_donors, revert_targets = make_donors()

def revert(revert_targets):
    revert_donors = make_donors()
    
    revert_primers = make_primers()

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

def get_targets_temp(args, sequence):
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
                    #targets.add((seq_i, orientation, start, end, filt_seq, side, filt_spacer, filt_pam))
                    targets.add((orientation, start, end, filt_seq, side, filt_spacer, filt_pam))
    return sorted(targets) # becomes a list

def get_targets(args, contigs, features):
    """
    Searches within the annotated features on all contigs for targets
    that match the args.motifs criteria, then filters them.
    
    Returns list of valid gRNA sites with the following format
      [(feature, contig, start, end, seq, side, spacer, pam), ...]
    """
    # Ideally, this code would procedurally write to the query file
    # as new sequences were being added, thus limiting the amount of memory
    # used to store sequences.
    
    #compiled_motif_regexes = []
    #for spacers, pams, side in args.parsed_motifs:
    #    compiled_motif_regexes.append(nucleotides.compile_motif_regex(spacers, pams, side, anchored=False))
    
    # Find unique gRNA sites within each feature
    targets = set()
    # Use a sliding window to make a list of queries
    for feature in features:
        
        feature_contig, feature_start, feature_end, feature_strand = features[feature]
        if (feature_end == None):
            feature_end = len(contigs[feature_contig])
        
        # Make sure the contig the feature is on is present in the FASTA
        if feature_contig in contigs:
            # print(feature, features[feature], file=sys.stderr)
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
                #for spacers, pams, side in args.parsed_motifs:
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
                            targets.add((feature, feature_contig, orientation, real_start, real_end, filt_seq, side, filt_spacer, filt_pam))
    return sorted(targets) # becomes a list


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
        
        # Generate the query list: [(feature, contig, start, end, sequence), ...]
        targets = get_targets(args, contigs, features)
        
        # Index the reference FASTA
        index_file = index_reference(args)
        
        # Write the query list to FASTA
        query_file = utils.generate_excision_query(os.path.join(args.folder, 'excision-query.fasta'), targets)
        
        # Use selected alignment program to find all matches in the genome
        ex_sam_file = align(query_file, index_file, args)
        
        # Open the SAM file
        ex_alignments = load_sam_file_test(os.path.join(args.folder, 'excision-query.sam'), args, contigs)
        
        # Calculate off-target/guide scores for each algorithm
        for s in ex_alignments:
            s.score_off_targets(args, homologs, features)
        
        # make list of all sequences to calculate Azimuth score on
        queries = []
        for s in ex_alignments:
            queries.append((s.contig_sequence, s.contig_target, s.contig_pam, s.contig_upstream, s.contig_downstream))
        for C in algorithms.batched_single_algorithms:
            batch_scores = C.calculate(queries)
            for i, s in enumerate(ex_alignments):
                s.score[C.name] = batch_scores[i]
        
        # Generate the FASTA with the final scores
        excision_spacers_file = utils.generate_excision_spacers(os.path.join(args.folder, 'excision-spacers.fasta'), ex_alignments, sep=':')
        
        # Print the sequences
        #for s in alignments:
        #    print(s)
        #    for a in s.alignments:
        #        print('  ', a)
        
        # Discard potential gRNAs that have mismatches with their target site
        #if (args.case == "invariant-lower"):
        #    pass
        #elif (args.case == "invariant-upper"):
        #    pass
        
        revert_targets, excise_donors = generate_excise_donor(args, features, contigs)
        
        revert_query_file = utils.generate_reversion_query(os.path.join(args.folder, 'reversion-query.fasta'), revert_targets)
        excise_dDNA_file = utils.generate_donor(os.path.join(args.folder, 'excision-dDNAs.fasta'), excise_donors)
        
        
        
        
        
        # Use selected alignment program to find all matches in the genome
        re_sam_file = align(revert_query_file, index_file, args)
        
        # Open the SAM file
        re_alignments = load_reversion_sam_file_test(os.path.join(args.folder, 'reversion-query.sam'), args, contigs, excise_donors)
        
        # Calculate off-target/guide scores for each algorithm
        for s in re_alignments:
            s.score_off_targets(args, homologs, features)
        
        # make list of all sequences to calculate Azimuth score on
        queries = []
        for s in re_alignments:
            queries.append((s.contig_sequence, s.contig_target, s.contig_pam, s.contig_upstream, s.contig_downstream))
        for C in algorithms.batched_single_algorithms:
            batch_scores = C.calculate(queries)
            for i, s in enumerate(re_alignments):
                s.score[C.name] = batch_scores[i]
        
        
        
        
        # Generate the FASTA with the final scores
        reversion_spacers_file = utils.generate_excision_spacers(os.path.join(args.folder, 'reversion-spacers.fasta'), re_alignments, sep=':')
        
        
        # Print the sequences
        for s in re_alignments:
            print(s)
            for a in s.alignments:
                print('  ', a)
        
        
        
        # Print time taken for program to complete
        print('Runtime: {}s'.format(time.time()-start))

def test(args):
    """Code to test the classes and functions in 'source/__init__.py'"""
    # Echo the command line parameters
    print(args)
    
    sys.exit()
    
    # Get timestamp
    start = time.time()
    
    # Load the FASTA file specified on the command line
    print("=== FASTA ===")
    contigs = utils.load_fasta_file(args.fasta)
    print(list(contigs.keys())[:5])
    
    # Test SAM file parsing
    print("=== SAM ===")
    try:
        alignments = load_sam_file_test(os.path.join(args.folder, 'excision-query.sam'), args, contigs)
        for s in alignments:
            print(s)
            for a in s.alignments:
                print('  ', a)
    except FileNotFoundError:
        print('Skipping...')
    
    # Open and parse the GFF file specified on the command line
    # returns a dictionary:
    #  features[ID] = (contig, start(bp), end(bp), frame)
    print("=== GFF ===")
    features = utils.load_gff_file(args.gff, args.features, args.tag)
    for k in list(features.keys())[:10]:
        print(k, features[k])
    
    # Test code to find all similar oligonucleotides in the FASTA
    print("=== Align ===")
    target = 'TCCGGTACAKTGAKTTGTAC'
    regex = nucleotides.build_regex(target, max_errors=2)
    matches = nucleotides.find_target_matches(regex, contigs, overlap=True)
    for m in matches:
        print(m)
        for seq in nucleotides.disambiguate_iupac(m[4]):
            print(seq, len(seq), hsuzhang.hsuzhang_score(target, seq))
    
    # Test Hsu score
    print("=== Hsu 2013 ===")
    a = 'CGATGGCTWGGATCGATTGAC'
    b = 'AAGTGCTCTTAAGAGAAATTC'
    c = 'ATGSCTCGGATCGATTGAC'
    print(hsuzhang.calcHitScore(a, a), hsuzhang.hsuzhang_score(a, a), hsuzhang.hsuzhang_score(a, a, True))
    print(hsuzhang.calcHitScore(a, b), hsuzhang.hsuzhang_score(a, b), hsuzhang.hsuzhang_score(a, b, True))
    print(hsuzhang.calcHitScore(a, c), hsuzhang.hsuzhang_score(a, c), hsuzhang.hsuzhang_score(a, c, True))
    
    # Test Doench 2014 score:
    print("=== Doench 2014 ===")
    pam = 'AGG'
    gRNAa = 'CGATGGCTTGGATCGATTGA'
    gRNAb = 'CGTTGGCTTGGATCGATTGA'
    gRNAc = 'CGATGGCTTCGATCGATTGA'
    gRNAd = 'CGATGGCTTCGAGCGATTGA'
    print(doench.on_target_score_2014(gRNAa, pam))
    print(doench.on_target_score_2014(gRNAb, pam))
    print(doench.on_target_score_2014(gRNAc, pam))
    print(doench.on_target_score_2014(gRNAd, pam))
    
    # Test Doench 2016 score
    print("=== Doench 2016 ===")
    pam = 'AGG'
    gRNAa = 'CGATGGCTTGGATCGATTGA'
    gRNAb = 'CGTTGGCTTGGATCGATTGA'
    gRNAc = 'CGATGGCTTCGATCGATTGA'
    gRNAd = 'CGATGGCTTCGAGCGATTGA'
    print(doench.on_target_score_2016(gRNAa, gRNAa, pam))
    print(doench.on_target_score_2016(gRNAa, gRNAb, pam))
    print(doench.on_target_score_2016(gRNAa, gRNAc, pam))
    print(doench.on_target_score_2016(gRNAa, gRNAd, pam))
    
    # Print time taken for program to complete
    print('Runtime: {}s'.format(time.time()-start))