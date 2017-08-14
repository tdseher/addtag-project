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

# Import non-standard packages
import regex

# Import included AddTag-specific modules
from . import utils
from . import nucleotides
from . import scores
from . import algorithms
from . import aligners

# Define meta variables
__author__ = "Thaddeus D. Seher (@tdseher) & Aaron Hernday"
__date__ = utils.load_git_date()
__fullversion__ = utils.load_git_version()
__version__ = __fullversion__[:7]
__revision__ = utils.load_git_revision()
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
  short    {__version__}
  full     {__fullversion__}
  revision {__revision__}
  date     {__date__}

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

glossary:
  RGN         RNA-guided nuclease, like Cas9, Cas3, and Cpf1.
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
  scaffold    The sequence within the gRNA that is responsible for Cas9 binding.
              Does not include the 20 nt spacer/targeting sequence that is used
              to guide Cas9 to target DNA.
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
  gRNA        Guide RNA is a single molecule with two regions: the spacer and
              the scaffold. It is a synthetic fusion of the endogenous bacterial
              crRNA and tracrRNA sequences, and provides both targeting
              specificity and scaffolding/binding ability for Cas9 nuclease.
              Does not exist in nature. Also referred to as sgRNA.
  sgRNA       Synthetic guide RNA, or single guide RNA (synonymous with 'gRNA').
              Combines the tracrRNA and crRNA, which are separate molecules,
              into a single RNA construct, simplifying the components needed to
              use CRISPR/Cas9 for genome editing (for plasmid or IVT
              expression). A linker loop sequence is included between the two.
  Cas         CRISPR-associated family of genes, which typically couple
              a nuclease, helicase, or polymerase domain with a 
              poly-nucleotide binding domain.
  Cas9        Cas9 family nucleases.
  eCas9       Any engineered Cas9 variant. Usually non-synonymous substitutions
              are placed at one or more residues predicted to interact with the
              non-target DNA strand.
                For instance, SpCas9-HF1 (high-fidelity variant 1) has alanine
                substitutions at four residues in SpCas9, identified from
                crystal structures, in order to disrupt non-specific contacts
                with the phosphate backbone of the target DNA strand (which
                interacts with the gRNA).
  dCas9       Catalytically 'dead' Cas9 protein, that drive RNA-DNA
              hybridization but fail to cleave the target DNA.
  FokI-dCas9  dCas9 fused with the dimerization-dependent FokI nuclease domain:
              creates a dimeric RNA-guided FokI-dCas9 nuclease (RFN)
              architecture requiring recognition of extended double-length
              target sites for efficient cleavage. Amino-terminal fusions of
              FokI to dCas9 can recognize two 20-nucleotide 'half-sites' in a
              'PAM-out' orientation separated by a 13-18 bp spacer and can
              efficiently cleave in this intervening region.
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
  Cpf1        Cpf1 family nucleases.
  Sp          Streptococcus pyogenes.
  St          Streptococcus thermophilus.
  Sm          Streptococcus mutans.
  Sa          Staphylococcus aureus.
  Nm          Neisseria meningitidis.
  Cj          Campylobacter jejuni.
  Td          Treponema denticola.
  Bl          Brevibacillus laterosporus.
  Pm          Pasteurella multocida.
  Fn          Francisella novicida.
  Ec          Escherichia coli.
  
motifs:
  Below are common SPACER>PAM arrangements (* = putative):
       5'-Motif-3'        Protein     System                          Citation
          N{{20}}>NGG       SpCas9      Streptococcus pyogenes (Sp)     ?
       N{{17,20}}>NGG       SpCas9      Streptococcus pyogenes          Fu, et al (2014)
          N{{20}}>NGA       SpCas9 VQR  H. sapiens/E.coli               Kleinstiver, et al (2015a)
          N{{20}}>NGNG      SpCas9 EQR  H. sapiens/E.coli               Kleinstiver, et al (2015a)
          N{{20}}>NGAG      SpCas9 EQR  H. sapiens/E.coli               Kleinstiver, et al (2015a)
          N{{20}}>NGCG      SpCas9 VRER H. sapiens/E.coli               Kleinstiver, et al (2015a)
          N{{20}}>NAAG      SpCas9 QQR  Streptococcus pyogenes QQR1     Anders, et al (2016)
          N{{20}}>NAG       SpCas9      H. sapiens; Cell-free           Hsu, et al (2013)
          N{{20}}>NRG       SpCas9      H. sapiens; Cell-free           Hsu, et al (2013); Karvelis, et al (2015)
         GN{{19}}>NRG       SpCas9      H. sapiens; Cell-free           Hsu, et al (2013)
  G{{,2}}N{{19,20}}>NGG       Cas9p       Plants                          Ma & Liu (2016)
         RN{{19}}>NGG       Cas9p       Plants                          Ma, et al (2015)
        RYN{{19}}>NGG       Cas9p       Plants                           + Ma & Liu (2016)
         N{{20?}}>NNAGAAW   StCas9      Streptococcus thermophilus (St) Horvath, et al (2008); Cong et al (2013)
          N{{20}}>NNAAAAW   StCas9      Streptococcus thermophilus      Fonfara, et al (2013)
       N{{20,23}}>NNAGAA    StCas9      H. sapiens/E.coli               Kleinstiver, et al (2015a)
          N{{20}}>NGGNG     StCas9      Streptococcus thermophilus      Horvath, et al (2008)
          N{{20}}>NHRBMAW   StCas9      Streptococcus thermophilus      Karvelis, et al (2015)
          N{{20}}>NGGWG     StCas9      Saccharomyces cerevisiae        Xu, et al (2015)
          N{{20}}>NGG       SmCas9      Streptococcus  mutans (Sm)      Fonfara, et al (2013)
          N{{21}}>NNGRRT    SaCas9      Staphylococcus aureus (Sa)      Ran, et al (2015)
       N{{21,23}}>NNGRRT    SaCas9      H. sapiens/E.coli               Kleinstiver, et al (2015a)
       N{{21,23}}>NNNRRT    SaCas9 KKH  Staphylococcus aureus KKH       Kleinstiver, et al (2015b)
          N{{20}}>NNNNGATT  NmCas9      Neisseria meningitidis (Nm)     ?
          N{{20}}>NNNNGMTT  NmCas9      Neisseria meningitidis          Hou, et al (2013)
          N{{20}}>NNNNACA   CjCas9      Campylobacter jejuni (Cj)       Fonfara, et al (2013)
          N{{20}}>NNNNRYAC  cjCas9      Campylobacter jejuni            ?
         N{{20?}}>NAAAAC    TdCas9      Treponema denticola (Td)        Zhang (unpublished)
       N{{18,21}}>NGGNCNDD  BlCas9      Brevibacillus laterosporus (Bl) Karvelis, et al (2015)
          N{{20}}>NNNNCND   BlCas9      Brevibacillus laterosporus      Karvelis, et al (2015)
          N{{20}}>NNNNCNDD  BlCas9      Brevibacillus laterosporus      Karvelis, et al (2015)
          N{{20}}>GNNNCNNA  PmCas9      Pasteurella multocida (Pm)      Fonfara, et al (2013)
          N{{20}}>NG        FnCas9      Francisella novicida (Fn) (*)   Fonfara, et al (2013)
           TTTN<N{{20,23}}  Cpf1        Acidaminococcus/Lachnospiraceae ?
            TTN<N{{20,23}}  FnCpf1      Francisella novicida (*)        ?
             AW<GN{{31,32}} EcCas3      Escherichia coli (Ec) (*)       Swarts, et al (2012)
            AWG<N{{32,33}}  EcCas3      Escherichia coli (*)            Swarts, et al (2012)

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
  folder/protection-primers.fasta   Primers for amplifying protection dDNAs
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
        'parsed_motif',
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
    
    def __init__(self, sequence, target, pam, motif, parsed_motif, contig, start, end, orientation, upstream='', downstream=''):
        # Most attributes derived from SAM output
        self.sequence = sequence
        self.target = target
        self.pam = pam
        self.motif = motif
        self.parsed_motif = parsed_motif
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
            c_score = C.calculate(this, parsed_motif=self.parsed_motif)
            self.score[C.name] = c_score
            if C.postfilter:
                if (C.minimum <= c_score <= C.maximum):
                    postfilter.append(True)
                else:
                    postfilter.append(False)
        for C in algorithms.paired_algorithms:
            c_score = C.calculate(parent, this, parsed_motif=self.parsed_motif)
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
            'motif=' + self.motif,
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
        return sorted(set(x[0] for x in self.locations))
    
    def get_location_features(self):
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
                        
                        targets.add((orientation, start, end, t_upstream, t_downstream, filt_seq, side, filt_spacer, filt_pam, args.motifs[i], tuple([tuple(x) if isinstance(x, list) else x for x in args.parsed_motifs[i]])))
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
    
    def get_trims(self, features):
        """
        Get the length of the upstream and downstream trims for each location
        """
        trims = []
        for loc in self.locations:
            contig, start, end, strand = features[loc[0]]
            #mAT = loc[4]
            left_trim = start - loc[3][1]
            right_trim = loc[5][0] - end
            #trims.append(len(loc[4]) + loc[5][0]-l[3][1])
            trims.append((left_trim, right_trim))
        return trims
    
    def get_inserts_and_trims(self, features):
        """
        Returns list of insert size, and us/ds trims with the following format:
          [(insert, us-trim, ds-trim), ...]
        """
        return_list = []
        for loc in self.locations:
            contig, start, end, strand = features[loc[0]]
            mAT = loc[4]
            left_trim = start - loc[3][1]
            right_trim = loc[5][0] - end
            return_list.append((mAT, left_trim, right_trim))
        return return_list

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
    
    def __init__(self, feature, contig, orientation, start, end, upstream, downstream, sequence, side, spacer, pam, motif, parsed_motif):
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
        aligned_parsed_motif = args.parsed_motifs[args.motifs.index(aligned_motif)]
        if ((aligned_target != None) and (aligned_pam != None)):
            a = Alignment(
                aligned_sequence,
                aligned_target,
                aligned_pam,
                aligned_motif,
                aligned_parsed_motif,
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
    
    def get_contigs(self):
        return sorted(set(x[1] for x in self.locations))
    
    def get_location_features(self):
        return sorted(set(x[0] for x in self.locations))
    
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
                        a_features = self.get_features(features, a.contig, a.start, a.end) # contig, start, end
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
                                cls(feature, feature_contig, orientation, real_start, real_end, t_upstream, t_downstream, filt_seq, side, filt_spacer, filt_pam, args.motifs[i], args.parsed_motifs[i])
                                # Maybe add to ReversionDonor here

class ReversionTarget(Target):
    prefix = 'reTarget'
    # Non-redundant dicts
    sequences = {} # key = nucleotide sequence, value = ReversionTarget object
    indices = {} # key = reTarget-234, value = ReversionTarget object
    
    @classmethod
    def get_targets(cls):
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
    parser.add_argument("-v", "--version", action='version', version='{__program__} {__version__} (revision {__revision__})'.format(**globals()))
    
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
        default=["N{20}>NGG"],
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
    parser.add_argument("--dDNA_gDNA_ratio", metavar="N", type=int, default=1000,
        help="Ratio of donor DNA to genomic DNA for calculating off-target scores")
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
#    parser.add_argument("--min_donor_distance", metavar="N", type=int, default=36,
#        help="The minimum distance in bp a difference can exist from the edge of donor DNA") # homology with genome
    parser.add_argument("--max_consecutive_ts", metavar="N", type=int, default=4,
        help="The maximum number of Ts allowed in generated gRNA sequences")
    parser.add_argument("--max_number_sequences_reported", metavar="N", type=int, default=5,
        help="The maximum number of sequences to report for each step")
    parser.add_argument("--min_weight_reported", metavar="N", type=float, default=0.01,
        help="Only gRNA-dDNA pairs with at least this much weight will be reported.")
    # program currently will only search 'both' strands
    #parser.add_argument("--strands", type=str, choices=["+", "-", "both"], default="both",
    #    help="Strands to search for gRNAs")
    
    # Add command line arguments for the additional hard constraints:
    #  Only report potential targets that have no off targets with mismatches within 8, 12, N nt from 3' end
    parser.add_argument("--processors", metavar="N", type=int, default=(os.cpu_count() or 1),
        help="Number of processors to use when performing pairwise sequence alignments")
    
    #available_aligners = ['addtag', 'blast+', 'blat', 'bowtie', 'bowtie2', 'bwa', 'cas-offinder']
    available_aligners = list(map(lambda x: x.name, aligners.aligners))
    parser.add_argument("--aligner", type=str, choices=available_aligners, default='bowtie2',
        help="Program to calculate pairwise alignments. Please note that the 'addtag' internal aligner is very slow.")
    # Other aligners to consider: 'rmap', 'maq', 'shrimp2', 'soap2', 'star', 'rhat', 'mrsfast', 'stampy'
    
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
    
    # Parse the motifs
    # Compile regex for motifs
    # Note that any logging these functions perform will not be written to disk,
    # as no 'handler' has been specified.
    args.parsed_motifs = []
    args.compiled_motifs = []
    for motif in args.motifs:
        spacers, pams, side = parse_motif(motif) # Parse the motif
        args.parsed_motifs.append((spacers, pams, side)) # Add to args
        args.compiled_motifs.append(nucleotides.compile_motif_regex(spacers, pams, side, anchored=False)) # Add to args
    
    # Add 'selected_aligner' to hold the actual aligner object
    for a in aligners.aligners:
        if a.name == args.aligner:
            args.selected_aligner = a
            break
    
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
    
    # Discard potential gRNAs that have mismatches with their target site
    #if (args.case == "invariant-lower"):
    #    pass
    #elif (args.case == "invariant-upper"):
    #    pass
    
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
    
#    # Apply all prefilters
#    this = (sequence, target, pam, upstream, downstream)
#    prefilter = []
#    for C in algorithms.single_algorithms:
#        c_score = C.calculate(this)
#        if C.prefilter:
#            if (C.minimum <= c_score <= C.maximum):
#                prefilter.append(True)
#            else:
#                prefilter.append(False)
#    # parent = (sequence, target, pam, upstream, downstream)
#    #for C in algorithms.paired_algorithms:
#    #    c_score = C.calculate(parent, this)
#    #    if C.prefilter:
#    #        if (C.minimum <= c_score <= C.maximum):
#    #            prefilter.append(True)
#    #        else:
#    #            prefilter.append(False)
#    
#    # If alignment meets all prefilter criteria, then set it as True
#    prefilter = all(prefilter) # all([]) returns True
    
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

def get_exTarget_homologs(features, homologs):
    """Get ExcisionTarget objects for each homologous feature group"""
    ext_dict = {}
    if homologs:
        groups = set()
        for f in features:
            groups.add(tuple(sorted(homologs[f])))
        for g in groups:
            ext_dict[g] = set()
            for name, obj in ExcisionTarget.indices.items():
                obj_features = set(utils.flatten(x[0].split(',') for x in obj.locations))
                if (len(obj_features.intersection(g)) == len(g)):
                    ext_dict[g].add(obj)
    return ext_dict

def get_exTarget_allele_specific(features):
    """Gets allele-specific ExcisionTarget objects"""
    ext_dict = {}
    for f in features:
        ext_dict[f] = set()
        for name, obj in ExcisionTarget.indices.items():
            obj_features = set(utils.flatten(x[0].split(',') for x in obj.locations))
            if ((len(obj_features) == 1) and (f in obj_features)):
                ext_dict[f].add(obj) # Store with the feature as key
    return ext_dict

def get_reTarget_homologs(features, homologs):
    """Get ReversionTarget objects for each homologous feature group"""
    ret_dict2 = {}
    
    if homologs:
        groups = set()
        for f in features:
            #contig, start, end, strand = features[feature]
            # Convert to tuple
            groups.add(tuple(sorted(homologs[f])))
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

def get_reTarget_allele(features):
    """Gets all ReversionTarget objects for each feature"""
    ret_dict = {}
    for f in features:
        ret_dict[f] = set()
        for name, obj in ReversionTarget.indices.items():
            obj_features = set(utils.flatten(x[0].split(',') for x in obj.locations))
            if (f in obj_features):
                ret_dict[f].add(obj) # Store with the feature as key
    return ret_dict

def get_reTarget_allele_specific(features):
    """Gets allele-specific ReversionTarget objects"""
    ret_dict = {}
    for f in features:
        ret_dict[f] = set()
        for name, obj in ReversionTarget.indices.items():
            obj_features = set(utils.flatten(x[0].split(',') for x in obj.locations))
            if ((len(obj_features) == 1) and (f in obj_features)):
                ret_dict[f].add(obj) # Store with the feature as key
    return ret_dict

def rank_donors(donor_list):
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

def rank_targets(target_list):
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

def get_best_table(args, features, homologs, feature2gene):
    """
    Identify and print the best spacers and dDNAs for each feature, given
    each mAT insert size, and us/ds trim length.
    Thus, the user can see the best spacer for any given combination of these.
    """
    #header = ['gene', 'features', 'insert', 'mAT', 'translations', '(us, ds) trim', 'weight', 'OT:Hsu-Zhang', 'OT:CFD', 'Azimuth', 'reTarget name', 'reTarget sequence', 'ExDonors']
    header = ['gene', 'features', 'us-trim:mAT:ds-trim', 'translations', 'weight', 'OT:Hsu-Zhang', 'OT:CFD', 'Azimuth', 'reTarget name', 'reTarget sequence', 'exDonors']
    print('\t'.join(header))
    
    # Get best ReversionTargets by calculating their weights, and also getting
    # the ExcisionDonors that correspond to the ReversionTarget
    ret_dict2 = get_reTarget_homologs(features, homologs)
    for feature_homologs in sorted(ret_dict2):
        # Get these two things which are common to all the top hits
        gene = feature2gene[feature_homologs[0]] # Get the gene name
        csfeatures = ','.join(feature_homologs) # Add the features as comma-separated list
        
        outputs = {}
        
        # Print the top N for each insert size and trim
        for weight, obj in rank_targets(ret_dict2[feature_homologs]):
            othz = round(obj.off_targets['Hsu-Zhang'], 2)
            otcfd = round(obj.off_targets['CFD'], 2)
            azimuth = round(obj.score['Azimuth'], 2)
            
            # Get the ExcisionDonor objects for this ki-spacer, and weigh them
            rds = rank_donors(obj.get_donors())
            # filter out all but the top-weighted ones
            rds = [x for x in rds if (x[0] == rds[0][0])]
            
            exdonors = ','.join(map(lambda x: x[1].name, rds))
            
            key0 = sorted(set(utils.flatten(map(lambda x: x[1].get_inserts_and_trims(features), rds))))
            key0s = ','.join('{}:{}:{}'.format(x[1], x[0], x[2]) for x in key0)
            key1 = sorted(set([(len(x[0]), x[1], x[2]) for x in key0])) # replace mAT with length
            key1s = ','.join(map(str, key1))
            translations = None
            
            sline = [gene, csfeatures, key0s, translations, weight, othz, otcfd, azimuth, obj.name, obj.spacer+'|'+obj.pam, exdonors]
            
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
    ret_dict = get_reTarget_allele_specific(features)
    for feature in sorted(ret_dict):
        gene = feature2gene[feature] # Get the gene name
        outputs = {}
        
        for weight, obj in rank_targets(ret_dict[feature]):
            othz = round(obj.off_targets['Hsu-Zhang'], 2)
            otcfd = round(obj.off_targets['CFD'], 2)
            azimuth = round(obj.score['Azimuth'], 2)
            
            rds = rank_donors(obj.get_donors())
            rds = [x for x in rds if (x[0] == rds[0][0])]
            
            exdonors = ','.join(map(lambda x: x[1].name, rds))
            
            key0 = sorted(set(utils.flatten(map(lambda x: x[1].get_inserts_and_trims(features), rds))))
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
    
    header = ['gene', 'features', 'weight', 'OT:Hsu-Zhang', 'OT:CFD', 'Azimuth', 'exTarget name', 'exTarget sequence', 'reDonors']
    print('\t'.join(header))
    
    # Print the best homozygous ExcisionTargets for each feature set
    ext_dict2 = get_exTarget_homologs(features, homologs)
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
        for weight, obj in rank_targets(ext_dict2[feature_homologs]):
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
    ext_dict = get_exTarget_allele_specific(features)
    for feature in sorted(ext_dict):
        gene = feature2gene[feature] # Get the gene name
        
        red_list = set()
        for name, obj in ReversionDonor.indices.items():
            if feature in obj.get_location_features():
                red_list.add(obj)
        red_list = sorted(red_list, key=lambda x: int(x.name.split('-')[1]))
        
        outputs = []
        
        for weight, obj in rank_targets(ext_dict[feature]):
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
    

def get_best(args, features, homologs):
    """Function that returns the best spacers and dDNAs for each feature"""
    
    display_num = 5
    
    # Print best ReversionTargets calculated and their corresponding ExcisionDonors
    ret_dict2 = get_reTarget_homologs(features, homologs)
    for k in ret_dict2:
        logging.info(str(k) + ' ' + str(len(ret_dict2[k])))
        # Print the top 5
        for rank, obj in rank_targets(ret_dict2[k])[:display_num]:
            logging.info(' ' + str(rank) + ' ' + str(obj))
            # Get the ExcisionDonor objects for this ki-spacer, and rank them
            rds = rank_donors(obj.get_donors())
            # filter out all but the top-ranked ones
            rds = [x for x in rds if (x[0] == rds[0][0])]
            for gap, exd_obj in rds:
                logging.info('   ' + str(gap) + ' ' + str(exd_obj.get_trims(features)) + ' ' + str(exd_obj))
    
    # Print best ExcisionTargets (not necessarily homozygous) for each feature
    # and the ReversionDonor
    for feature in sorted(features):
        logging.info(feature)
        et_list = []
        for name, obj in ExcisionTarget.indices.items():
            if feature in obj.get_location_features():
                et_list.append(obj)
        for rank, obj in rank_targets(et_list)[:display_num]:
            logging.info('  ' + str(rank) + ' ' + str(obj))
        
        red_list = []
        for name, obj in ReversionDonor.indices.items():
            if feature in obj.get_location_features():
                red_list.append(obj)
        for obj in red_list:
            logging.info('  ' + str(obj))

def old_get_best(args, features, contigs):
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

def main():
    """Function to run complete AddTag analysis"""
    
    # Obtain command line arguments and parse them
    args = parse_arguments()
    
    # Get timestamp for analysis beginning
    start = time.time()
    
    # Create the project directory if it doesn't exist
    os.makedirs(args.folder, exist_ok=True)
    
    # Create the logger
    logging.basicConfig(filename=os.path.join(args.folder, 'log.txt'), level=logging.INFO, format='%(message)s') # format='%(levelname)s %(asctime)s: %(message)s'
    
    # Echo the command line parameters to STDOUT and the log
    print(args, flush=True)
    logging.info(args)
    
    # Load the FASTA file specified on the command line
    contigs = utils.load_fasta_file(args.fasta)
    
    # Open and parse the GFF file specified on the command line
    features = utils.load_gff_file(args.gff, args.features, args.tag)
    
    # Make index of homologs
    if args.feature_homologs:
        homologs, feature2gene = utils.load_homologs(args.feature_homologs)
    else:
        homologs, feature2gene = None, None
    
    # Merge features?
    #features = merge_features(features)
    
    # Search features within contigs for targets that match the motifs
    ExcisionTarget.get_targets(args, contigs, features)
    
    # Write the query list to FASTA
    ex_query_file = ExcisionTarget.generate_query_fasta(os.path.join(args.folder, 'excision-query.fasta'))
    
    # Generate excision dDNAs and their associated reversion gRNA spacers
    ExcisionDonor.generate_donors(args, features, contigs)
    ReversionTarget.get_targets()
    ex_dDNA_file = ExcisionDonor.generate_fasta(os.path.join(args.folder, 'excision-dDNAs.fasta'))
    re_query_file = ReversionTarget.generate_query_fasta(os.path.join(args.folder, 'reversion-query.fasta'))
    
    # Generate reversion dDNAs and write them to FASTA
    ReversionDonor.generate_donors(args, features, contigs)
    re_dDNA_file = ReversionDonor.generate_fasta(os.path.join(args.folder, 'reversion-dDNAs.fasta'))
    
    # Index args.fasta for alignment
    #index_file = index_reference(args)
    genome_index_file = args.selected_aligner.index(args.fasta, os.path.basename(args.fasta), args.folder, args.processors)
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
    ExcisionTarget.load_sam(exq2gDNA_sam_file, args, contigs)
    ExcisionTarget.load_sam(exq2exdDNA_sam_file, args, ExcisionDonor.get_contig_dict())
    
    # Calculate off-target/guide scores for each algorithm
    logging.info("ExcisionTarget after SAM parsing and off-target scoring")
    for et_seq, et_obj in ExcisionTarget.sequences.items():
        et_obj.score_off_targets(args, homologs, features)
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
    ReversionTarget.load_sam(req2gDNA_sam_file, args, contigs)
    ReversionTarget.load_sam(req2exdDNA_sam_file, args, ExcisionDonor.get_contig_dict())
    ReversionTarget.load_sam(req2redDNA_sam_file, args, ReversionDonor.get_contig_dict())
    
    # Calculate off-target/guide scores for each algorithm
    logging.info("ReversionTarget after SAM parsing and off-target scoring")
    for re_seq, re_obj in ReversionTarget.sequences.items():
        re_obj.score_off_targets(args, homologs, features)
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
    get_best(args, features, homologs)
    get_best_table(args, features, homologs, feature2gene)
    
    # Print time taken for program to complete
    logging.info('{} finished'.format(__program__))
    logging.info('Runtime: {}s'.format(time.time()-start))
