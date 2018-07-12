#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/oligos/addtagprimer.py

# Import standard packages
import sys
import math
import os
import subprocess

# import non-standard package
import regex

# Code to analyze PCR primers
# IDT uses mfold_util 4.5 for hairpins (UNAFold)

# JalView uses RNAAliFold

# Primer design possibilities for wild-type, mintag, addtag, mutant/transgene
#     ─genome┐┌─us_homology─┐┌─feature─┐┌─ds_homology─┐┌genome─
# Fo  ──> - - - - - - - - - - - - - - - - - - - - - - - - - <──  Ro
# Fo  ──> - - - - - - - - - - <──  Rf
#                             Ff  ──> - - - - - - - - - - - <──  Ro
#                         Ff  ──>- -<──  Rf
# Sample FASTA headers for these
# >oligo-1 pair=oligo-2 location=feature strand=+ feature=NNNNNNN amplicon_size=N
# >oligo-2 pair=oligo-1 location=feature strand=- feature=NNNNNNN amplicon_size=N
#
# >oligo-3 pair=oligo-4 location=upstream strand=+ feature=NNNNNNN amplicaon_size=N
# >oligo-4 pair=oligo-3 location=feature strand=- feature=NNNNNNN amplicaon_size=N
#
# >oligo-5 pair=oligo-6 location=feature strand=+ feature=NNNNNNN amplicaon_size=N
# >oligo-6 pair=oligo-5 location=downstream strand=- feature=NNNNNNN amplicaon_size=N



# Primer design possibilities for flanktag
#     ─genome┐┌─us_homology─┐┌─uptag─┐┌─bartag─┐┌─dntag─┐┌─ds_homology─┐┌genome─
# Fo  ──> - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -<──  Ro
# Fo  ──> - - - - - - - - - - <──  Ru
# Fo  ──> - - - - - - - - - - - - - - - <──  Rb
# Fo  ──> - - - - - - - - - - - - - - - - - - - - - <──  Rd
#                          Fu  ──>- - - - - - - - - - - - - - - - - - - - - -<──  Ro
#                                      Fb  ──>- - - - - - - - - - - - - - - -<──  Ro
#                                                Fd  ──>- - - - - - - - - - -<──  Ro
#                           Fu  ──> - - - - - - - -<──  Rd
# Sample FASTA headers for these
# >oligo-7 pair=oligo-8 location=upstream strand=+ feature=NNNNNNN amplicaon_size=N
# >oligo-8 pair=oligo-7 location=uptag strand=- feature=NNNNNNN amplicaon_size=N

class Sequence(object):
    def __init__(self, sequence, concentration, molecule='DNA'):
        self.sequence = sequence # ACGT...
        self.molecule = molecule # 'DNA' or 'RNA'
        self.concentration = concentration # concentration in molarity

class Geometry(object):
    def __init__(self, sequences):
        self.sequences = sequences # [Sequence(), ...]
        self.delta_G = self.get_delta_G()
        self.delta_H = self.get_delta_H()
        self.delta_S = self.get_delta_S()
        self.probability = None
        self.melting_temperature = self.get_melting_temperature()

    def get_delta_G(self, sodium=0.05, magnesium=0.0):
        pass
    def get_delta_H(self):
        pass
    def get_delta_S(self):
        pass
    def get_melting_temperature(self):
        pass

def get_geometries(seqs, sodium=0.05, magnesium=0.0):
    pass

def make_primer():
    # primer length should be 17-28 nt long
    # 3' GC clamp
    # %GC should be between 40-60
    # Tm should be between 50-70
    # min delta-G should be -3
    # should avoid runs of 3-or-more Cs and Gs in 3' end
    # 3' ends of primers should NOT be complementary
    #    (as even a little dimerization will inhibit target annealing)
    
    pass

def get_best_pair(primer, seq):
    """
    Searches the sequence 'seq' for the best primers to pair with 'primer'
    """
    # amplicon size: min=200, opt=Npone, max=1000
    # number results: 5
    # Exclude case-masked characters
    pairs = []
    return pairs

def get_primer_pair(seq):
    """
    seqrch input sequence for the best primer pairs
    """
    # amplicon size: min=200, opt=None, max=1000
    # number results: 5
    
    # Exclude case-masked characters
    pairs = []
    return pairs

def get_primers(seq):
    """
    Search input sequence for list of good primers
    """
    # number results: 5
    # Exclude case-masked characters
    primers = []
    return primers

def score_primer_pair():
    # Check primer pair
    # make sure difference between primer pair melting temps should be minimized (<5 C)
    return 0.0

def score_primer():
    """
    Takes input configuration, and outputs a score stating how good the primer should be.
    Higher scores are better.
    Use this score to rank primers.
    """
    # primer considerations
    
    # Primer length: typically 18-22 nt
    # deta-G: You want to avoid structures with delta-G < -5 kcal/mol
    # GC content (should be between 40-60%)
    # successively penalize long stretches of AAAA, CCCC, GGGG, or TTTT
    
    # check melting temperature
    #  typically 52-68 C
    #  better: 55-65 C
    # penalize score if temp is too low or too high
    
    # %GC         min=35, opt=50, max=65
    # Tm          min=55, opt=62, max=65
    # size        min=18, opt=22, max=26
    # 3' GC clamp min=0
    
    return 0.0

def get_length(seq):
    """Returns length of the oligonucleotide"""
    return len(seq)

def get_repeats(seq):
    # Primers should not have palindromes
    # primers should not have inverted repeat sequences.
    # Such sequences cause primers to either dimerize or form hairpins that interfere with proper annealing to the template.
    return None

def get_hairpins(seq):
    """
    Identify any potential hairpins and some statistics about each
    Returns
     alignment
     delta-G
     self-melting temperature?
    """
    # if there is a hairpin, we can expect a self-dimer also
    return None

def get_self_dimers(seq):
    """
    Identify any potential self-dimers, and return statistics about each
    Returns
     maximum delta-G
     list of
      alignment
      delta-G
      paired melting temperature?
      number base pairs
    """
    # If there are self-dimers, then we should also expect a hairpin
    # Also called homo-dimers
    # The delta-G is calculated by taking into account the longest stretch of
    # complementary bases. These pairs of complementary bases are represented
    # by a solid line. Dotted lines represent additional complementary bases
    # for that dimer structure, but their presence does not impact calculated
    # delta G values. Actual delta-G values may vary based on presence of
    # additional complementary bases. The Maximum Delta G value refers to
    # the free energy of the oligo sequence binding to its perfect complement.
    
    # Example return
    #  Input Sequence:   5'- ACGTATCGAGCAACTAGCG -3'
    #  Maximum Delta G:  -36.38 kcal/mole
    #
    #  Delta G:    -6.76 kcal/mole
    #  Base Pairs:  4
    #  Alignment:
    #   5'      ACGTATCGAGCAACTAGCG
    #                ||||          
    #   3' GCGATCAACGAGCTATGCA
    #
    #
    #  Delta G:    -4.16 kcal/mole
    #  Base Pairs:  4
    #  5' ACGTATCGAGCAACTAGCG
    #                  ||||             
    #  3'            GCGATCAACGAGCTATGCA
    #
    #
    #  Delta G:  -3.61 kcal/mole
    #  Base Pairs:  2
    #  5' ACGTATCGAGCAACTAGCG
    #      ||      ::      :: 
    #  3'  GCGATCAACGAGCTATGCA
    
    # alignment ideas...
    #  strong   |
    #  moderate !
    #  weak     :
    #  fragile  .
    #  top      ^  '
    #  bottom   v  ,
    #  cross    X
    #   
    return None

def get_dimers(seq1, seq2):
    """
    Identify potential dimers between seq1 and seq2
    Returns
     alignment
     delta-G
     paired melting temperature?
    """
    
    # primers should NOT be palindromes?
    #  if (seq1 == rc(seq2)):
    #      # these are palindromes
    #  ACGT = rc(ACGT)
    #  
    # primers should NOT be complements to one another
    #  if (seq1 == c(seq2)):
    #      # These are complements (no reversing)
    #     ACGT & TGCA will anneal to each other
    #   

    
    return None

def get_list_dimers(seq, oligos):
    """
    Identify potential dimers between seq and each other oligo,
    and return statistics about each
    Returns
     alignment
     delta-G
     paired melting temperature?
    """
    dimers = []
    for o in oligos:
        dimers.append(get_dimers(seq, o))
    return dimers

def get_multiplex_dimers(oligos):
    """
    Identify potential dimers between all oligos
    Returns 
     alignment
     delta-G
     paired meltaing temparature?
    """
    dimers = []
    for o1 in oligos:
        for o2 in oligos:
            dimers.append(get_dimers(o1, o2))
    return dimers

def get_gc(seq):
    """
    Returns GC content
    """
    return None

def get_runs(seq):
    """
    Count the number of times a base repeats:
     eg CCCC, GGGG, AAAA, TTTT
    """
    # Avoid runs of over 3 nucleotides (i.e., CCCC). Long stretches of G, in particular, give problems.
    # long poly-G or poly-C should be avoided (as it can promote non-specific annealing)
    # long Poly-A or poly-T should be avoided (as it can "breathe" and open the primer-template complex?)
    # polypyrimidine(T, C) and polypurine(A, G) would be avoided. (as it leads to an odd shape of the double helix.)
    return None

def get_melting_temperature(seq):
    """
    Calculate melring temperature of oligonucleotide
    """
    # Annealing temp should typically be 5 C less than melting temperature
    
    # Thermodynamic parameters
    # Algorithm: SantaLucia 1999/1998 (Recommended)
    #            Modified Breslauer 1986 (Phusion, Phire, DyNAzyme))
    #            Marmur
    #            Wallace
    # Defaults: Benchling  Primer3  IDT
    # DNA:      250 nM     50 nM    250 nM
    # MG++:     0 mM       1.5 mM   0 mM   # divalent cation/salt
    # Na+/K+:   50 mM      50 mM    50 mM  # monovalent cation/salt
    # dNTP:     0 mM       0.6 mM   0 mM   # nucleotide triphosphate
    # 
    #
    # Marmur and Wallace formulae Tm estimation only take into account the
    # number of GC and AT nucleotides. The position of the nucleotides in
    # the primer is not considered (primer 1-4: same Tm; primer 5-6: longer,
    # higher Tm).
    
    # Both Breslauer and SantaLucia nearest-neighbour thermodynamics factor in
    # the nucleotide environment. GC rich islands lead to higher Tm estimates
    # with the upwards trend being much more pronounced when the Breslauer
    # calculations are used (primer 1, 2).
    
    return 0.0

def get_marmur(seq):
    # Melting temperature (Tm) estimation 
    # Marmur formula: Tm = 4 x GC + 2 x AT
    # not recommended for more than 13nt; assumes 50mM monovalent cations
    # Marmur J and Doty P (1962) J Mol Biol 5:109-118; PMID 14470099
    
    # For sequences less than 14 nucleotides the formula is:
    #  Tm= (wA+xT) * 2 + (yG+zC) * 4
    # where w,x,y,z are the number of the bases A,T,G,C in the sequence, respectively.
    
    tm = 0
    for nt in seq:
        if (nt in ['A', 'T', 'a', 't']): # Does not handle IUPAC
            tm += 2
        elif (nt in ['G', 'C', 'g', 'c']): # Does not handle IUPAC
            tm += 4
    return tm

def get_wallace(seq):
    # Melting temperature (Tm) estimation 
    # Wallace formula: Tm = 64.9 +41*(yG+zC-16.4)/(wA+xT+yG+zC)
    # Wallace RB et al. (1979) Nucleic Acids Res 6:3543-3557, PMID 158748
    # online tool (http://www.basic.northwestern.edu/biotools/oligocalc.html)
    # using Wallace formula for oligos >13
    
    # where w,x,y,z are the number of the bases A,T,G,C in the sequence, respectively.
    # For sequences longer than 13 nucleotides, the equation used is
    #  Tm= 64.9 +41*(yG+zC-16.4)/(wA+xT+yG+zC)
    # When degenerated nucleotides are included in the primer sequence
    # (Y,R,W,S,K,M,D,V,H,B or N), those nucleotides will be internally
    # substituted prior to minimum and maximum Tm calculation.
    #
    # Example:
    #  Primer sequence:                            CTCTRYCTWSCTCTCT
    #  Sequence for minimum Tm calculation:        CTCTATCTAGCTCTCT
    #  Sequence for maximum Tm calculation:        CTCTGCCTAGCTCTCT
    #                                                  ^^  ^^
    
    A = seq.count('A')
    C = seq.count('C')
    G = seq.count('G')
    T = seq.count('T')
    tm = 64.9 + 41*(G + C - 16.4)/(A + T + G + C) # Does not account for IUPAC ambiguities
    return tm

def get_biophp(seq):
    # supposedly...
    # from table at http://www.ncbi.nlm.nih.gov/pmc/articles/PMC19045/table/T2/ (SantaLucia, 1998)
    # BioPython claims these values are from Allawi and SantaLucia (1997), Biochemistry 36: 10581-10594
    # enthalpy values (dH)
    array_h = {
        'AA': -7.9, 'TT': -7.9,
        'AC': -8.4, 'GT': -8.4,
        'AG': -7.8, 'CT': -7.8,
        'CA': -8.5, 'TG': -8.5,
        'CC': -8.0, 'GG': -8.0,
        'GA': -8.2, 'TC': -8.2,
        'AT': -7.2,
        'CG':-10.6,
        'GC': -9.8,
        'TA': -7.2,
    }
    # entropy values (dS)
    array_s = {
        'AA': -22.2, 'TT': -22.2,
        'AC': -22.4, 'GT': -22.4,
        'AG': -21.0, 'CT': -21.0,
        'CA': -22.7, 'TG': -22.7,
        'CC': -19.9, 'GG': -19.9,
        'GA': -22.2, 'TC': -22.2,
        'AT': -20.4,
        'CG': -27.2,
        'GC': -24.4,
        'TA': -21.3,
    }
    tm = 0
    return tm

def get_idt(seq):
    # see
    #  https://www.idtdna.com/calc/Analyzer/Home/Definitions
    
    tm = 0.0
    return tm

def get_promega(seq):
    # Promega Formula
    # retrieved from http://www.promega.com/a/apps/biomath/?calc=tm
    # on 07-28-2017
    # 
    # The most sophisticated Tm calculations take into account the exact sequence and base stacking parameters, not just the base composition(1,2,3).
    # 
    # The equation used is:
    #  Tm = dH (kcal / ((degrees C)*Mol)) / (dS + R ln ([primer] / 2)) - 273.15 (degrees C)
    # 
    # dH = the enthalpy of base stacking interactions adjusted for helix initiation factors (3,4).
    # 
    # dS = the entropy of base stacking adjusted for helix initiation factors (3,4) and for the contributions of salts to the entropy of the system (3).
    # 
    # R = the universal gas constant: 1.987 Cal / (degrees C) * Mol 
    # 
    # Most melting temperature calculations do not take into account the effects of magnesium on helix stability. Therefore, most empirical guidelines used to design experiments will not apply when the magnesium effects are included. We have included the option to consider magnesium in the equation if it is desirable but have not included it in the default setting. Including magnesium will generally raise the theoretical melting temperature by about 5-8 (degrees C) for oligonucleotides in a 1.5 mM Mg++ solution (5,6).
    # 
    # 1. Rychlik, W. and Rhoads, R.E. (1989) Nucl. Acids Res. 17, 8543.
    # 2. Borer P.N. et al. (1974) J. Mol. Biol. 86, 843.
    # 3. SantaLucia, J. (1998) Proc. Nat. Acad. Sci. USA 95, 1460.
    # 4. Allawi, H.T. and SantaLucia, J. Jr. (1997) Biochemistry 36, 10581.
    # 5. von Ahsen N. et al. (1999) Clin. Chem. 45, 2094.
    # 6. Nakano S. et al. (1999) Nucl. Acids Res. 27, 2957.
    
    dh = 1 # in degrees C
    ds = 1 # in degrees C
    r = 1
    primer = 1 # Concentration of the primer
    tm = dh/(ds+r*math.log(primer/2)) - 273.15 # tm is in (degrees C)
    
    return tm

def get_meh1(seq):
    # Units in nmol/OD260
    return None

def get_meh2(seq):
    # units in ug/OD260:
    return None

def get_extinction_coefficient(seq):
    # Units in L/(mol * cm)
    return None

def get_molecular_weight(seq):
    """
    Returns the molecular weight of the oligonucleotide
    """
    
    # Molecular Weight (Anhydrous)
    # Molecular weight (MW) is the sum of the atomic masses of the constituent
    # atoms for 1 mole of oligonucleotide.
    #
    # The anhydrous molecular weight represents the pure oligo free of any of
    # the counter ions or water molecules that are normally weakly bound to
    # an oligo after synthesis. This calculation gives with the molecular
    # weight measured by mass spectroscopy.
    # 
    # Molecular weight of an oligomer is a sum of the weights of individual
    # bases and chemical modifications. Oligos are typically synthesized
    # without 5'-phosphate group, which must be subtracted.
    # 
    # Anhydrous MW = sum(Individual Base MW) + sum(Individual Mod MW) - PO2H + H2
    #  where PO2H = 63.980 and H2 = 2.016
    # 
    # Molecular weights of DNA bases:
    #  MW(dA) = 313.209
    #  MW(dC) = 289.184
    #  MW(dG) = 329.208
    #  MW(dT) = 304.196
    #  MW(dU) = 290.169
    #  MW(dI) = 314.194
    # 
    # RNA bases:
    #  The molecular weight of an RNA nucleotide is the weight of a DNA
    #  nucleotide + 15.999, accounting for the additional oxygen atom present
    #  (Example: rA is dA (313.209) + 15.999 = 329.208). When determining the
    #  weight of uracil (rU) start with dU and not thymine (dT).
    # 
    # 2'-O-Methyl bases:
    #  The molecular weight of an 2'-O-methyl RNA nucleotide is the weight of
    #  a DNA nucleotide + 30.026, accounting for the additional methoxy group
    #  (-OCH3) present (For example, mA is dA (313.209) + 30.026 = 343.235).
    #  When determining the weight of uracil (mU) start with dU and not
    #  thymine (dT).
    # 
    # LNA bases:
    #  The molecular weight of an LNA nucleotide is the weight of a DNA
    #  nucleotide + 28.011, accounting for the bridging oxygen and carbon from
    #  the 2' carbon to the 4' carbon (-OCH2-) present
    #  (For example, +A is dA (313.209) + 28.011 = 341.220). The exception is
    #  LNA C which contains an additional methyl group off of the 5-carbon
    #  (14.026).
    # 
    # Phosphorothioated bases:
    #  The molecular weight of a phosphorothioate nucleotide is the weight of a
    #  nucleotide (DNA, RNA, 2'-O-methyl, LNA) + 16.061, accounting for the
    #  substitution of one sulfur atom for a non-bridging oxygen atom in the
    #  phosphodiester backbone
    #  (For example, A* is A (313.209) + 16.061 = 329.270). Phosphorothioate
    #  modification refers to substitutions affecting the internucleoside
    #  linkages and does not involve the free 3'- or 5'- ends. Thus a 20-mer
    #  phosphorothioate oligonucleotide has 19 phosphorothioate linkages.
    # 
    # See mixed base definitions and modifications for their effects on molecular weight.
    
    
    # units in g/mol
    return 0.0

def get_strong_tail(seq):
    # 3' end of the primer (0-5 nt from the end)
    # The 3' end of the primer should be an exact match to the template DNA.
    # GC clamp:
    #  should not contain more than 3 C/Gs
    #  should not end with T/A (i.e. it should end with G, C, GC, CG) (what about CC or GG?)
   return None

def get_tail_delta_g(seq):
    """Calculate delta-G for the 3' end of the sequence"""
    #3' end stability: the maximum delta-G value (less negative (i.e. the maximum) delta-G will result in less false-priming)
    return 0.0

def get_delta_g(seq):
    """calculate delta-G for entire sequence"""
    return 0.0

def check_head_adapter(seq):
    # 5' end
    # 5' Y-adapter: 5' tails can be readily added to primers without impacting primer annealing.
    return None

def test():
    """Code to test the functions and classes"""
    pass

if (__name__ == "__main__"):
    test()
