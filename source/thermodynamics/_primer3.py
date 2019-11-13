#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/thermodynamics/_primer3.py

# Import standard packages
import sys
import os
import logging
import math
import importlib.util

# Import non-standard packages
import regex

logger = logging.getLogger(__name__)

# import included AddTag-specific modules
if (__name__ == "__main__"):
    from oligo import Oligo, Structure
else:
    from .oligo import Oligo, Structure

# Treat modules in PACKAGE_PARENT as in working directory
if (__name__ == "__main__"):
    # Relative path for package to import
    PACKAGE_PARENT = '..'
    # Obtain path of currently-running file
    SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
    # Convert to absolute path, and add to the PYTHONPATH
    sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))
    
    from nucleotides import rc
else:
    from ..nucleotides import rc

class Primer3(Oligo):
    logger = logger.getChild(__qualname__)
    
    def __init__(self):
        super().__init__(
            name="Primer3",
            authors=['Untergasser, Andreas', 'Cutcutache, Ioana', 'Koressaar, Triinu', 'Ye, Jian', 'Faircloth, Brant C.', 'Remm, Maido', 'Rozen, Steven G.'],
            title='Primer3--new capabilities and interfaces',
            journal='Nucleic Acids Research',
            issuing='40(15):e115',
            year=2012,
            doi='https://doi.org/10.1093/nar/gks596'
            #citation="Untergasser, et al. Primer3--new capabilities and interfaces. Nucleic Acids Research 40(15): e115 (2012)"
        )
        spec = importlib.util.find_spec('primer3')
        if spec:
            self.available = True
    
    @classmethod
    def find_structures(cls, folder, seq1, seq2=None, sodium=0.05, magnesium=0.0, temperature=25, concentration=0.00000025, **kwargs):
        """
        Should return the list of 'Structure' objects with delta-G, deltaH, deltaS, and Tm values.
        
        Accepts 1 or 2 input sequences. Automatically runs either:
         * Hairpin     (1 input sequence: A=seq1, UNAFold run on A)
         * Homodimer   (2 identical input sequences: A=seq1=seq2, UNAFold run on A & A)
         * Heterodimer (2 input sequences: A=seq1 B=seq2, UNAFold run on A & B)
        """
        import primer3

        mv_conc = sodium*1000 # 50.0 # in mM
        dv_conc = magnesium*1000 # 0.0 # in mM
        dntp_conc = 0.6 # in mM
        dna_conc = concentration*1000*1000*1000 # 250.0 # in nM
        #temperature = 25 # keep as-is
        
        if (seq1 == seq2): # Homodimer calculation
            t = primer3.calcHomodimer(seq1, mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc, dna_conc=dna_conc, temp_c=temperature)
        elif (seq2 == None): # Hairpin calculation
            t = primer3.calcHairpin(seq1, mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc, dna_conc=dna_conc, temp_c=temperature)
        else: # Heterodimer calculation, Tm calculation [seq1, rc(seq1)]
            t = primer3.calcHeterodimer(seq1, seq2, mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc, dna_conc=dna_conc, temp_c=temperature)
        
        if t.structure_found:
            s = Structure(seq1, seq2, t.dg/1000, t.dh/1000, t.ds, t.tm, sodium, magnesium, temperature, concentration)
        else:
            s = Structure(seq1, seq2, math.inf, math.inf, math.inf, math.inf, sodium, magnesium, temperature, concentration)
        
        return [s]

    # TODO: Either remove or re-incorporate using Primer3 to create 'PrimerPair' objects
    def scan_sequence(self, seq, primer_size=(18,26), amplicon_size=(50,60)):
        if (__name__ == "__main__"):
            from oligo import Primer, PrimerPair
        else:
            from .oligo import Primer, PrimerPair

        import primer3

        number_records = 20
        
        mv_conc = 50.0 # in mM
        dv_conc = 0.0 # in mM
        dntp_conc = 0.6 # in mM
        dna_conc = 250.0 # in nM
        
        temperature = 25
        
        
        seq_args = {
            'SEQUENCE_ID': 'TEST',
            'SEQUENCE_TEMPLATE': seq
        }
        
        global_args = {
            # Parameters for design
            #'PRIMER_TASK': 'generic',
            #  generic
            #  check_primers
            #  pick_primer_list
            #  pick_sequencing_primers
            #  pick_cloning_primers
            #  pick_discriminative_primers
            'PRIMER_PICK_LEFT_PRIMER':               1,
            'PRIMER_PICK_INTERNAL_OLIGO':            0,
            'PRIMER_PICK_RIGHT_PRIMER':              1,
            
            'PRIMER_NUM_RETURN':        number_records, # in output records
            'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT':  1,
            
            # Parameters for LEFT/RIGHT oligos
            'PRIMER_OPT_SIZE':                      20, # in nt
            'PRIMER_MIN_SIZE':          primer_size[0], # in nt
            'PRIMER_MAX_SIZE':          primer_size[1], # in nt
            'PRIMER_MAX_POLY_X':                     4, # in nt
            'PRIMER_MAX_NS_ACCEPTED':                0, # in nt
            'PRIMER_MAX_END_GC':                     3, # in nt
            'PRIMER_GC_CLAMP':                       1, # in nt (1 or more trailing G or C nt)
            'PRIMER_PRODUCT_OPT_SIZE':               0, # in nt (0 means don't prefer any one size)
            #'PRIMER_PRODUCT_SIZE_RANGE': [(31,40),(41,50),(51,60),(61,70)], # in nt
            'PRIMER_PRODUCT_SIZE_RANGE': [amplicon_size], # in nt
            
            'PRIMER_OPT_TM':                      60.0, # in degrees C
            'PRIMER_MIN_TM':                      55.0, # in degrees C
            'PRIMER_MAX_TM':                      65.0, # in degrees C
            'PRIMER_PAIR_MAX_DIFF_TM':             2.0, # in degrees C
            
            'PRIMER_OPT_GC_PERCENT':              50.0, # in percent
            'PRIMER_MIN_GC':                      40.0, # in percent
            'PRIMER_MAX_GC':                      60.0, # in percent
            
            'PRIMER_SALT_MONOVALENT':          mv_conc, # in mM
            'PRIMER_SALT_DIVALENT':            dv_conc, # in mM
            'PRIMER_DNA_CONC':                dna_conc, # in nM (Not the concentration of oligos in the reaction mix but of those annealing to template.)
            'PRIMER_DNTP_CONC':              dntp_conc, # in mM
            
            'PRIMER_MAX_SELF_ANY':                   8, # alignment score
            'PRIMER_MAX_SELF_END':                   3, # alignment score
            'PRIMER_PAIR_MAX_COMPL_ANY':             8, # alignment score
            'PRIMER_PAIR_MAX_COMPL_END':             3, # alignment score
            
            'PRIMER_MAX_SELF_ANY_TH':             45.0, # degrees C
            'PRIMER_MAX_SELF_END_TH':             35.0, # degrees C
            'PRIMER_PAIR_MAX_COMPL_ANY_TH':       45.0, # degrees C
            'PRIMER_PAIR_MAX_COMPL_END_TH':       35.0, # degrees C
            'PRIMER_MAX_HAIRPIN_TH':              35.0, # degrees C
            
            # parameters for INTERNAL oligos
            'PRIMER_INTERNAL_OPT_SIZE':             20, # in nt
            'PRIMER_INTERNAL_MIN_SIZE': primer_size[0], # in nt
            'PRIMER_INTERNAL_MAX_SIZE': primer_size[1], # in nt
            'PRIMER_INTERNAL_MAX_POLY_X':            4, # in nt
            'PRIMER_INTERNAL_MAX_NS_ACCEPTED':       0, # in nt
            
            'PRIMER_INTERNAL_OPT_TM':             60.0, # in degrees C
            'PRIMER_INTERNAL_MIN_TM':             55.0, # in degrees C
            'PRIMER_INTERNAL_MAX_TM':             65.0, # in degrees C
            
            'PRIMER_INTERNAL_OPT_GC_PERCENT':     50.0, # in percent
            'PRIMER_INTERNAL_MIN_GC':             40.0, # in percent
            'PRIMER_INTERNAL_MAX_GC':             60.0, # in percent
            
            'PRIMER_INTERNAL_SALT_MONOVALENT': mv_conc, # in mM
            'PRIMER_INTERNAL_SALT_DIVALENT ':  dv_conc, # in mM
            'PRIMER_INTERNAL_DNTP_CONC':     dntp_conc, # in mM
            'PRIMER_INTERNAL_DNA_CONC':       dna_conc, # in nM
            
            'PRIMER_INTERNAL_MAX_SELF_ANY':          8, # alignment score
            'PRIMER_INTERNAL_MAX_SELF_END':          3, # alignment score
            
            'PRIMER_INTERNAL_MAX_SELF_ANY_TH':    45.0, # degrees C
            'PRIMER_INTERNAL_MAX_SELF_END_TH':    35.0, # degrees C
            'PRIMER_INTERNAL_MAX_HAIRPIN_TH':     24.0, # degrees C
        }
        primers = primer3.bindings.designPrimers(seq_args, global_args)
        records = []
        found_records = primers['PRIMER_PAIR_NUM_RETURNED']
        for i in range(min(found_records, number_records)):
            n = str(i)
            headers = ['^PRIMER_PAIR_'+n+'_', '^PRIMER_LEFT_'+n+'_', '^PRIMER_RIGHT_'+n+'_']
            rr = '|'.join(headers)
            records.append({})
            for p in primers:
                m = regex.search(rr, p)
                if m:
                    records[-1][p] = primers[p]
            #print(p, primers[p])
        
        outputs = []
        for i, r in enumerate(records):
            n = str(i)
            
            prefix = 'PRIMER_LEFT_'+n+'_'
            p_seq = r[prefix+'SEQUENCE']
            left_seq = p_seq
            p_pos = seq.find(p_seq) # 0-based indexing
            p_hairpin = primer3.calcHairpin(p_seq, mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc, dna_conc=dna_conc, temp_c=temperature)
            p_homodimer = primer3.calcHomodimer(p_seq, mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc, dna_conc=dna_conc, temp_c=temperature)
            p_rc = primer3.calcHeterodimer(p_seq, rc(p_seq), mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc, dna_conc=dna_conc, temp_c=temperature)
            p_gc = (p_seq.count('C') + p_seq.count('G'))/len(p_seq)
            o_left = Primer(p_seq, p_pos, '+', p_hairpin, p_homodimer, p_rc, p_gc)
            
            prefix = 'PRIMER_RIGHT_'+n+'_'
            p_seq = r[prefix+'SEQUENCE']
            right_seq = p_seq
            p_pos = seq.find(rc(p_seq)) # 0-based indexing
            p_hairpin = primer3.calcHairpin(p_seq, mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc, dna_conc=dna_conc, temp_c=temperature)
            p_homodimer = primer3.calcHomodimer(p_seq, mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc, dna_conc=dna_conc, temp_c=temperature)
            p_rc = primer3.calcHeterodimer(p_seq, rc(p_seq), mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc, dna_conc=dna_conc, temp_c=temperature)
            p_gc = (p_seq.count('C') + p_seq.count('G'))/len(p_seq)
            o_right = Primer(p_seq, p_pos, '-', p_hairpin, p_homodimer, p_rc, p_gc)
            
            p_heterodimer = primer3.calcHeterodimer(left_seq, right_seq, mv_conc=mv_conc, dv_conc=dv_conc, dntp_conc=dntp_conc, dna_conc=dna_conc, temp_c=temperature)
            o_het = PrimerPair(o_left, o_right, p_heterodimer)
            
            
            # ThermoResult object:
            #  dg               deltaG (Gibbs free energy) of the structure (cal/mol)
            #  dh               deltaH (entropy) of the structure (cal/mol)
            #  ds               deltaS (enthalpy) of the structure (cal/K*mol)
            #  structure_found  Whether or not a structure (hairpin, dimer, etc) was found as a result of the calculation.
            #  tm               Melting temperature of the structure in deg. C
            
            
            #print(i, left_seq, right_seq, min([left_hairpin.dg, right_hairpin.dg, left_homodimer.dg, right_homodimer.dg, heterodimer.dg]))
            
            outputs.append(o_het)
        
        return outputs

def old_test():
    """Code to test the classes and functions in 'source/oligos/_primer3.py'"""
    
    C = Primer3()
    print("===", C.name, "===")
    seq = 'TTCGTGTAGGATCACACCCGTTCCAAGATGTATAATCAGGAGACTCTTACGGTTACGAGGGACCCTCATCCAAGGACTCTAGGTGCAAAGTAACCGGTGG' # 2 pairs
    
    primer_pairs = C.scan_sequence(seq)
    for pp in primer_pairs:
        #print(pp, pp.forward_primer.strand, pp.reverse_primer.strand)
        print(pp)
    
    print(seq)
    for pp in primer_pairs:
        print(' '*pp.forward_primer.position + pp.get_formatted())

def test():
    """Code to test the classes and functions in 'source/oligos/_unafold.py'"""
    
    C = Primer3()
    print("===", C.name, "===")
    
    a = 'GAAATCGCTTAGCGCGAACTCAGACCAT'
    b = 'CCTAGCTATTTAATAAATC'
    c = 'TTCTCCACTTCCATCACCGT'
    
    tempdir = '.'
    
    print('Hairpin: {}'.format(repr(a)))
    for s in C.find_structures(tempdir, a, None):
        print('', s)
    print('Homodimer: {} {}'.format(repr(a), repr(a)))
    for s in C.find_structures(tempdir, a, a):
        print('', s)
    print('Heterodimer: {} {}'.format(repr(a), repr(b)))
    for s in C.find_structures(tempdir, a, b):
        print('', s)
    print('Reverse-complements: {} {}'.format(repr(a), repr(rc(a))))
    for s in C.find_structures(tempdir, a, rc(a)):
        print('', s)

if (__name__ == "__main__"):
    test()
