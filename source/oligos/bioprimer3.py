#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/oligos/addtagprimer.py

# Import standard packages
import sys
import os

# import non-standard package
import Bio



# See http://biopython.org/DIST/docs/api/


Bio.Emboss.Applications.Primer3Commandline
Bio.Emboss.Applications.PrimerSearchCommandline

Bio.Emboss.Primer3 # Code to parse output from the EMBOSS eprimer3 program.
Bio.Emboss.Primer3.Primers
Bio.Emboss.Primer3.Record
Bio.Emboss.Primer3.parse
Bio.Emboss.Primer3.read



Bio.Emboss.PrimerSearch # Code to interact with the primersearch program from EMBOSS.
Bio.Emboss.PrimerSearch.Amplifier
Bio.Emboss.PrimerSearch.InputRecord
Bio.Emboss.PrimerSearch.OutputRecord
Bio.Emboss.PrimerSearch.read

# Code to use Primer3 to find primers for confirming AddTag




##### Example 1 ########

# We want to create a run of the primer3 program with our parameters
from Bio.Emboss.Applications import Primer3Commandline
primer_cl = Primer3Commandline()
primer_cl.set_parameter("-sequence", "in.pr3")
primer_cl.set_parameter("-outfile", "out.pr3")
primer_cl.set_parameter("-productsizerange", "350,10000")
primer_cl.set_parameter("-target", "%s,%s" % (start, end))
# start and end are the middle region we want to design primers around

# Biopython has a single way to run a program and get the output.
# Python interacts readily with the operating system, making it easy to run other programs from it.
#from Bio.Application import generic_run
#result, messages, errors = generic_run(primer_cl)
stdout, stderr = primer_cl()

# Primer3 generates an output file that looks like:
# 
# # PRIMER3 RESULTS FOR CLB11512.1_789.buhui
#                    Start Len Tm GC% Sequence
#   1 PRODUCT SIZE: 227
#     FORWARD PRIMER 728 20 59.91 50.00 TTCACCTACTGCAAACGCAC
#     REVERSE PRIMER 935 20 59.57 50.00 TTGGTACGTTGTCCATCTCG
# 
# A ton of these files are not easy to deal with.

# We use the Biopython parser to parse these files.
open_outfile = open("out.pr3", "r")
from Bio.Emboss.Primer import Primer3Parser
parser = Primer3Parser()
primer_record = parser.parse(open_outfile)
# The result is that we get the information into a python ready format that we can readily output.

# We write the forward and reverse sequences along with the sequence name to a comma separated value file.

primer = primer_record.primers[0]
print("%s,%s,%s" % (sequence_name, primer.forward_seq, primer.reverse_seq))
# The result is an output full of primers that you can then deal with.




###### Example 2 ############

from Bio.Emboss import Primer3
inputFile = "./wherever/your/outputfileis.out"
with open(inputFile) as fileHandle:
    record = Primer3.parse(fileHandle)
    # XXX check is len>0
    primers = record.next().primers
    numPrimers = len(primers)
    # you should have access to each primer, using a for loop
    # to check how to access the data you care about. For example:



# https://github.com/biopython/biopython/issues/485
#### Example 3 #####

#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from Bio import SeqIO
from Bio.Emboss.Applications import Primer3Commandline
from Bio.Seq import Seq
from Bio import Entrez
Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are

entrez_id="aj627603"

handle = Entrez.efetch(db="nucleotide", id=entrez_id, rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
SeqIO.write(record, ".cache/"+entrez_id+".fasta", "fasta")

primer_cl = Primer3Commandline(sequence=".cache/"+entrez_id+".fasta",auto=True)
primer_cl.outfile = "out.pr3"
primer_cl.numreturn = 3
primer_cl.osize = 20
primer_cl.maxsize = 26
primer_cl.opttm = 60
primer_cl.mintm = 52
primer_cl.mingc = 35
primer_cl.maxgc = 75
primer_cl.psizeopt = 200

primer_cl()

def test():
    """Code to test the classes and functions in 'source/oligos/_unafold.py'"""
    
    C = BioPrimer3()
    print("===", C.name, "===")
    seq = 'TTCGTGTAGGATCACACCCGTTCCAAGATGTATAATCAGGAGACTCTTACGGTTACGAGGGACCCTCATCCAAGGACTCTAGGTGCAAAGTAACCGGTGG' # 2 pairs
    
    primer_pairs = C.scan_sequence(seq)
    for pp in primer_pairs:
        #print(pp, pp.forward_primer.strand, pp.reverse_primer.strand)
        print(pp)
    
    print(seq)
    for pp in primer_pairs:
        print(' '*pp.forward_primer.position + pp.get_formatted())

if (__name__ == "__main__"):
    test()