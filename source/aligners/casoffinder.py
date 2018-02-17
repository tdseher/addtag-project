#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/aligners/casoffinder.py

# List general Python imports
import sys
import os
import subprocess
import logging
logger = logging.getLogger(__name__)

if (__name__ == "__main__"):
    from aligner import Aligner
else:
    from .aligner import Aligner

class CasOFFinder(Aligner):
    def __init__(self):
        super().__init__("cas-offinder", "Bae, Park, & Kim", 2014,
            citation="Bae, Park, & Kim. Cas-OFFinder: a fast and versatile algorithm that searches for potential off-target sites of Cas9 RNA-guided endonucleases. Bioinformatics 30(10), 1473-1475 (2014).",
            input='casoffinder',
            output='casoffinder',
            truncated=False,
        )
    
    def index(self, fasta, output_filename, output_folder, threads, *args, **kwargs):
        """
        Cas-OFFinder does not require separate indexing of the reference genome.
        """
        return ''
    
    def align(self, query, subject, output_filename, output_folder, threads, *args, **kwargs):
        """
        """
        # Folder to store all Cas-OFFinder input files
        casoffinder_folder = os.path.join(output_folder, 'casoffinder')
        os.makedirs(casoffinder_folder, exist_ok=True)
        
        qname = os.path.splitext(os.path.basename(query))[0]
        out_file_paths = {}
        
        # Create a single input for each motif to search
        for parsed_motif in kwargs['parsed_motifs']:
            spacers, pams, side = parsed_motif
            for i, s in enumerate(spacers):
                for j, p in enumerate(pams):
                    if (side == '>'):
                        motif = s+p
                    else:
                        motif = p+s
                    
                    # Build special input file for Cas-OFFinder
                    # /path/to/reference/fasta/folder
                    # NNNNNNNNNNNNNNNNNNNNNRG
                    # GGTGTTGAGAGTGCCATCGANNN 5
                    # GTCGTTTTGATTTGGATGTANNN 5
                    # TCAAAATCGTGTTGCGAAATNNN 5
                    input_filename_path = os.path.join(casoffinder_folder, qname+'.casoffinder-input-'+str(i)+str(j))
                    output_filename_path = os.path.join(casoffinder_folder, qname+'.casoffinder-output-'+str(i)+str(j))
                    out_file_paths[(s,p)] = output_filename_path
                    with open(input_filename_path, 'w') as flo:
                        # Add path to reference FASTA
                        print(os.path.dirname(subject), file=flo)
                    
                        # Expand abbreviated motif into the full IUPAC string
                        print(motif, file=flo)
                        
                        # Print all the sequences to search the genome for, as well as the maximum difference
                        for seq in sequences:
                            print(seq + ' ' + str(5), file=flo)
                    
                    # Options to run Cas-OFFinder with
                    options = [
                        input_filename_path,
                        'C',
                        output_filename_path,
                    ]
        
        
        
        
        
        return self.process('cas-offinder', output_filename_path, options)

def test():
    """Code to test the alignment"""
    print("=== Cas-OFFinder ===")
    C = CasOFFinder()
    C.index('')
    C.align('', '')

if (__name__ == '__main__'):
    test()
