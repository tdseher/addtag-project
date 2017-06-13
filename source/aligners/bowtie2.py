#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/aligners/bowtie2.py

# List general Python imports
import sys
import os
import subprocess

# import non-standard package
import regex

# import AddTag-specific packages
from .. import utils

if (__name__ == "__main__"):
    from aligner import Aligner
else:
    from .aligner import Aligner

class Bowtie2(Aligner):
    def __init__(self):
        super().__init__("Bowtie 2", "Langmead & Salzberg", 2012,
            citation="Langmead & Salzberg. Fast gapped-read alignment with Bowtie 2. Nature Methods 9, 357-359 (2012)"
        )
    
    def index(self, genome_fasta, *args, **kwargs):
        pass
    
    def align(self, index, output, *args, **kwargs):
        pass

def align(output_sam_file, input_query_file, index, threads=(os.cpu_count() or 1), folder=os.getcwd(), options=None):
    """Aligns sequences using Bowtie 2
    
    Description of arguments:
      'name'    - the basename to give the generated files
      'folder' - generated files will be placed here
      'sep'     - the delimiter to use for generating the sequence headers
      'index'   - the index generated by running bowtie2-build on a FASTA
    """
    
    #name = os.path.splitext(os.path.basename(input_query_file))[0]
    #sam_file = os.path.join(folder, name + '.sam')
    #error_file = os.path.join(folder, name + '.err')
    error_file = os.path.splitext(output_sam_file)[0] + '.err'
    
    # if options is specified, then append to the list of options that will be called
    fixed_options_dict = {
        '-p': threads,
        '-S': output_sam_file, #sam_file,
        '-x': index,
        '-U': input_query_file,
        '-k': 20,   # specifies the maximum number of alignments per sequence to return
        '-N': 1,    # Sets the number of mismatches to allowed in a seed alignment
                    # during multiseed alignment. Can be set to 0 or 1. Setting
                    # this higher makes alignment slower (often much slower) but
                    # increases sensitivity.
        '-L': 10,   # Sets the length of the seed substrings to align during
                    # multiseed alignment. Smaller values make alignment slower
                    # but more sensitive.
        '-D': 20,   # Up to <int> consecutive seed extension attempts can "fail"
                    # before Bowtie 2 moves on, using the alignments found so far.
                    # A seed extension "fails" if it does not yield a new best or
                    # a new second-best alignment. This limit is automatically
                    # adjusted up when -k or -a are specified. Default: 15.'
        '-R': 3,    # the maximum number of times Bowtie 2 will "re-seed" reads
                    # with repetitive seeds. When "re-seeding," Bowtie 2 simply
                    # chooses a new set of reads (same length, same number of
                    # mismatches allowed) at different offsets and searches for
                    # more alignments. A read is considered to have repetitive
                    # seeds if the total number of seed hits divided by the
                    # number of seeds that aligned at least once is greater
                    # than 300. Default: 2.
        '-i': 'S,1,0', # split query read into seed index every 1 bp.
        '--n-ceil': 'L,3,0', # Maximum number of Ns is 3
        '-f': None, # query reads are FASTA
        '--end-to-end': None, # disallow soft masking
        '--mp': '3,2', # Set mismatch penalty to 3 Default: '6,2'
        '--rdg': '4,2', # Sets the read gap open (<int1>) and extend (<int2>)
                        # penalties. A read gap of length N gets a penalty of
                        # <int1> + N * <int2>. Default: '5,3'.
        '--rfg': '4,2', # Sets the reference gap open (<int1>) and extend (<int2>)
                        # penalties. A reference gap of length N gets a penalty of
                        # <int1> + N * <int2>. Default: '5,3'.
        '--score-min': 'L,-0.9,-0.9', # Sets function governing minimum alignment
                                      # score needed for an alignment to be considered
                                      # good enough to report, where 'read length' = L:
                                      #   score(L, a, b) = a + b * L
                                      # Default: 'L,-0.6,-0.6'
                                      # Switching to 'L,-0.9,-0.9' returns roughly 9 times
                                      # more alignments than the default
    }
    fixed_options = utils.flatten(fixed_options_dict.items(), remove_none=True)
    
    flat_options = []
    if options:
        if isinstance(options, str):
            flat_options = regex.split(r'\s+', options)
        elif isinstance(options, list):
            flat_options = options
        elif isinstance(options, dict):
            flat_options = utils.flatten(options.items(), remove_none=True)
    
    # Throw an exception if the user tries to set a reserved option
    for o in fixed_options:
        if o in flat_options:
            raise Exception("Command option {!r} is reserved by AddTag".format(o))
    
    command_list = ['bowtie2'] + fixed_options + flat_options
    command_list = list(map(str, command_list))
    command_str = ' '.join(command_list)
    
    # Perform alignment
    print('command: {!r}'.format(command_str))
    with open(error_file, 'w+') as flo:
        cp = subprocess.run(command_list, shell=False, check=True, stdout=flo, stderr=subprocess.STDOUT)
    
    print('FASTA aligned: {!r}'.format(output_sam_file))
    return output_sam_file

def index_reference(fasta, tempdir=os.getcwd(), threads=(os.cpu_count() or 1), options=None, folder=None):
    '''Call bowtie2-build on non-compressed FASTA file.
    
    The indexed FASTA will be stored in 'folder' within 'tempdir' if 'folder'
    is not specified, then the 'fasta' basename will be used as both 'folder'
    and index name.
    '''
    # Find the FASTA basename
    name = os.path.splitext(os.path.basename(fasta))[0]
    
    if not folder:
        folder = name
    # Make the directory if it does not yet exist
    try:
        os.mkdir(os.path.join(tempdir, folder))
    except FileExistsError:
        pass
    
    index_file = os.path.join(tempdir, folder, name)
    error_file = index_file + '.err'
    
    # if options is specified, then append to the list of options that will be called
    fixed_options_dict = {
        '--threads': threads,
        '--offrate': 1,
    }
    fixed_options = utils.flatten(fixed_options_dict.items())
    
    flat_options = []
    if options:
        if isinstance(options, str):
            flat_options = regex.split(r'\s+', options)
        elif isinstance(options, list):
            flat_options = options
        elif isinstance(options, dict):
            flat_options = utils.flatten(options.items())
    
    # Throw an exception if the user tries to set a reserved option
    for o in fixed_options:
        if o in flat_options:
            raise Exception("Command option {!r} is reserved by AddTag".format(o))
    
    # Index the reference, storing the STDOUT and STDERR to the same file
    # If exit status is nonzero, then raise a 'subprocess.CalledProcessError'
    
    # Syntax:
    #   subprocess.run(args, *, stdin=None, input=None, stdout=None,
    #   stderr=None, shell=False, timeout=None, check=False, encoding=None,
    #   errors=None)
    # Run the command described by args. Wait for command to complete,
    # then return a CompletedProcess instance.
    
    # Letting Python handle redirections
    command_list = ['bowtie2-build'] + fixed_options + flat_options + [fasta, index_file]
    command_list = list(map(str, command_list))
    command_str = ' '.join(command_list)
    print('command: {!r}'.format(command_str))
    with open(error_file, 'w+') as flo:
        cp = subprocess.run(command_list, shell=False, check=True, stdout=flo, stderr=subprocess.STDOUT)
    
    # Using OS shell to handle redirections
    #command_list = ['bowtie2-build'] + fixed_options + flat_options + [fasta, index_file, '>', error_file, '2>&1']
    #command_list = list(map(str, command_list))
    #command_str = ' '.join(command_list)
    #cp = subprocess.run(command_str, shell=True, check=True) # works perfectly
    
    print('FASTA file indexed: {!r}'.format(fasta))
    print('Bowtie 2 index: {!r}'.format(index_file))
    
    # Return path to indexed reference
    return index_file
    
def cleanup():
    '''remove all files generated'''
    pass

def old_test():
    """Code to test the classes and functions in 'source/bowtie2.py'"""
    
    print("=== index_reference ===")
    sequences = [
        ('Ca22chr2A_C_albicans_SC5314',  392771,   392791, 'TTACGATCCATTAACTCCTC'),
        ('Ca22chr2A_C_albicans_SC5314',  400860,   400880, 'ACCTCCATTATTTGGTGAGT'),
        ('Ca22chr5A_C_albicans_SC5314', 1150500,  1150520, 'CATGACTAATGGAATCATGG'),
        ('Ca22chrRA_C_albicans_SC5314',   74620,    74640, 'ATATTGTAGTTGTTAAATAC'),
        ('Ca22chr1A_C_albicans_SC5314',    5554,     5574, 'TATTCTCCWGCTTCTTCTTT'),
        ('Ca22chr1A_C_albicans_SC5314', 2902109,  2902129, 'AGTTCTCTCTCTCTYTCTCT'),
    ]
    #contigs = utils.load_fasta_file(sys.argv[1])
    #for s in sequences:
    #    for c in contigs:
    #        i = contigs[c].find(s[3])
    #        if i>=0:
    #            print(s[3], c, i)
    ref = index_reference(sys.argv[1])
    
    print("=== align ===")
    sam = align('temp_ca_alignment', sequences, ref)

def test():
    """Code to test the alignment"""
    print("=== Bowtie 2 ===")
    C = Bowtie2()
    C.index('')
    C.align('', '')

if (__name__ == '__main__'):
    test()