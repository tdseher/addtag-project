#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/man.py

# List general Python imports
import sys
import os
import subprocess

__text__ = '''\
% ADDTAG(1) Version 1
% Thaddeus Seher (programmer); Aaron Hernday (PI)
% February 4, 2017
'''

def get_readme_text(path):
    """Gets the text from 'README.md'"""
    text = ''
    with open(os.path.join(path, 'README.md'), 'r') as flo:
        text = flo.read()
    return text

def get_help_text(path):
    """Runs 'addtag --help'"""
    commands = [os.path.join(path, 'addtag'), '--help']
    cp = subprocess.run(commands, shell=False, check=True, stdout=subprocess.PIPE)
    return cp.stdout.decode()

def run_pandoc(path, text):
    commands = ['pandoc', '-s', '-t', 'man', '-o', os.path.join(path, 'addtag.1')]
    #cp = subprocess.run(commands, shell=False, check=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    cp = subprocess.run(commands, shell=False, check=True, input=text.encode(), stdout=subprocess.PIPE)
    #stdout, stderr = cp.communicate(text)
    #stdout.decode()

def convert_readme(text):
    # Raise level 3 headers '###' to level 1 '#'
    # Strip URLs?
    # parse the authors section to automatically get the authors
    pass

def convert_usage(text):
    # Remove usage overview
    # Remove DNA/RNA/protein diagrams
    # For each usage section header, add '#' so it becomes a man section header
    # Format each motif
    # Format outputs
    # Format arguments, appropriately bold/underlining
    pass

def main():
    '''
    Script that takes the 'README.md' and the output of 'addtag --help'
    and formats them with the proper markdown. Then it passes them into
    pandoc (http://pandoc.org/demos.html, http://pandoc.org/demo/pandoc.1.md)
    which outputs the properly-formatted 'addtag.1' man file
      program    type  input        output
    $ pandoc -s -t man addtag.1.md -o addtag.1
    '''
    
    addtag_path = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__), '..')))
    readme_text = get_readme_text(addtag_path)
    help_text = get_help_text(addtag_path)
    
    all_text = __text__ + "\n" + readme_text + "\n" + help_text
    
    run_pandoc(addtag_path, all_text)

if (__name__ == '__main__'):
    main()
