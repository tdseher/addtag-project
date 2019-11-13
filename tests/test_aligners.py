#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# tests/test_aligners.py

import os

from source import aligners

def test_aligner():
    # TODO: Create input data programmatically instead of loading pre-made files
    data_folder = os.path.dirname(os.path.abspath(__file__))
    subject_file = 'test_subject.fasta'
    query_file = 'test_query.fasta'
    subject_path = os.path.join(data_folder, subject_file)
    query_path = os.path.join(data_folder, query_file)

    threads = 4

    # Index and align
    for aligner in aligners.aligners:
        if aligner.available:
            print("===", aligner.name, "===")
            
            out_folder = 'test_{}_aligner'.format(aligner.name)
            os.makedirs(out_folder, exist_ok=True)
            
            index_path = aligner.index(subject_path, 'test_{}_index'.format(aligner.name), out_folder, threads=threads)
            print('index_path = {}'.format(index_path))
            
            #if (aligner.input == 'fastq'):
            #    query_file = 'test_query.fastq'
            #    query_path = os.path.join(data_folder, query_file)
            
            alignment_path = aligner.align(query_path, index_path, 'test_{}_alignment'.format(aligner.name), out_folder, threads=threads)
            print('alignment_path = {}'.format(alignment_path))
            
            print('records = [')
            for record in aligner.load(alignment_path):
                print('  {}'.format(record))
            print(']')

# End
