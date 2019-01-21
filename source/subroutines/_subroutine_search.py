#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/subroutines/_subroutine_search.py

# Import non-standard packages
import regex

# Import included AddTag-specific modules
from . import subroutine
from .. import utils
from .. import nucleotides

class SearchParser(subroutine.Subroutine):
    def __init__(self, subparsers):
        self.subparsers = subparsers
        
        self.name = 'search'
        self.description = (
            "description:" "\n"
            "  Search input FASTA for DNA sequence, and output a GFF file for use as input" "\n"
            "  for the 'generate' subroutine." "\n"
        )
        self.help = "Search for positions of DNA in FASTA to make GFF file."
        self.epilog = (
            "example:" "\n"
            "  Running AddTag with the following arguments:" "\n"
            "   $ python3 {__program__} {__subroutine__} --fasta genome.fasta --query AGTCGCCAAC" "\n"
            "     CCAATTAGGAG > search.gff" "\n"
            "  " "\n"
            "  Will produce a GFF3 output file 'search.gff' like the following:" "\n"
            "     chr1A	addtag	search	928299	928308	.	-	.	ID=search_0_0" "\n"
            "     chr1B	addtag	search	928338	928347	.	-	.	ID=search_0_1" "\n"
            "     chr4A	addtag	search	3652	3661	.	-	.	ID=search_0_2" "\n"
            "     chr4B	addtag	search	3652	3661	.	-	.	ID=search_0_3" "\n"
            "     chr3A	addtag	search	819381	819391	.	+	.	ID=search_1_0" "\n"
            "     chr3B	addtag	search	819367	819377	.	+	.	ID=search_1_1" "\n"
            "  " "\n"
            "  You can use the '--identifier' option, as follows:" "\n"
            "   $ python3 {__program__} {__subroutine__} --fasta genome.fasta --query AGTCGCCAAC" "\n"
            "     CCAATTAGGAG --identifier feature1 feature2 > search.gff" "\n"
            "  " "\n"
            "  This assigns the identifier to each query:" "\n"
            "     chr1A	addtag	search	928299	928308	.	-	.	ID=feature1_0" "\n"
            "     chr1B	addtag	search	928338	928347	.	-	.	ID=feature1_1" "\n"
            "     chr4A	addtag	search	3652	3661	.	-	.	ID=feature1_2" "\n"
            "     chr4B	addtag	search	3652	3661	.	-	.	ID=feature1_3" "\n"
            "     chr3A	addtag	search	819381	819391	.	+	.	ID=feature2_0" "\n"
            "     chr3B	addtag	search	819367	819377	.	+	.	ID=feature2_1" "\n"
        ).format(**dict(list(subroutine.__dict__.items()) + list({"__subroutine__": self.name}.items()))) # key:value pairs in the latter dict will replace any instances in the prior dict
        
        self.define_parser()
        self.define_arguments()
    
    def define_arguments(self):
        # Add mandatory arguments
        required_group = self.parser.add_argument_group('required arguments')
        required_group.add_argument("--fasta", required=True, nargs="+", metavar="*.fasta", type=str,
            help="FASTA files with contigs to find unique gRNA sites. Input FASTA \
            must not be compressed. User can decide whether ambiguous bases can \
            be chosen for gRNA sites. All FASTA sequences should have unique \
            primary headers (everything between the '>' symbol and the first \
            whitespace should be unique). You should include FASTA of the genome \
            and any plasmids.")
        
        required_group.add_argument("--query", required=True, nargs="+", metavar="SEQUENCE",
            type=str, help="One or more sequences to search.")
        
        # Add optional arguments
        self.parser.add_argument("--identifier", metavar="FEATURE", type=str, nargs="+",
            help="Identifier for input query sequences.")
        
        self.parser.add_argument("--tag", metavar="TAG", type=str, default="ID",
            help="GFF3 attribute tag.")
    
    def compute(self, args):
        """
        Search input FASTA for arbitrary IUPAC sequence
        Print out a GFF3 that can be used as an input to "generate"
        """
        if not args.identifier:
            args.identifier = ["search_"+str(x) for x in range(len(args.query))]
        if (len(args.query) != len(args.identifier)):
            raise Exception("The number of sequences specified in '--query' does not match the number of features in '--identifier'.")
        
        myflags = regex.ENHANCEMATCH | regex.IGNORECASE # regex.ENHANCEMATCH|regex.IGNORECASE|regex.BESTMATCH
        
        contig_sequences = utils.load_multiple_fasta_files(args.fasta)
        
        for i, qseq in enumerate(args.query):
            qpat = nucleotides.build_regex_pattern(qseq) #max_substitutions=0, max_insertions=0, max_deletions=0, max_errors=0, capture=True
            c = regex.compile(qpat, flags=myflags)
            
            rc_qseq = nucleotides.rc(qseq)
            rc_qpat = nucleotides.build_regex_pattern(rc_qseq)
            rc_c = regex.compile(rc_qpat, flags=myflags)
            
            j = 0
            
            for ctg_name, ctg_seq in contig_sequences.items():
                for m in c.finditer(ctg_seq):
                    # contig_id, source, feature_type, start, end, score, strand, frame, attributes
                    contig_id = ctg_name
                    source = "addtag"
                    feature_type = "search"
                    start = m.start()+1
                    end = m.end()
                    score = '.' # <------ replace with E-value or percent identity
                    strand = '+'
                    frame = '.'
                    attributes = args.tag+'='+args.identifier[i] + '_' + str(j)
                    print("\t".join(str(x) for x in [contig_id, source, feature_type, start, end, score, strand, frame, attributes]))
                    j += 1
                
                for m in rc_c.finditer(ctg_seq):
                    contig_id = ctg_name
                    source = "addtag"
                    feature_type = "search"
                    start = m.start()+1
                    end = m.end()
                    score = '.' # <------ replace with E-value or percent identity
                    strand = '-'
                    frame = '.'
                    attributes = args.tag+'='+args.identifier[i] + '_' + str(j)
                    print("\t".join(str(x) for x in [contig_id, source, feature_type, start, end, score, strand, frame, attributes]))
                    j += 1
        
        # End 'compute()'
        