#!/usr/bin/env python3

"""Program for identifying unique endogenous gRNA sites
and creating unique synthetic gRNA sites."""

# Import standard packages
import sys
import argparse
import textwrap
import time
import string
import math

# Import non-standard packages
import regex

# Define meta variables
__author__ = "Thaddeus D. Seher (@tdseher), & Aaron Hernday"
__date__ = "2017-02-03"
__version__ = "1"

__description__ = """\
Program for identifying unique endogenous gRNA sites
and creating unique synthetic gRNA sites. \
Copyright (c) {__date__} {__author__}. All rights reserved. \
""".format(**locals())



def parse_arguments():
    # Create the parser
    parser = argparse.ArgumentParser(
        description="""Counts how many times each read maps with certain \
                    orientations in a set of SAM files. \
                    Copyright 2015 Thaddeus Seher.""",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Add mandatory arguments
    parser.add_argument("fasta", type=str, help="FASTA file containing all contigs to derive gRNA sites for. All FASTA sequences should have unique primary headers (everything between the '>' and the first ' ' should be unique).")
    parser.add_argument("gff", type=str, help="GFF file specifying chromosomal features")
    
    # Add optional arguments
    parser.add_argument("--min_contig_edge_distance", metavar="N", type=int, default=500, help="Minimum distance from contig edge a site can be found")
    #parser.add_argument("--features", metavar="FEATURE,FEATURE", type=str, default="ORF", help="Features to design gRNA sites against")
    parser.add_argument("--features", metavar="FEATURE", type=str, nargs="+", default=["ORF"], help="Features to design gRNA sites against")
    parser.add_argument("--min_donor_length", metavar="N", type=int, default=80, help="The minimum length of the final computed donor DNA for each site")
    parser.add_argument("--max_donor_length", metavar="N", type=int, default=90, help="The maximum length of the final computed donor DNA for each site")
    parser.add_argument("--min_feature_edge_distance", metavar="N", type=int, default=24, help="The minimum distance a gRNA site can be from the edge of the feature. If negative, the maximum distance a gRNA site can be outside the feature.")
    parser.add_argument("--min_donor_insertions", metavar="N", type=int, default=2, help="The uniqueness of final donor DNA compared to the rest of the genome")
    parser.add_argument("--min_donor_deletions", metavar="N", type=int, default=2, help="The uniqueness of final donor DNA compared to the rest of the genome")
    parser.add_argument("--min_donor_mismatches", metavar="N", type=int, default=2, help="The uniqueness of final donor DNA compared to the rest of the genome")
    parser.add_argument("--min_donor_differences", metavar="N", type=int, default=3, help="The uniqueness of final donor DNA compared to the rest of the genome")
    parser.add_argument("--min_donor_distance", metavar="N", type=int, default=36, help="The minimum distance in bp a difference can exist from the edge of donor DNA")
    parser.add_argument("--min_target_length", metavar="N", type=int, default=23, help="The minimum length of the 'target'/'spacer'/gRNA site")
    parser.add_argument("--max_target_length", metavar="N", type=int, default=28, help="The maximum length of the 'target'/'spacer'/gRNA site")
    parser.add_argument("--strands", type=str, choices=["+", "-", "+/-"], default="+/-", help="Strands to search for gRNAs")
    
    return parser.parse_args()

def list_pam_sites():
    # Code taken from
    #  https://github.com/maximilianh/crisporWebsite/crispor.py
    pamDesc = [ ('NGG','20bp-NGG - Cas9 Streptococcus Pyogenes and Cas9-HF1'),
         ('TTTN','TTTN-23bp - Cpf1 Acidaminococcus / Lachnospiraceae'),
         #('TTN','TTN-23bp - Cpf1 F. Novicida'), # Jean-Paul: various people have shown that it's not usable yet
         ('NGA','20bp-NGA - Cas9 S. Pyogenes mutant VQR'),
         ('NGCG','20bp-NGCG - Cas9 S. Pyogenes mutant VRER'),
         ('NNAGAA','20bp-NNAGAA - Cas9 S. Thermophilus'),
         ('NGGNG','20bp-NGGNG - Cas9 S. Thermophilus'),
         ('NNGRRT','21bp-NNG(A/G)(A/G)T - Cas9 S. Aureus'),
         ('NNNNGMTT','20bp-NNNNG(A/C)TT - Cas9 N. Meningitidis'),
         ('NNNNACA','20bp-NNNNACA - Cas9 Campylobacter jejuni'),
       ]

def read_fasta_file(args):
    """Load contig sequences from file into dict()
    Primary sequence headers must be unique
    """
    contigs = {}
    with open(args.fasta, 'r') as flo:
        name = None
        for line in flo:
            line = line.rstrip()
            if line.startswith('>'):
                name = regex.split(r'\s+', line[1:], 1)[0]
                contigs[name] = ''
            else:
                # Handle malformatted FASTA
                if ((name == None) or (name == "")):
                    raise ValueError('FASTA file malformatted')
                else:
                    contigs[name] += line
    
    return contigs

def lcs(string1, string2):
    """Find the longest common substring between two strings"""
    import difflib
    
    matcher = difflib.SequenceMatcher(None, string1, string2, True)
    match = matcher.find_longest_match(0, len(string1), 0, len(string2))
    # Match(a=0, b=15, size=9)
    return match
    # print(string1[match.a: match.a + match.size])
    # print(string2[match.b: match.b + match.size])

def rc(seq, kind="dna"):
    """Returns the reverse-complement of a string containing DNA or RNA characters"""
    if (kind == "dna"):
        complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB') # exclude ws, WS
    elif (kind == "rna"):
        complements = string.maketrans('acgurymkbdhvACGTRYMKBDHV', 'ugcayrkmvhdbTGCAYRKMVHDB') # exclude ws, WS
    else:
        raise ValueError("'" + str(kind) + "' is an invalid argument for rc()")
    return seq.translate(complements)[::-1]

def on_target_score_2014(seq):
    """Function to calculate the sgRNA on-target efficacy score,
    as defined in Doench et al 2014"""
    # Code retrieved from:
    #  https://github.com/maximilianh/crisporWebsite/doenchScore.py
    # This is a 13 line python function to calculate the sgRNA on-target efficacy score from the article
    # "Rational design of highly active sgRNAs for CRISPR-Cas9-mediated gene inactivation"
    # by J Doench et al. 2014
    # The authors' web tool is available at http://www.broadinstitute.org/rnai/public/analysis-tools/sgrna-design
    # Thanks to Cameron Mac Pherson at Pasteur Paris for fixing my original version. Maximilian Haeussler 2014
    
    # Code for testing this function:
    # print "expected result:", 0.713089368437
    # print on_target_score_2014("TATAGCTGCGATCTGAGGTAGGGAGGGACC")
    # print "expected result:", 0.0189838463593
    # print on_target_score_2014("TCCGCACCTGTCACGGTCGGGGCTTGGCGC")
    
    # pasted/typed table from PDF and converted to zero-based positions
    params = [
        (1,'G',-0.2753771),(2,'A',-0.3238875),(2,'C',0.17212887),(3,'C',-0.1006662),
        (4,'C',-0.2018029),(4,'G',0.24595663),(5,'A',0.03644004),(5,'C',0.09837684),
        (6,'C',-0.7411813),(6,'G',-0.3932644),(11,'A',-0.466099),(14,'A',0.08537695),
        (14,'C',-0.013814),(15,'A',0.27262051),(15,'C',-0.1190226),(15,'T',-0.2859442),
        (16,'A',0.09745459),(16,'G',-0.1755462),(17,'C',-0.3457955),(17,'G',-0.6780964),
        (18,'A',0.22508903),(18,'C',-0.5077941),(19,'G',-0.4173736),(19,'T',-0.054307),
        (20,'G',0.37989937),(20,'T',-0.0907126),(21,'C',0.05782332),(21,'T',-0.5305673),
        (22,'T',-0.8770074),(23,'C',-0.8762358),(23,'G',0.27891626),(23,'T',-0.4031022),
        (24,'A',-0.0773007),(24,'C',0.28793562),(24,'T',-0.2216372),(27,'G',-0.6890167),
        (27,'T',0.11787758),(28,'C',-0.1604453),(29,'G',0.38634258),(1,'GT',-0.6257787),
        (4,'GC',0.30004332),(5,'AA',-0.8348362),(5,'TA',0.76062777),(6,'GG',-0.4908167),
        (11,'GG',-1.5169074),(11,'TA',0.7092612),(11,'TC',0.49629861),(11,'TT',-0.5868739),
        (12,'GG',-0.3345637),(13,'GA',0.76384993),(13,'GC',-0.5370252),(16,'TG',-0.7981461),
        (18,'GG',-0.6668087),(18,'TC',0.35318325),(19,'CC',0.74807209),(19,'TG',-0.3672668),
        (20,'AC',0.56820913),(20,'CG',0.32907207),(20,'GA',-0.8364568),(20,'GG',-0.7822076),
        (21,'TC',-1.029693),(22,'CG',0.85619782),(22,'CT',-0.4632077),(23,'AA',-0.5794924),
        (23,'AG',0.64907554),(24,'AG',-0.0773007),(24,'CG',0.28793562),(24,'TG',-0.2216372),
        (26,'GT',0.11787758),(28,'GG',-0.69774)
    ]
    
    intercept =  0.59763615
    gcHigh = -0.1665878
    gcLow = -0.2026259
    
    score = intercept
    
    guideSeq = seq[4:24]
    gcCount = guideSeq.count("G") + guideSeq.count("C")
    if gcCount <= 10:
        gcWeight = gcLow
    if gcCount > 10:
        gcWeight = gcHigh
    score += abs(10-gcCount)*gcWeight
    
    for pos, modelSeq, weight in params:
        subSeq = seq[pos:pos+len(modelSeq)]
        if subSeq==modelSeq:
            score += weight
    return 1.0/(1.0+math.exp(-score))

def get_mm_scores(file_path):
    """Load the mm scores defined by Doench et al (2016)"""
    # by defauly, file_path should equal 'mismatch_score.pkl'
    import fractions
    mm_scores = pickle.load(open(file_path,'rb'))
    for k in mm_scores:
        print(k, mm_scores[k], fractions.Fraction(mm_scores[k]).limit_denominator())
    # This looks to be 5-dimensional data:
    #  1  2  3 4/5                 4  5
    # rG:dA, 5 0.3                 3/10
    # rC:dT,12 0.7142857140000001  5/ 7
    # rC:dT,18 0.538461538         7/13
    # rU:dG,14 0.28571428600000004 2/ 7
    # rC:dA,13 0.7                 7/10
    # rC:dC,18 0.133333333         2/15
    # ...etc...
    return mm_scores

def on_target_score_2016():
    """
    Also called CFD score.
    
    The on-target score is from Doench, Fusi et al. (2016)
    doi:10.1038/nbt.3437 and measures activity at the target location.
    Higher scores give higher confidence that the guide will be active.
    Scores range from 0-100, and should be used to rank guides relative to
    each other.
    
    The on-target score represents the cleavage efficiency of Cas96.
    You can think of the score as the probability a given gRNA will be in
    top 20% of cleavage activity. Note that the scoring system is not linear,
    and only 5% of gRNAs receive a score of 60 or higher.
    """
    
    # Unpickle mismatch scores and PAM scores
    # def get_mm_pam_scores():
    try:
        mm_scores = pickle.load(open('mismatch_score.pkl','rb'))
        pam_scores = pickle.load(open('pam_scores.pkl','rb'))
        #return (mm_scores,pam_scores)
    except: 
        raise Exception("Could not find file with mismatch scores or PAM scores")
    
   # Calculates CFD score
   # def calc_cfd(wt,sg,pam):
    #mm_scores,pam_scores = get_mm_pam_scores()
    score = 1
    sg = sg.replace('T','U')
    wt = wt.replace('T','U')
    s_list = list(sg)
    wt_list = list(wt)
    for i,sl in enumerate(s_list):
        if wt_list[i] == sl:
            score*=1
        else:
            key = 'r'+wt_list[i]+':d'+revcom(sl)+','+str(i+1)
            score*= mm_scores[key]
    score*=pam_scores[pam]
    return (score)

def off_target_score():
    """The off-target score is from Hsu et al. (2013) doi:10.1038/nbt.2647
    and measures how specific the guide is to the target location. Guides with
    scores over 50 are good enough to be candidates. Higher scores are better.
    Scores range from 0-100, and should be used to rank guides relative
    to each other.
    
    The off-target score tells you the inverse probability of Cas9 off-target
    binding. A higher score means the sequence has less chance to bind to
    sequences in the rest of the genome.
    
    Returns the off-target score and the list of potential off-target sequences
    and their genomic coordinates.
    """
    pass

def scores():
    pass
    # Metrics for "efficiency"
    #  Fusi
    #   Range: 0-100. Boosted Regression Tree model, trained on data produced
    #   by Doench et al (881 guides, MOLM13/NB4/TF1 cells + unpublished
    #   additional data). Delivery: lentivirus. See Fusi et al.:
    #   http://biorxiv.org/content/early/2015/06/26/021568
    #  Chari
    #   Range: 0-100. Support Vector Machine, converted to rank-percent,
    #   trained on data from 1235 guides targeting sequences that were also
    #   transfected with a lentivirus into human 293T cells. See Chari et al.:
    #   http://www.nature.com/nmeth/journal/v12/n9/abs/nmeth.3473.html
    #  Xu
    #   Range ~ -2 - +2. Aka 'SSC score'. Linear regression model trained on
    #   data from >1000 genes in human KBM7/HL60 cells (Wang et al) and mouse
    #   (Koike-Yusa et al.). Delivery: lentivirus. Ranges mostly -2 to +2.
    #   See Xu et al.: http://genome.cshlp.org/content/early/2015/06/10/gr.191452.115
    #  Doench
    #   Range: 0-100. Linear regression model trained on 880 guides transfected
    #   into human MOLM13/NB4/TF1 cells (three genes) and mouse cells
    #   (six genes). Delivery: lentivirus. The Fusi score can be considered an
    #   updated version this score, as their training data overlaps a lot.
    #   See Doench et al.: http://www.nature.com/nbt/journal/v32/n12/full/nbt.3026.html
    #  Wang
    #   Range: 0-100. SVM model trained on human cell culture data on guides
    #   from >1000 genes. The Xu score can be considered an updated version of
    #   this score, as the training data overlaps a lot. Delivery: lentivirus.
    #   See Wang et al.: http://http//www.ncbi.nlm.nih.gov/pmc/articles/PMC3972032/
    #  Moreno-Mateos
    #   Range: mostly 0-100. Linear regression model, trained on data from 1000
    #   guides on >100 genes, from zebrafish 1-cell stage embryos injected with
    #   mRNA. See Moreno-Mateos et al.:
    #   http://www.nature.com/nmeth/journal/v12/n10/full/nmeth.3543.html
    #  Housden
    #   Range: ~ 1-10. Weight matrix model trained on data from Drosophila mRNA
    #   injections. See Housden et al.:
    #   http://stke.sciencemag.org/content/8/393/rs9.long
    #  Prox GC, -GC
    #   This column shows two heuristics based on observations rather than
    #   computational models: Ren et al 2014 (http://www.cell.com/cell-reports/fulltext/S2211-1247(14)00827-4)
    #   obtained the highest cleavage in Drosophila when the final 6bp
    #   contained >= 4 GCs, based on data from 39 guides. Farboud et al.
    #   (http://www.genetics.org/content/early/2015/02/18/genetics.115.175166.abstract)
    #   obtained the highest cleavage in C. elegans for the 10 guides that
    #   ended with -GG, out of the 50 guides they tested.
    #   This field contains + if the final GC count is >= 4 and GG if the guide ends with GG.
    #   -GC citation: https://mcb.berkeley.edu/labs/meyer/publicationpdfs/959.full.pdf
    
    # Out-of-frame
    #  Range: 0-100. Predicts the percentage of clones that will carry
    #  out-of-frame deletions, based on the micro-homology in the sequence
    #  flanking the target site. See Bae et al.
    # (http://www.nature.com/nmeth/journal/v11/n7/full/nmeth.3015.html)
    
    # two types of off-target scores
    #  CFD off-target score
    #  MIT off-target score
    #  Hsu specificity score
    #   The higher the specificity score, the lower are off-target effects in the genome.
    #   The specificity score ranges from 0-100 and measures the uniqueness
    #   of a guide in the genome. See Hsu et al. Nat Biotech 2013.
    #   (http://dx.doi.org/10.1038/nbt.2647) We recommend values >50, where possible.
    
    # Histogram of off-targets:
    #  For each number of mismatches, the number of off-targets is indicated.
    #  Example:
    #   1-3-20-50-60    This means 1 off-target with 0 mismatches, 3 off-targets with 1 mismatch, 20 off-targets with 2 mismatches, etc.
    #   0-2-5-10-20     These are the off-targets that have no mismatches in the 12 bp adjacent to the PAM. These are the most likely off-targets.
    #   
    #   Off-targets are considered if they are flanked by one of the motifs NGG, NAG or NGA.
    

def generate_excise_target(args, feature):
    """Finds the gRNA sequence to cut within the specified places
    on the feature. May target either the + or - strand.
    """
    pass

def generate_revert_target(args, feature):
    """
    Creates the gRNA sequence that targets the excise donor DNA oligo.
    May target either the + or - strand.
    """
    pass

def generate_excise_donor(args, feature, revert_target):
    """Creates the DNA oligo with the structure:
    homology--unique gRNA--homology
    that excises the target feature
    """
    pass

def generate_revert_donor(args, feature):
    """
    Use template DNA sequence to create oligo that will be used for fixing
    the DSB
    """
    # Optionally provide a list of restriction enzyme targeting sites
    # on either side for easier cloning
    pass

def format_output():
    # output should take the following format
    # contig    gRNA-start-pos   gRNA-end-pos    gRNA-strand   features    gRNA-sequence    PAM     on-target-score    off-target-score
    pass

def process(args):
    # The Cas9 cuts 3-4bp upstream of the PAM sequence.
    
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

def main():
    # Get timestamp
    start = time.time()
    
    # Obtain command line arguments
    args = parse_arguments()
    print(args, file=sys.stderr)
    
    # Convert input files to memory structures
    contigs = read_fasta_file(args)
    
    
    # Print time taken for program to complete
    print(time.time()-start, file=sys.stderr)

if (__name__ == "__main__"):
    main()
