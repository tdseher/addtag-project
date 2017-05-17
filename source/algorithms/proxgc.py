#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/proxgc/__init__.py

# Metric for "efficiency"
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