#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/xu/__init__.py

# Metric for "efficiency"
#  Han Xu
#   Range ~ -2 - +2. Aka 'SSC score'. Linear regression model trained on
#   data from >1000 genes in human KBM7/HL60 cells (Wang et al) and mouse
#   (Koike-Yusa et al.). Delivery: lentivirus. Ranges mostly -2 to +2.
#   See Xu et al.: http://genome.cshlp.org/content/early/2015/06/10/gr.191452.115
# 
#   Apparently also called: Spacer Scoring for CRISPR (SSC)
#   See Xu et al, Gen Res 2015, PMID 26063738, http://crispr.dfci.harvard.edu/SSC/

# This program scans spacer sequence of CRISPR to predict effectiveness of guide RNA.

# Either port the algorithm to Python, or run SSC externally, and import the results




