#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/oof/__init__.py

# Metric for "efficiency"
# Out-of-frame
#  Range: 0-100. Predicts the percentage of clones that will carry
#  out-of-frame deletions, based on the micro-homology in the sequence
#  flanking the target site. See Bae et al.
# (http://www.nature.com/nmeth/journal/v11/n7/full/nmeth.3015.html)
#
# Microhomology-based choice of Cas9 nuclease target sites
# Sangsu Bae, Jiyeon Kweon, Heon Seok Kim & Jin-Soo Kim 
# Nature Methods volume 11, pages705–706 (2014)

# For an implementation example, see "cas-designer" (http://www.rgenome.net/cas-designer/portable) Python script
