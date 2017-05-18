#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/algorithm.py

class Aligner(object): # Name of the subclass
    """
    General class that should be subclassed when adding a new alignment program
    """
    def __init__(self, 
        name,
        author,
        year,
        citation=None,
        output='sam',
    ):
        """
        Specify general information regarding this new instance of the 
        Aligner class.
        """
        self.name = name         # Unique name for the algorithm (str). No other Aligner objects should have this name.
        self.author = author     # Author of the algorithm (str)
        self.year = year         # Year algorithm published (int)
        self.citation = citation # Citation (None, str)
        self.output = output     # Designate the output format the alignment will be in, so AddTag can select the correct parser
        
    def index(self):
        pass
    
    def align(self):
        pass
