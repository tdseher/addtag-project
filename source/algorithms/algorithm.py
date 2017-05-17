#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/algorithm.py

class Algorithm(object): # Name of the subclass
    """
    General class that should be subclassed when adding a new scoring type
    """
    def __init__(self, 
        name,
        author,
        year,
        citation=None,
        off_target=False,
        prefilter=False,
        postfilter=False,
        minimum=0.0,
        maximum=100.0,
        default=None
    ):
        """
        Specify general information regarding this new instance of the 
        Algorithm class.
        """
        self.name = name             # Unique name for the algorithm (str). No other Algorithm objects should have this name.
        self.author = author         # Author of the algorithm (str)
        self.year = year             # Year algorithm published (int)
        self.citation = citation     # Citation (None, str)
        self.off_target = off_target # Calculate off-target/Guide/Efficiency score (True, False)
        self.prefilter = prefilter   # Use algorithm to filter sequences before alignment
        self.postfilter = postfilter # Use algorithm to filter sequences after alignment
        self.minimum = minimum       # Minimum score to be included in pre- and post- alignment filters (float)
        self.maximum = maximum       # Maximum score to be included in pre- and post- alignment filters (float)
        self.default = default       # If defined, then this will be the on-target default score (None, float)
    
    def calculate(self, *args, **kwargs):
        """
        Defines interface for scoring.
        Overload this method.
        
        See SingleSequenceAlgorithm and PairedSequenceAlgorithm for arguments.
        Returns output of self.score()
        """
        
        return self.score(args, kwargs)
    
    def score(self, *args, **kwargs):
        """
        Overload this method.
        Returns a float
        """
        return 0.0
    
    def __repr__(self):
        """
        Return the string representation of the Algorithm
        """
        return self.__class__.__name__ + '(' + ', '.join(['name='+repr(self.name), 'author='+repr(self.author), 'year='+repr(self.year)]) + ')'

class SingleSequenceAlgorithm(Algorithm):
    def calculate(self, intended, *args, **kwargs):
        """
        Provides universal interface for single sequence scoring.
        Calculate the score for the input sequence.
        
        Overload this method.
        
        Please include **kwargs when you overload this function to catch
        unused arguments passed in.
        
        Returns:
         float
        Arguments:
         sequence, target, pam, upstream, downstream (for on-target)
        Optional arguments:
         iupac (True, False)
         distance (int)
        """
        sequence, target, pam, upstream, downstream = intended
        
        return self.score(target, pam)

class PairedSequenceAlgorithm(Algorithm):
    def calculate(self, intended, potential, *args, **kwargs):
        """
        Provides universal interface for scoring a pair of sequences.
        Calculate the score for the input sequences.
        
        Overload this method.
        
        Please include **kwargs when you overload this function to catch
        unused arguments passed in.
        
        Returns:
         float
        Arguments:
         sequence, target, pam, upstream, downstream (for on-target)
         sequence, target, pam, upstream, downstream (for off-target)
        Optional arguments:
         iupac (True, False)
         distance (int)
        """
        on_sequence, on_target, on_pam, on_upstream, on_downstream = intended
        off_sequence, off_target, off_pam, off_upstream, off_downstream = potential
        
        return self.score(on_target, off_target)

class BatchedSingleSequenceAlgorithm(Algorithm):
    def calculate(self, batch, *args, **kwargs):
        """
        Provides universal interface for scoring single sequences all at once.
        Calculate the score for each input sequence.
        
        Overload this method.
        
        Please include **kwargs when you overload this function to catch
        unused arguments passed in.
        
        Returns:
         list of floats
        Arguments:
         [(sequence, target, pam, upstream, downstream), ...]
         
        Optional arguments:
         iupac (True, False)
         distance (int)
        """
        scores = []
        for sequence, target, pam, upstream, downstream in batch:
            scores.append(self.score(target, pam))
        return scores
