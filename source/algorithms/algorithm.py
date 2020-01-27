#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/algorithm.py

# TODO: Some Algorithms require specific Python+module versions. Thus, wrappers should use virtual environments.
#       For instance: Azimuth requires Python 2.7, but Elevation requires 3.6 (?).
#       Another example: DeepCpf1/CINDEL requires Keras with 'Theano' backend,
#       but ???? requires 'TensorFlow'.
#       To reconcile these, AddTag will need to use wrappers that instantize Python virtual environments
#       See code in '_subroutine_setup.py' for how to do this

# TODO: If Algorithm requires a wrapper (Because it uses a different Python version, for instance),
#       Then this should be indicated in the '__init__()' method
#       Along with the name of the virtual environment 'venv' to load.

class Algorithm(object): # Name of the subclass
    """
    General class that should be subclassed when adding a new scoring type
    """
    def __init__(
        self,
        name,
        authors,
        title,
        journal,
        issuing,
        year,
        doi,
        #citation=None,
        off_target=False,
        on_target=False,
        prefilter=False,
        postfilter=False,
        minimum=0.0,
        maximum=100.0,
        default=None,
        rgns=None,
    ):
        """
        Specify general information regarding this new instance of the 
        Algorithm class.
        """
        self.name = name             # Unique name for the algorithm (str). No other Algorithm objects should have this name.
        self.authors = authors       # Author of the algorithm (str)
        self.title = title
        self.journal = journal
        self.issuing = issuing
        self.year = year             # Year algorithm published (int)
        self.doi = doi
        #self.citation = citation     # Citation (None/str)
        self.off_target = off_target # Use algorithm to calculate off-target/Guide/Efficiency score (True/False)
        self.on_target = on_target   # Algorithm calculates on-target score (True/False)
        self.prefilter = prefilter   # Use algorithm to filter sequences before alignment
        self.postfilter = postfilter # Use algorithm to filter sequences after alignment
        self.minimum = minimum       # Minimum score to be included in pre- and post- alignment filters (float)
        self.maximum = maximum       # Maximum score to be included in pre- and post- alignment filters (float)
        self.default = default       # If defined, then this will be the on-target default score (None/float)
        self.rgns = rgns
        
        self.available = self.is_available()
    
    def is_available(self):
        """
        Determines if the prerequisites for the Algorithm have been met.
        For instance, it will query the PATH to determine if the algorithm software is available.

        Overload this method
        :return: True or False
        """
        return False
    
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
    
    def weight(self, *args, **kwargs):
        """
        Overload this method.
        Returns a float, usually between 0.0 and 1.0
        Argument: a score
        Returns: amount of weight this score should have.
        """
        # TODO: Make weight configurable via a config file?
        # By default, we want the score to be unweighted, we return 1.0
        return 1.0
    
    def __repr__(self):
        """
        Return the string representation of the Algorithm
        """
        return self.__class__.__name__ + '(' + ', '.join(['name='+repr(self.name), 'authors='+repr(self.authors), 'year='+repr(self.year)]) + ')'

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
        # Example:
        # motif       TTTN<N{19}/.{4}\
        # genome      AGCCAACGACTTTAGCTAGCTAAAGGACCTATGCCCATTACATGCCG
        # sequence              TTTAGCTAGCTAAAGGACCTATGCCCA
        # upstream    AGCCAACGAC
        # pam                   TTTA
        # target                    GCTAGCTAAAGGACCTATG
        # cut sites                                   ><  ><
        # sense cuts                                  ><
        # antisense cuts                                  ><
        # double-strand cuts
        # downstream                                       TTACATGCCG
        
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
