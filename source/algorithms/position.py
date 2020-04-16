#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/algorithms/position.py

# Import non-standard packages
if (__name__ == "__main__"):
    from algorithm import SingleSequenceAlgorithm
else:
    from .algorithm import SingleSequenceAlgorithm

# TODO: Finish adding the Algorithm that weighs the Target based on its position in the Feature
#       If the user prefers the Target to be toward the beginning of the Feature
#       (left if Feature orientation=='+', and right if Feature orientation=='-'),
#       Then Targets toward the beginning will have higher weights than those toward the end.
#       If the user wants to disrupt the beginning of a coding sequence (adding a STOP codon), then the resultant 
#       product would have a more severe mutation than if the disruption were at the end of the CDS.

# TODO: Will need to add the POSITION of the Target within the Feature to the Algorithm inputs.

class Position(SingleSequenceAlgorithm):
    def __init__(self):
        super().__init__(
            name="Position",
            authors=['Seher, Thaddeus D.'],
            title='',
            journal='',
            issuing='',
            year=2017,
            doi='',
            off_target=False,
            prefilter=True,
            postfilter=False,
            minimum=0.0,
            maximum=100.0,
            default=None
        )
    
    def calculate(self, intended, *args, **kwargs):
        '''
        Expects 'position=N', 'feature_length=N', 'feature_orientation=+/-'
        
        :param intended: 
        :param args: 
        :param kwargs: 
        :return: 
        '''
        #sequence, target, pam, upstream, downstream = intended
        if 'position' in kwargs:
            return self.score(kwargs['feature_position'], kwargs['feature_length'], kwargs['feature_orientation'])
        else:
            return None # Using 'None' as the default value might raise errors elsewhere
    
    def score(self, f_position, f_length, f_orientation):
        """
        Caclulates the position in terms of percent.
        """
        #if (0 <= f_position < f_length):
        if (f_orientation == '+'):
            return 100*f_position/f_length
        else:
            return 100*(f_length-f_position)/f_length
        #else:
        #    return None

def test():
    a = ('', 'AAAATTAACTATAGGTAAAG', 'TGG', '', '')
    b = ('', 'AACATCAACTCTAGCTAACG', 'CGG', '', '')
    c = ('', 'AACATCACCTCTGGCTAACG', 'CGG', '', '')
    d = ('', 'ACCACCAACTCTAGCTGACG', 'CGG', '', '')
    e = ('',  'CCACCAACTCTAGCTGACG', 'CGG', '', '')
    f = ('',                    'N', 'CGG', '', '')
    g = ('',                    'H', 'CGG', '', '')
    h = ('', 'AGTCAAAGGAATAGAGAAAC', 'CCAAAC', '', '')
    
    print("=== Position ===")
    C = Position()
    print(C.calculate(a))
    print(C.calculate(a, feature_position=-10, feature_length=1000, feature_orientation='+'))
    print(C.calculate(a, feature_position=10000, feature_length=1000, feature_orientation='+'))
    print(C.calculate(a, feature_position=10, feature_length=1000, feature_orientation='+'))
    print(C.calculate(b, feature_position=100, feature_length=1000, feature_orientation='+'))
    print(C.calculate(c, feature_position=500, feature_length=1000, feature_orientation='+'))
    print(C.calculate(d, feature_position=750, feature_length=1000, feature_orientation='+'))
    print(C.calculate(e, feature_position=10, feature_length=1000, feature_orientation='-'))
    print(C.calculate(f, feature_position=100, feature_length=1000, feature_orientation='-'))
    print(C.calculate(g, feature_position=500, feature_length=1000, feature_orientation='-'))
    print(C.calculate(h, feature_position=750, feature_length=1000, feature_orientation='-'))

if (__name__ == "__main__"):
    test()
