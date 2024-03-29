#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/donors.py

# Import standard packages
import logging

# Import non-standard packages
import regex

# Import included AddTag-specific modules
from . import utils
from . import nucleotides
from . import algorithms
from .motifs import OnTargetMotif, OffTargetMotif

logger = logging.getLogger(__name__)

class Alignment(object):
    """Class representing an alignment"""
    __slots__ = [
        'sequence',
        'target',
        'pam',
        'motif',
        'contig',
        'start',
        'end',
        'orientation',
        'upstream',
        'downstream',
        'postfilter',
        'score',
        'action',
    ]
    
    def __init__(self, sequence, target, pam, motif, contig, start, end, orientation, upstream='', downstream=''):
        # Most attributes derived from SAM output
        self.sequence = sequence
        self.target = target
        self.pam = pam
        self.motif = motif
        self.contig = contig
        self.start = start
        self.end = end
        self.orientation = orientation
        
        self.upstream = upstream
        self.downstream = downstream
        
        # Variable to hold whether or not this alignment should be included
        # in the off-target scoring
        self.postfilter = None
        
        # Variable to hold the scores
        self.score = {}
        self.action = 'None'
        
        # TODO: Add the following Algorithms:
        #        lidentities, ridentities, r_score, bae, chari, oof, proxgc, wang, xu
        #self.ridentities = nucleotides.ridentities(self.contig_target, aligned_target)
        #self.r_scores = {}
        #for i in [4, 8, 12, 16]:
        #    self.r_scores[i] = scores.r_score(self.contig_target, aligned_target, i)
    
    def calculate_scores(self, parent, *args, **kwargs):
        # parent = (sequence, target, pam, upstream, downstream)
        this = (self.sequence, self.motif.parsed_list[2], self.target, self.pam, self.upstream, self.downstream)
        postfilter = []
        # TODO: Make this ONLY calculate score if it was selected by user in command line arguments
        for C in algorithms.single_algorithms:
            c_score = C.calculate(this, parsed_motif=self.motif.parsed_list, **kwargs)
            self.score[C.name] = c_score
            if C.postfilter:
                if (C.minimum <= c_score <= C.maximum):
                    postfilter.append(True)
                else:
                    postfilter.append(False)
        for C in algorithms.paired_algorithms:
            c_score = C.calculate(parent, this, parsed_motif=self.motif.parsed_list, **kwargs)
            self.score[C.name] = c_score
            if C.postfilter:
                if (C.minimum <= c_score <= C.maximum):
                    postfilter.append(True)
                else:
                    postfilter.append(False)
        
        # If alignment meets all postfilter criteria, then set it as True
        self.postfilter = all(postfilter) # all([]) returns True
    
    def get_upstream_sequence(self, length, contigs):
        """Returns upstream sequence"""
        return ''
    
    def get_downstream_sequence(self, length, contigs):
        """Returns downstream sequence"""
        return ''
    
    def __repr__(self):
        return self.__class__.__name__ + '(' + ' '.join([
            self.target + '|' + self.pam,
            self.contig + ':' + self.orientation + ':' + str(self.start) + '..' + str(self.end),
            'action=' + str(self.action),
            'motif=' + self.motif.motif_string,
            'postfilter=' + str(self.postfilter)] +
            [x + '=' + str(round(self.score[x], 2)) for x in self.score]
            ) + ')'

class Location(object):
    '''
    Class representing the place a Target aligns on a contig
    '''
    __slots__ = [
        'feature',
        'contig',
        'orientation',
        'start',
        'end',
        'upstream',
        'sequence',
        'downstream',
    ]
    # TODO: Consider removing 'feature' from this class
    def __init__(self, feature, contig, orientation, start, end, upstream, sequence, downstream):
        """
        Create an object that stores a location on a contig
        """
        self.feature = feature
        self.contig = contig
        self.orientation = orientation
        self.start = start
        self.end = end
        self.upstream = upstream
        self.sequence = sequence
        self.downstream = downstream
    
    def __eq__(self, other):
        return (
            (self.feature == other.feature) and
            (self.contig == other.contig) and
            (self.orientation == other.orientation) and
            (self.start == other.start) and
            (self.end == other.end) # and
            #(self.upstream == other.upstream) and
            #(self.sequence == other.sequence) and
            #(self.downstream == other.downstream)
        )
        
    # __lt__(), __gt__(), __le__(), __ne__(), __ge__(), __eq__()
    
    def __repr__(self):
        """
        Return a string containing a printable representation of the Location object.
        """
        return self.__class__.__name__ + '(' + ', '.join([
            'feature={}'.format(self.feature),
            'position={}:{}:{}..{}'.format(self.contig, self.orientation, self.start, self.end),
            ]) + ')'

class Target(object):
    """Data structure defining a gRNA Target"""
    prefix = 'Target'
    sequences = {} # key=nucleotide sequence, value=ExcisionTarget/ReversionTarget object
    indices = {} # key=exTarget-102, value=ExcisionTarget/ReversionTarget object
    #equivalents = {} # key=Target.name, value=set(Target, Target, ...)
    equivalents = {} # key=sequence, value=set(sequence, sequence, ...)
    
    logger = logger.getChild('Target')
    
    # TODO: Strongly consider making 'Target' objects use 'Location' objects instead of just tuples
    # TODO: Strongly consider letting 'Target' objects store/calculate scores for several different sequences
    #       (that are all within some range of homology)
    # TODO: Allow each 'Target' to have a separate off-target score for '8-mer', '12-mer', and 'full-mer'
    def __init__(self, feature, contig, orientation, start, end, upstream, downstream, sequence, side, spacer, pam, motif, parsed_motif, add=True):
        """Create a structure for holding individual sequence information"""
        location = (feature, contig, orientation, start, end, upstream, downstream)
        self.locations = set()

        self.sequence = sequence
        self.side = side
        self.spacer = spacer
        self.pam = pam
        self.motif = motif
        self.parsed_motif = parsed_motif
        
        # Final weight calculation
        self.weight = 0.0
        self.rank = None

        # List to store alignments
        self.alignments = []

        # Get the index number
        self.index = len(self.sequences)
        self.name = self.prefix + '-' + str(self.index)

        if add:
            the_target = self.sequences.setdefault(self.sequence, self)
            the_target.locations.add(location)
            self.indices.setdefault(the_target.name, the_target)

        # Scores for this sequence only (not PairedSequenceAlgorithm)
        self.score = {}

        # We need an off-target score for each algorithm and each potential dDNA
        # Thus it will be structured as a dict nested within a dict:
        # self.off_targets = {'CFD': {'exDonor-0': 10.0, 'exDonor-1': 12.0}}
        self.off_targets = {}
        if (len(self.locations) > 0):
            self.calculate_default_scores()
    
    @classmethod
    def load_spacers(cls, args, filenames):
        '''
        Populate '*Target.sequences' and '*Target.indices' according to input FASTA files
        '''

        # The 'excision-targets.fasta' file has the following format:
        #  >exTarget-17 motif=N{17}|N{3}>NGG locations=2 alignments=2/12 on-target=69.01 off-target=100.0 pam=GGG
        #  CCTCGAGCACTTCCACTGTG
        # The PAM motif should be obtained from either the 'motif=...' attribute,
        # or from the '--motifs ...' command line option

        for filename in filenames:
            with open(filename, 'r') as flo:
                # For now, we will ignore the 'motif=...' attribute
                contigs = utils.old_load_fasta_file(filename)

                for header, sequence in contigs.items():
                    # Search 'sequence' for any number of SPACER>PAM motif matches
                    # Get the 'spacer', 'side', and 'pam' for each 'motif' for this 'sequence'
                    # Pass through 'prefilter'
                    targets = cls.get_targets(args, sequence)

                    for t in targets:
                        # Create the *Target object, and add to *Target.sequences and *Target.indices
                        cls(
                            feature="Feature",
                            contig=header,
                            orientation=t[0],
                            start=t[1],
                            end=t[2], # should be equal to len(sequence)
                            upstream=t[3],
                            downstream=t[4],
                            squence=t[5],
                            side=t[6],
                            spacer=t[7],
                            pam=t[8],
                            motif=t[9], # also called 'motif_string'
                            parsed_motif=t[10]
                        )
        return None

    @classmethod
    def load_alignment(cls, filename, args, contigs):
        """
        Read in alinment file (SAM/BLASTN).
        sep is the separator for the header. Positions are converted to 0-index
        Creates a list of Sequence objects
        """
        # Iterate through all alignment records
        for record in args.selected_aligner.load(filename):
            # Record(
            #     query_name, subject_name,
            #     query_sequence, subject_sequence,
            #     query_position, subject_position,
            #     query_length, subject_length,
            #     flags, cigar, score, evalue, length
            # )
            
            target = cls.indices[record.query_name] # key=exTarget-519
            
            # Get orientation
            alignment_orientation = utils.sam_orientation(record.flags)
            
            # Get alignment position
            alignment_contig = record.subject_name
            alignment_start, alignment_end = record.subject_position
            
            # Reverse-complement if needed
            alignment_contig_sequence = contigs[alignment_contig]
            alignment_sequence = alignment_contig_sequence[alignment_start:alignment_end]
            alignment_upstream = alignment_contig_sequence[alignment_start-10:alignment_start]
            alignment_downstream = alignment_contig_sequence[alignment_end:alignment_end+10]
            #actual_sequence = record.query_sequence # Should be column: 9, SEQ, query SEQuence on the same strand as the reference
            if (alignment_orientation == '-'):
                alignment_sequence = nucleotides.rc(alignment_sequence)
                alignment_upstream, alignment_downstream = nucleotides.rc(alignment_downstream), nucleotides.rc(alignment_upstream)
                #actual_sequence = nucleotides.rc(actual_sequence)
            
            # No actual check made here to see if length of query matches length of subject, and they still conform to motif.
            # This is done inside the 'target.add_alignment()' method
            
            target.add_alignment(
                args,
                alignment_sequence, # aligned_sequence (as when matched with reference, thus may be revcomp of initial query)
                alignment_contig, # aligned_contig
                alignment_start, # aligned_start
                alignment_end, # aligned_end
                alignment_orientation, # aligned_orientation (+/-)
                alignment_upstream,
                alignment_downstream,
            )
        
        cls.logger.info(cls.__name__ + ' alignment file parsed: {!r}'.format(filename))
    
    @classmethod
    def generate_query_fasta(cls, filename, sep=':'):
        """
        Creates a FASTA file (That becomes header for SAM file)
        >exTarget-id feature:contig:orientation:start..end ...
        """
        with open(filename, 'w') as flo:
            for sequence, obj in sorted(cls.sequences.items(), key=lambda x: int(x[1].name.split('-')[1])):
                #for location in obj.locations:
                #    print(' '.join(['>'+obj.format_location(location), obj.name]), file=flo)
                #    print(sequence, file=flo)
                
                # A bug in Bowtie2 makes it so that if a FASTA hader has this:
                #   >contig motif=N{20}>NGG
                # Then the query name will be 'NGG'.
                # Because it erroneously thinks the header is to the right of the second '>' symbol
                
                #print(' '.join(['>'+obj.name] + sorted([obj.format_location(x, sep) for x in obj.locations]) + ['motif='+obj.motif]), file=flo)
                print(' '.join(['>'+obj.name] + sorted([obj.format_location(x, sep) for x in obj.locations])), file=flo)
                print(sequence, file=flo)
                
        cls.logger.info(cls.__name__ + ' query FASTA generated: {!r}'.format(filename))
        return filename

    @classmethod
    def generate_query_fastq(cls, filename, sep=':', qualscore_char='I'):
        """
        Creates a FASTQ file (That becomes header for a SAM record)
          @exTarget-id feature:contig:orientation:start..end ...
          ...
          +
          ...
        """
        with open(filename, 'w') as flo:
            for sequence, obj in sorted(cls.sequences.items(), key=lambda x: int(x[1].name.split('-')[1])):
                print(' '.join(['@'+obj.name] + sorted([obj.format_location(x, sep) for x in obj.locations])), file=flo)
                print(sequence, file=flo)
                print('+', file=flo)
                print(''.join(qualscore_char for x in sequence))

        cls.logger.info(cls.__name__ + ' query FASTQ generated: {!r}'.format(filename))
        return filename
    
    @classmethod
    def generate_spacers_fasta(cls, filename, sep=':'):
        """
        Creates a FASTA file
        >exTarget-id motif=STR locations=N alignments=N/N on-target=N off-target=N
        """
        with open(filename, 'w') as flo:
            for sequence, obj in cls.sequences.items():
                print(' '.join([
                    '>' + obj.name,
                    'motif=' + obj.motif,
                    'locations=' + str(len(obj.locations)),
                    'alignments=' + str(len([a for a in obj.alignments if a.postfilter])) + '/' + str(len(obj.alignments)),
                    'on-target=' + str(round(obj.score['Azimuth'], 2)),
                    'off-target=' + str(round(obj.off_targets['Hsu-Zhang'], 2)), # TODO: Change this to be only the Algorithms specified by the user on the command line
                    'weight=' + str(obj.weight),
                    'rank=' + str(obj.rank),
                    'pam=' + obj.pam,
                ]), file=flo)
                print(obj.spacer, file=flo)
        cls.logger.info(cls.__name__ + ' spacers FASTA generated: {!r}'.format(filename))
        return filename
    
    def format_location(self, location, sep=':'):
        feature, contig, orientation, start, end, upstream, downstream = location
        return sep.join([feature, contig, orientation, '..'.join([str(start), str(end)])])
    
    def format_sequence(self):
        """
        Returns nucleotides with proper 'SPACER>PAM' or 'PAM<SPACER' formatting.
        """
        if (self.side == '>'):
            return self.spacer + '>' + self.pam
        else:
            return self.pam + '<' + self.spacer
    
    @classmethod
    def target_filter(cls, args, sequence, side, target, pam, upstream, downstream):
        '''
        Filters the candidate gRNA sequence based on the following criteria:
         1) case: ignore, upper-only, lower-only, mixed-lower, mixed-upper, mixed-only
         2) ambiguous character expansion: exclusive, discard, keep, disambiguate
         3) SPACER>PAM check using regex (following disambiguation expansion)
         4) Prefilters: maximum consecutive Ts, %GC
        
        Returns list of validated sequences
        '''
        # TODO: The 'Target.target_filter()' method should also return WHY a potential target was rejected.
        #       Something like the following:
        #         [case=pass, ambiguities=pass, motifs=pass, polyT=fail, %GC=skipped]
        
        seq = sequence
        
        # Check the case of the potential gRNA sequence
        if (args.case == "upper-only"):
            if regex.search(r'[a-z]', seq):
                return [] # Reject this sequence because it has lower-case characters
        elif (args.case == "lower-only"):
            if regex.search(r'[A-Z]', seq):
                return [] # Reject this sequence because it has upper-case characters
        elif (args.case == "mixed-lower"):
            if not regex.search(r'[a-z]', seq):
                return []
        elif (args.case == "mixed-upper"):
            if not regex.search(r'[A-Z]', seq):
                return []
        elif (args.case == "mixed-only"):
            if not (regex.search(r'[a-z]', seq) and regex.search(r'[A-Z]', seq)):
                return [] # Reject this sequence because it does not have both lower-case and upper-case characters
        #elif (args.case == "ignore") # then do nothing
        #    pass
        
        # Molecule  5'-sequence-3'                                          Description     ignore  upper-only  lower-only  mixed-lower  mixed-upper  mixed-only
        # ========  ======================================================  ==============  ======  ==========  ==========  ===========  ===========  ==========
        # genome    ACCATAGGAATCCAGCGGCGATCTTAAaggaggatctaggtcgatagcggaata  -               -       -           -           -            -            -
        # spacer           GAATCCAGCGGCGATCTTAA                             All upper-case  keep    keep        discard     discard      keep         discard
        # spacer                     GCGATCTTAAaggaggatct                   Mixed case      keep    discard     discard     keep         keep         keep
        # spacer                               aggaggatctaggtcgatag         All lower-case  keep    discard     keep        keep         discard      discard
        
        # Discard potential gRNAs that have mismatches with their target site
        #if (args.case == "invariant-lower"):
        #    pass
        #elif (args.case == "invariant-upper"):
        #    pass
        
        # Convert input sequence to upper-case so it can be evaluated by the scoring algorithms
        seq = seq.upper()
        
        
        if (args.ambiguities == 'discard'):
            # If target sequence has any ambiguities, then discard it
            if regex.search(r'[^ATCGatcg]', seq):
                return []
            # Otherwise, proceed
            seqs1 = [seq]
        elif (args.ambiguities == 'disambiguate'):
            # Disambiguate sequences if necessary
            seqs1 = nucleotides.disambiguate_iupac(seq)
        elif (args.ambiguities == 'exclusive'):
            # If no ambiguous characters are found, hten discard it
            if not regex.search(r'[^ATCGatcg]', seq):
                return [] # Reject this sequence
            # Otherwise, disambiguate
            seqs1 = nucleotides.disambiguate_iupac(seq)
        else:
            # Do nothing if just 'keep'
            seqs1 = [seq]
        
        
        # After disambiguation, some spacers may no longer confine to their motifs.
        # Thus, we must remove targets that do not confine to at least one of
        # the defined SPACER>PAM motifs.
        #
        # For instance, if motif is AN{4}>NGG
        #     genome ACCTACATCWAGCTAGGCTCTAA
        #     spacer          WAGCT
        #        pam               AGG
        # expansion1          TAGCT           <-- Discard
        # expansion2          AAGCT           <-- Keep
        # 
        # Also, separate the SPACER and PAM motifs
        seqs2 = []
        sides2 = []
        targets2 = []
        pams2 = []
        for nt in seqs1:
            # for i in range(len(args.parsed_motifs)):
            for motif in OnTargetMotif.motifs:
            #for spacers, pams, side in args.parsed_motifs:
                #m = nucleotides.motif_conformation(nt, spacers, pams, side)
                #spacers, pams, side = args.parsed_motifs[i]
                spacers, pams, side = motif.parsed_list
                #compiled_regex = args.compiled_motifs[i]
                m = nucleotides.motif_conformation2(nt, side, motif.compiled_regex)
                if m:
                    seqs2.append(nt)
                    sides2.append(side) # TODO: write debug of this to test it is correct
                    targets2.append(m[0])
                    pams2.append(m[1])
                    #break # consider breaking here
        #seqs = temp_seqs
        
        # Remove targets with T{5,} (old code)
        #seqs = [ nt for nt in seqs if ('T'*(args.max_consecutive_ts+1) not in nt) ]
        
        # Apply all prefilters to each sequence
        seqs3 = []
        for i in range(len(seqs2)):
            prefilter_passes = []
            prefilter_choices = [C for C in algorithms.single_algorithms + algorithms.paired_algorithms + algorithms.batched_single_algorithms if C.prefilter] # Really, should only be SingleSequenceAlgorithms
            this = (seqs2[i], sides2[i], targets2[i], pams2[i], upstream, downstream) # For GC and PolyT, only 'target' is used
            for C in prefilter_choices:
                c_score = C. calculate(this)
                if (C.minimum <= c_score <= C.maximum):
                    prefilter_passes.append(True)
                else:
                    prefilter_passes.append(False)
            # Code to handle PairedSequenceAlgorithms
            # parent = (sequence, target, pam, upstream, downstream)
            #for C in algorithms.paired_algorithms:
            #    c_score = C.calculate(parent, this)
            #    if C.prefilter:
            #        if (C.minimum <= c_score <= C.maximum):
            #            prefilter_passes.append(True)
            #        else:
            #            prefilter_passes.append(False)
            
            # If alignment meets all prefilter criteria, then set it as True
            if all(prefilter_passes): # FYI: all([]) returns True
                #seqs3.append(seqs2[i])
                seqs3.append((seqs2[i], targets2[i], pams2[i]))
                
        
        
        
        # # Remove targets whose %GC is outside the chosen bounds
        # temp_seqs3 = []
        # temp_targets3 = []
        # temp_pams3 = []
        # for i in range(len(temp_seqs2)):
        #     if (args.target_gc[0] <= nucleotides.gc_score(temp_targets2[i]) <= args.target_gc[1]):
        #         temp_seqs3.append(temp_seqs2[i])
        #         temp_targets3.append(temp_targets2[i])
        #         temp_pams3.append(temp_pams2[i])
        # #seqs = temp_seqs2
        
        # [(seq, spacer, pam), (seq, spacer, pam)]
        #rets = []
        #for i in range(len(temp_seqs2)):
        #    rets.append((temp_seqs2[i], temp_targets2[i], temp_pams2[i]))
        #return rets
        return seqs3
    
    @classmethod
    def old_get_targets(cls, args, sequence):
        """
        Tries to match all OnTargetMotif motifs to the input sequence
        ANYWHERE in the sequence
        Returns list of matches (as a tuple) from all OnTargetMotif motifs.
        Does NOT evaluate how good each Target is (in terms of score).
        """
        targets = set()
        #for seq_i, sequence in enumerate(dDNAs):
        for orientation in ['+', '-']:
            if (orientation == '-'):
                sequence = nucleotides.rc(sequence)
            
            #for i in range(len(args.parsed_motifs)):
            for mymotif in OnTargetMotif.motifs: # <------ Do I need to add OffTargetMotif.motifs here as well?
                #spacers, pams, side = args.parsed_motifs[i]
                spacers, pams, side = mymotif.parsed_list
                #compiled_regex = args.compiled_motifs[i]
                #matches = nucleotides.motif_search(sequence, spacers, pams, side)
                matches = nucleotides.motif_search2(sequence, side, mymotif.compiled_regex)
                for seq, start, end, spacer, pam in matches:
                    if (orientation == '-'):
                        start, end = len(sequence) - end, len(sequence) - start
                    upstream = sequence[start-10:start]
                    downstream = sequence[end:end+10]
                    filtered_targets = cls.target_filter(seq, spacer, pam, upstream, downstream, args)
                    for filt_seq, filt_spacer, filt_pam in filtered_targets:
                        targets.add((orientation, start, end, upstream, downstream, filt_seq, side, filt_spacer, filt_pam, mymotif.motif_string, tuple([tuple(x) if isinstance(x, list) else x for x in mymotif.parsed_list])))
        return sorted(targets) # becomes a list
    
    @classmethod
    def get_targets(cls, args, sequence, start=0, end=None, protruded_targets=False):
        """
        Tries to match all OnTargetMotif motifs to the input sequence
        ANYWHERE in the sequence
        Returns list of matches (as a tuple) from all OnTargetMotif motifs.
        Does NOT evaluate how good each Target is (in terms of score).
        By default 'protruded_targets=False', this function will not find Targets
        that are partially outside the input sequence bounds.
        When 'progruded_targets=True', then Targets that only partially overlap
        the input Feature will be identified.
        """
        # TODO: Add argument 'search=True' that when false, will disable searching for multiple Targets.
        #       Thus, 'start=0' and 'end=len(sequence)'
        
        if (end == None):
            end = len(sequence)
        
        targets = set()
        
        for mymotif in OnTargetMotif.motifs:
            spacers, pams, side = mymotif.parsed_list
            if protruded_targets:
                max_target_length = max(map(len, spacers)) + max(map(len, pams))
                actual_start = start - (max_target_length-1)
                actual_end = end + (max_target_length-1)
            else:
                actual_start = start
                actual_end = end
            
            orientation = '+'
            search_seq = sequence[actual_start:actual_end]
            # Perform the string search
            matches = nucleotides.motif_search2(search_seq, side, mymotif.compiled_regex)
            
            for seq, m_start, m_end, spacer, pam in matches:
                s_start = actual_start + m_start
                s_end = actual_start + m_end
            
                upstream = sequence[s_start-10:s_start]
                downstream = sequence[s_end:s_end+10]
                
                filtered_targets = cls.target_filter(args, seq, side, spacer, pam, upstream, downstream)
                for filt_seq, filt_spacer, filt_pam in filtered_targets:
                    targets.add((orientation, s_start, s_end, upstream, downstream, filt_seq, side, filt_spacer, filt_pam, mymotif.motif_string, tuple([tuple(x) if isinstance(x, list) else x for x in mymotif.parsed_list])))
                    #cls.logger.info("DEBUG: {}, {}, {}, {}, {}, {}, {}, {}, {}".format(orientation, seq, mstart, mend, spacer, pam, mymotif.compiled_regex, sstart, send))
            
            orientation = '-'
            search_seq = nucleotides.rc(search_seq)
            # Perform the string search
            matches = nucleotides.motif_search2(search_seq, side, mymotif.compiled_regex)
            
            for seq, m_start, m_end, spacer, pam in matches:
                s_start = actual_end - m_end
                s_end = actual_end - m_start
                
                upstream = nucleotides.rc(sequence[s_end:s_end+10])
                downstream = nucleotides.rc(sequence[s_start-10:s_start])
                
                filtered_targets = cls.target_filter(args, seq, side, spacer, pam, upstream, downstream)
                for filt_seq, filt_spacer, filt_pam in filtered_targets:
                    targets.add((orientation, s_start, s_end, upstream, downstream, filt_seq, side, filt_spacer, filt_pam, mymotif.motif_string, tuple([tuple(x) if isinstance(x, list) else x for x in mymotif.parsed_list])))
                    #cls.logger.info("DEBUG: {}, {}, {}, {}, {}, {}, {}, {}, {}".format(orientation, seq, mstart, mend, spacer, pam, mymotif.compiled_regex, sstart, send))
        
        return sorted(targets) # becomes a list
    
    def calculate_default_scores(self):
        """Populate the default scores for this Sequence"""
        # loc = (feature, contig, orientation, start, end, upstream, downstream)
        loc = next(iter(self.locations)) # Pull an arbitrary location record # TODO: Make this more defined (deterministic)
        #parent = (self.sequence, self.target, self.pam, self.upstream, self.downstream)
        parent = (self.sequence, self.side, self.spacer, self.pam, loc[5], loc[6])
        
        # 'loc[0]' may be a regular 'Feature', or it may be a comma-separated list of features 'FeatureA,FeatureB'
        # At this time, there is no real way to know which Feature the position stored in 'loc[3]' refers to
        # so I will arbitrarily pick the first
        # This is BAD because it means Feature IDs cannot natively contain commas...
        # Apparently this only happens in ReversionTarget objects?
        
        # Code to add position within Feature to Algorithm input
        from . import feature
        f = feature.Feature.features[loc[0].split(',')[0]] # Arbitrarily select first element
        # Assuming that (f.contig == aligend_contig)
        # because this uses the Target.start and Feature.start, this metric is asymmetrical with regards to Feature.orientation
        f_position, f_length, f_orientation = loc[3]-f.start, f.end-f.start, f.strand
        
        for C in algorithms.single_algorithms:
            if (C.default != None):
                self.score[C.name] = C.default
            else:
                self.score[C.name] = C.calculate(parent, parsed_motif=self.parsed_motif, feature_position=f_position, feature_length=f_length, feature_orientation=f_orientation)
        for C in algorithms.paired_algorithms:
            if (C.default != None):
                self.score[C.name] = C.default
            else:
                self.score[C.name] = 0.0
    
    def add_alignment(self, args, aligned_sequence, aligned_contig, aligned_start, aligned_end, aligned_orientation, aligned_upstream, aligned_downstream):
        """Add a genomic position to the list of alignments"""
        # loc = (feature, contig, orientation, start, end, upstream, downstream)
        loc = next(iter(self.locations)) # Pull an arbitrary location record # TODO: Make this more defined (deterministic)
        #parent = (self.sequence, self.target, self.pam, self.upstream, self.downstream)
        parent = (self.sequence, self.side, self.spacer, self.pam, loc[5], loc[6])
        
        # Code to add position within Feature to Algorithm input
        from . import feature
        f = feature.Feature.features[loc[0].split(',')[0]] # Arbitrarily select first element
        # Assuming that (f.contig == aligend_contig)
        # because this uses the Target.start and Feature.start, this metric is asymmetrical with regards to Feature.orientation
        f_position, f_length, f_orientation = aligned_start-f.start, f.end-f.start, f.strand
        
        #aligned_target, aligned_pam, aligned_motif = self.split_spacer_pam(aligned_sequence, args)
        aligned_target, aligned_pam, aligned_motif = self.match_best_motif(aligned_sequence, args) # New code
        
        #aligned_parsed_motif = args.parsed_motifs[args.motifs.index(aligned_motif)]
        #aligned_parsed_motif = aligned_motif.parsed_list
        if ((aligned_target != None) and (aligned_pam != None)):
            a = Alignment(
                aligned_sequence,
                aligned_target,
                aligned_pam,
                aligned_motif,
                aligned_contig,
                aligned_start,
                aligned_end,
                aligned_orientation,
                upstream=aligned_upstream,
                downstream=aligned_downstream
            )
            a.calculate_scores(parent, feature_position=f_position, feature_length=f_length, feature_orientation=f_orientation)
            self.alignments.append(a)
        else:
            print('Cannot add alignment:', aligned_sequence, aligned_contig, aligned_start, aligned_end, aligned_orientation, file=sys.stderr)
    
    def split_spacer_pam(self, sequence, args):
        # TODO: See below FIXME
        # FIXME: the 'PAM-Identity' messes up when both '--motifs' and '--off_target_motifs' are specified
        #        For instance, an alignment that is GGG is sometimes evaluated against 'NAG' instead of 'NGG'
        #        This leads to a PAM-identity of '54.55' instead of '100.00'
        #        WHY is it evaluating against the off-target-motif and not the on-target motif?
        for motif in OnTargetMotif.motifs:
            spacers, pams, side = motif.parsed_list
            m = nucleotides.motif_conformation2(sequence, side, motif.compiled_regex)
            if m: # If the motif matches, then return it immediately
                return m[0], m[1], motif # spacer, pam, motif
        
        for motif in OffTargetMotif.motifs:
            spacers, pams, side = motif.parsed_list
            m = nucleotides.motif_conformation2(sequence, side, motif.compiled_regex)
            if m: # If the motif matches, then return it immediately
                return m[0], m[1], motif
    
        # If the motif does not match, then fudge it for an off-target/on-target motif
        # This will return the last-fudged motif, and not necessarily the best-fudged motif
        # *MAY* need to improve this code later
        for motif in (OnTargetMotif.motifs + OffTargetMotif.motifs)[::-1]:
            spacers, pams, side = motif.parsed_list
            if (side == '>'):
                l = max(map(len, pams))
                return sequence[:-l], sequence[-l:], motif # spacer, pam, motif
            else: #elif (side == '<'):
                l = max(map(len, pams))
                return sequence[:l], sequence[l:], motif # spacer, pam, motif
    
    def match_best_motif(self, sequence, args):
        '''
        Will try finding the motif that matches BEST with the input sequence. 
        '''
        motif_list = OnTargetMotif.motifs + OffTargetMotif.motifs
        spacer_list = []
        pam_list = []
        fuzzy_counts_list = []
        for motif in motif_list:
            #spacers, pams, side = motif.parsed_list
            spacer, pam, fuzzy_counts = motif.find_best_fit(sequence)
            spacer_list.append(spacer)
            pam_list.append(pam)
            fuzzy_counts_list.append(fuzzy_counts)
        
        # Sort by (first) sum, then deletions, insertions, and substitutions (last)
        order = sorted(zip(range(len(fuzzy_counts_list)), fuzzy_counts_list), key=lambda x: (sum(x[1]),)+x[1][::-1])
        
        # Select the best motif that has the LEAST differences (smallest fuzzy_counts)
        best_i = order[0][0]
        return spacer_list[best_i], pam_list[best_i], motif_list[best_i]
    
    def get_features(self):
        """Return (sorted) list of all feature names this Target maps to"""
        return sorted(set(x[0] for x in self.locations))
    
    def get_location_features(self):
        """Return (sorted) list of all feature names this Target maps to"""
        return sorted(set(x[0] for x in self.locations))
    
    def get_contigs(self):
        """Return (sorted) list of all contig names this Target maps to"""
        return sorted(set(x[1] for x in self.locations))
    
### Unused 04/08/2019 ###
#    def get_parent(self):
#        """Returns a parent tuple for an arbitrary location record"""
#        loc = next(iter(self.locations)) # Pull an arbitrary location record
#        #parent = (self.sequence, self.target, self.pam, self.upstream, self.downstream)
#        parent = (self.sequence, self.spacer, self.pam, loc[5], loc[6])
#        return parent
#########################
    
    @classmethod
    def score_batch(cls, args):
        """Performs BatchedSingleSequenceAlgorithm calculations on all sequences"""
        # Get immutable list of dict keys
        t_sorted = sorted(cls.indices.keys(), key=lambda x: int(x.split('-')[1]))
        
        # Make list to populate with values to pass into the calculate() methods
        queries = []
        for i, t_index in enumerate(t_sorted):
            t_obj = cls.indices[t_index]
            loc = next(iter(t_obj.locations)) # Pull an arbitrary location record
            parent = (t_obj.sequence, t_obj.side, t_obj.spacer, t_obj.pam, loc[5], loc[6]) # upstream, downstream
            # TODO: Add 'feature_position', 'feature_length', and 'feature_orientaiton'
            queries.append(parent)
        
        # Loop through the Algorithms
        for C in algorithms.batched_single_algorithms:
            if C.name in args.ontargetfilters:
                # Calculate
                batch_scores = C.calculate(queries) # motif=cls.motif
                
                # Assign the score to the appropriate Target
                for i, t_index in enumerate(t_sorted):
                    t_obj = cls.indices[t_index]
                    t_obj.score[C.name] = batch_scores[i]
    
    # Off-target score
    #  Hsu specificity score
    #   The higher the specificity score, the lower are off-target effects in the genome.
    #   The specificity score ranges from 0-100 and measures the uniqueness
    #   of a guide in the genome. See Hsu et al. Nat Biotech 2013.
    #   (http://dx.doi.org/10.1038/nbt.2647) We recommend values >50, where possible.
    def off_target_score(self, off_target_scores, on_target_scores=(100,)):
        """Calculate the off-target score using the Zhong lab MIT Guide Score algorithm
        
        'off_target_scores' should be a list/tuple, with each element ranging from 0-100.
        'on_target_scores' should be a list/tuple, with each element ranging from 0-100.
        Returns off target score ranging from 0-100.
        
        Returns the naive proportion of gRNAs that will bind to a single target
        out of all possible targets.
        
        Higher values indicate greater specificity to intended target. Higher
        scores are better. Lower values mean that the candidate gRNA may bind
        nonspecifically. Guides with scores over 50 are good enough to be candidates
        (about 50% of the gRNA binding events will be to the intended target).
        The off-target score tells you the inverse probability of Cas9 off-target
        binding. A higher score means the sequence has less chance to bind to
        off-target sequences in the rest of the genome.
        
        The off-target score is from Hsu et al. (2013) doi:10.1038/nbt.2647
        and measures how specific the guide is to the target location. Scores should
        be used to rank guides relative to each other.
        
        This "MIT Guide Score" is defined on http://crispr.mit.edu/about
        Also called "Efficiency score", I think (or Eff Score) (see CRISPOR).
        """
        return  100.0*sum(on_target_scores)/(sum(on_target_scores)+sum(off_target_scores))
    
    def score_off_targets(self, args, homologs):
        """
        Calculate Guide Score (off-target score) for all single/paired
        algorithms with the 'off_target=True' attribute.
        """
        # Import included AddTag-specific modules
        from . import feature
        from .donors import ExcisionDonor, ReversionDonor
        
        # Get list of algorithms whose scores should be used for off-target
        # calculations
        calculators = []
        for C in algorithms.single_algorithms:
            if C.off_target:
                calculators.append(C)
        for C in algorithms.paired_algorithms:
            if C.off_target:
                calculators.append(C)
        
        # Make empty lists for scores
        # It should have this format
        #  dict[algorithm][genome/exDonor/reDonor] = [score1, score2, score3, ...]
        # This is for keeping track of which scores should get the ratio multiplier
        on_targets = {}
        off_targets = {}
        for C in calculators:
            on_targets[C.name] = {'gDNA': [], 'dDNA': []}
            off_targets[C.name] = {'gDNA': [], 'dDNA': []}
        
        # Get list of on-target features
        on_target_features = set([x[0] for x in self.locations])
        if isinstance(self, ReversionTarget):
            on_target_features = set(utils.flatten([x.split(',') for x in on_target_features]))
        
        if homologs:
            for f in list(on_target_features):
                on_target_features.update(homologs.get(f, set()))
        
        self.logger.info('on_target_features' + str(on_target_features)) # FIXME: This info output is messed up in log file. Find out why.
        # Check each alignment
        for a in self.alignments:
            a_features = None
            temp_string = ''
            if a.postfilter:
                if isinstance(self, ExcisionTarget):
                    if a.contig.startswith('exDonor-'):
                        a_features = ExcisionDonor.indices[a.contig].get_features()
                        temp_string += 'case exT-exD'.rjust(18) + ' '
                        
                        for i, C in enumerate(calculators):
                            c_score = a.score[C.name]
                            if (C.minimum <= c_score <= C.maximum):
                                if ((len(a_features) == 0) or (len(on_target_features.intersection(a_features)) == 0)):
                                    # if alignment is NOT a target:
                                    # And is also an exDonor, then we will ignore it
                                    # As we are assuming mono-plex gRNA
                                    a.action = 'skip'
                                # If the features match, then this is an off-target
                                else:
                                    # Scale this by the ratio, as they are dDNA alignments
                                    off_targets[C.name]['dDNA'].append(c_score)
                                    a.action = 'off'
                                    
                    else: # This is a genomic contig and not a dDNA
                        a_features = feature.Feature.get_overlapping_features(a.contig, a.start, a.end)
                        #a_features = self.get_features(features, a.contig, a.start, a.end) # contig, start, end
                        temp_string += 'case exT-gDNA'.rjust(18) + ' '
                        
                        for i, C in enumerate(calculators):
                            c_score = a.score[C.name]
                            if (C.minimum <= c_score <= C.maximum):
                                # if alignment is NOT a target:
                                if ((len(a_features) == 0) or (len(on_target_features.intersection(a_features)) == 0)):
                                    off_targets[C.name]['gDNA'].append(c_score)
                                    a.action = 'off'
                                #if (len(on_target_features.intersection(a_features)) > 0):
                                # if alignment is an intended target
                                else:
                                    on_targets[C.name]['gDNA'].append(c_score)
                                    a.action = 'on'
                
                elif isinstance(self, ReversionTarget):
                    if a.contig.startswith('exDonor-'):
                        # reTarget vs exDonor:
                        #   if multiplex:
                        #     if exDonor == intended exDonor (for this feature), then on-target
                        #     if exDonor is the same feature, but NOT intended, then ignore
                        #     if exDonor is for another feature, then it is off-target
                        #   if monoplex:
                        #     if exDonor == intended exDonor (for this feature), then on-target
                        #     if exDonor is the same feature, but NOT intended, then ignore
                        #     if exDonor is for another feature, then ignore
                        #a_features = ExcisionDonor.indices[a.contig].get_features()
                        intended_contigs = self.get_contigs()
                        temp_string += 'case reT-exD'.rjust(18) + ' '
                        
                        for i, C in enumerate(calculators):
                            c_score = a.score[C.name]
                            if (C.minimum <= c_score <= C.maximum):
                                if a.contig in intended_contigs:
                                    on_targets[C.name]['gDNA'].append(c_score)
                                    a.action = 'on'
                                else:
                                    a.action = 'skip'
                                
                                #if ((len(a_features) == 0) or (len(on_target_features.intersection(a_features)) == 0)):
                                #    pass # do nothing for monoplex
                                #else:
                                #    on_targets[C.name]['gDNA'].append(c_score)
                    
                    elif a.contig.startswith('reDonor-'):
                        # reTarget vs reDonor:
                        #   multiplex: these are all off-target (use ratio)
                        #   monoplex: if reDonor has same feature as reTarget, then it is an off-target (use ratio)
                        a_features = ReversionDonor.indices[a.contig].get_features()
                        temp_string += 'case reT-reD'.rjust(18) + ' '
                        
                        for i, C in enumerate(calculators):
                            c_score = a.score[C.name]
                            if (C.minimum <= c_score <= C.maximum):
                                if ((len(a_features) == 0) or (len(on_target_features.intersection(a_features)) == 0)):
                                    a.action = 'skip'
                                else:
                                    # use ratio
                                    off_targets[C.name]['dDNA'].append(c_score)
                                    a.action = 'off'
                    else:
                        # reTarget vs genome: these are all off-target, regardless of the feature
                        temp_string += 'case reT-gDNA'.rjust(18) + ' '
                    
                        for i, C in enumerate(calculators):
                            c_score = a.score[C.name]
                            if (C.minimum <= c_score <= C.maximum):
                                off_targets[C.name]['gDNA'].append(c_score)
                                a.action = 'off'
            else:
                temp_string += 'failed post-filter' + ' '
            temp_string += ' '.join([a.action.rjust(4), str(a), str(a_features)])
            self.logger.info(temp_string)
        
        for i, C in enumerate(calculators):
            on_str = str(len(on_targets[C.name]['gDNA'])) + ' + (' + str(args.dDNA_gDNA_ratio)+')'+str(len(on_targets[C.name]['dDNA']))
            off_str = str(len(off_targets[C.name]['gDNA'])) + ' + (' + str(args.dDNA_gDNA_ratio)+')'+str(len(off_targets[C.name]['dDNA']))
            self.logger.info(' '.join([self.name, C.name, '('+on_str+')/(('+on_str+') + ('+off_str+'))']))
        
        # Perform off-target calculations
        for C in calculators:
            try:
                on_list = on_targets[C.name]['gDNA'] + args.dDNA_gDNA_ratio * on_targets[C.name]['dDNA']
                off_list = off_targets[C.name]['gDNA'] + args.dDNA_gDNA_ratio * off_targets[C.name]['dDNA']
                self.off_targets[C.name] = self.off_target_score(off_list, on_list)
            except ZeroDivisionError:
                self.off_targets[C.name] = 0.0 # This should never happen
    
    def __repr__(self):
        """Return a string containing a printable representation of the Target object."""
        return self.__class__.__name__ + '(' + ' '.join([
            self.name,
            self.format_sequence(),
            'motif=' + self.motif,
            'locations=' + str(len(self.locations)),
            'alignments=' + str(len([a for a in self.alignments if a.postfilter])) + '/' + str(len(self.alignments))] +
            [x + '=' + str(round(self.score[x], 2)) for x in self.score] + 
            ['OT:' + x + '=' + str(round(self.off_targets[x], 2)) for x in self.off_targets]
            ) + ')'

class ExcisionTarget(Target):
    prefix = 'exTarget'
    # Non-redundant dicts
    sequences = {} # key=nucleotide sequence, value=ExcisionTarget object
    indices = {} # key=exTarget-102, value=ExcisionTarget object
    #equivalents = {} # key=Target.name, value=set(Target, Target, ...)
    equivalents = {} # key=sequence, value=set(sequence, sequence, ...)
    
    logger = logger.getChild('ExcisionTarget')
    
    @classmethod
    def old_search_all_features(cls, args, features, contigs):
        for feature_name in features:
            feature_contig, feature_start, feature_end, feature_strand = features[feature_name]
            
            if feature_contig in contigs:
                cls.expand_feature(args, (feature_name, feature_contig, feature_start, feature_end, feature_strand), (feature_contig, contigs[feature_contig]))
            else:
                cls.logger.info("The contig '{}' for feature '{}' is not in the input FASTA.".format(feature_contig, feature_name))
        
    
    @classmethod
    def search_all_features(cls, args, contigs):
        # Import included AddTag-specific modules
        from . import feature
        
        num_before = len(cls.sequences)
        cls.logger.info('{} starts with {} Target objects'.format(cls.__name__, num_before))
        
        for feature_name, f in feature.Feature.features.items():
            cls.logger.info("Searching Feature '{}' for ExcisionTarget objects.".format(feature_name))
            
            # Search for targets in the feature
            contig_sequence = contigs[f.contig]
            #feature_sequence = contig_sequence[f.start:f.end]
            targets = cls.get_targets(args, contig_sequence, start=f.start, end=f.end, protruded_targets=args.protruded_targets) # Does both orientations (+/-)
            cls.logger.info("  found {} potential Targets".format(len(targets)))
            
            # Create ExcisionTarget objects for each found target
            for t in targets:
                t_orientation = t[0]
                t_start = t[1]
                t_end = t[2]
                t_upstream = t[3]
                t_downstream = t[4]
                t_sequence = t[5]
                t_side = t[6]
                t_spacer = t[7]
                t_pam = t[8]
                t_motif_string = t[9]
                t_motif_parsed_list = t[10]
                
                # These Target coordinate offsets are already converted in the 'Target.get_targets()' function
                #if (t_orientation == '+'):
                #    real_start = f.start + t_start
                #    real_end = f.start + t_end
                #else:
                #    contig_length = len(contig_sequence)
                #    real_start = contig_length - t_end + f.start
                #    real_end = contig_length - t_start + f.start
                
                ##OLD cls(feature, feature_contig, orientation, real_start, real_end, t_upstream, t_downstream, filt_seq, side, filt_spacer, filt_pam, motif.motif_string, motif.parsed_list)
                ##REFERENCE CODE targets.add((orientation, start, end, upstream, downstream, filt_seq, side, filt_spacer, filt_pam, mymotif.motif_string, tuple([tuple(x) if isinstance(x, list) else x for x in mymotif.parsed_list])))
                #cls(f.name, f.contig, t_orientation, real_start, real_end, t_upstream, t_downstream, t_sequence, t_side, t_spacer, t_pam, t_motif_string, t_motif_parsed_list)
                cls(f.name, f.contig, t_orientation, t_start, t_end, t_upstream, t_downstream, t_sequence, t_side, t_spacer, t_pam, t_motif_string, t_motif_parsed_list)
                ## t = (orientation, start, end, t_upstream, t_downstream, filt_seq, side, filt_spacer, filt_pam, args.motifs[i], args.parsed_motifs[i])
                ##SIMILAR ReversionTarget(obj_features, obj.name, t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], t[10])
        
        num_after = len(cls.sequences)
        cls.logger.info('{} ends with {} Target objects'.format(cls.__name__, num_after))
        cls.logger.info('A total of {} Target objects have been added'.format(num_after-num_before))
        
        # END 'search_all_features()'
    
    @classmethod
    def old_expand_feature(cls, args, f, contig_sequence, expansion_size=100, minimum_targets_per_feature=10):
        
        feature_start = f.start
        feature_end = f.end
        
        contig_length = len(contig_sequence)
        
        if (feature_end == None):
            feature_end = len(contig_sequence)
        
        targets = []
        while ((len(targets) < minimum_targets_per_feature) and ((feature_start, feature_end) != (0, contig_length))):
            feature_sequence = contig_sequence[feature_start:feature_end]
            targets = Target.get_targets(args, feature_sequence) # Does both orientations (+/-)
            
            feature_start = max(0, feature_start-expansion_size)
            feature_end = min(contig_length, feature_end + expansion_size)
        
        # Create ExcisionTarget objects for each found target
        for t in targets:
            t_orientation = t[0]
            t_start = t[1]
            t_end = t[2]
            t_upstream = t[3]
            t_downstream = t[4]
            t_sequence = t[5]
            t_side = t[6]
            t_spacer = t[7]
            t_pam = t[8]
            t_motif_string = t[9]
            t_motif_parsed_list = t[10]
            
            if (t_orientation == '+'):
                real_start = feature_start + t_start
                real_end = feature_start + t_end
            else:
                real_start = contig_length - t_end + feature_start
                real_end = contig_length - t_start + feature_start
            
            #OLD cls(feature, feature_contig, orientation, real_start, real_end, t_upstream, t_downstream, filt_seq, side, filt_spacer, filt_pam, motif.motif_string, motif.parsed_list)
            #REFERENCE CODE targets.add((orientation, start, end, upstream, downstream, filt_seq, side, filt_spacer, filt_pam, mymotif.motif_string, tuple([tuple(x) if isinstance(x, list) else x for x in mymotif.parsed_list])))
            cls(f.name, f.contig, t_orientation, real_start, real_end, t_upstream, t_downstream, t_sequence, t_side, t_spacer, t_pam, t_motif_string, t_motif_parsed_list)
            # t = (orientation, start, end, t_upstream, t_downstream, filt_seq, side, filt_spacer, filt_pam, args.motifs[i], args.parsed_motifs[i])
            #SIMILAR ReversionTarget(obj_features, obj.name, t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], t[10])
    
    # @classmethod
    # def get_targets(cls, args, contigs, features):
    #     """
    #     Searches within the annotated features on all contigs for targets
    #     that match the args.motifs criteria, then filters them.
    #     
    #     The constructor call adds the valid gRNA sites to 'ExcisionTarget.sequences' and 'ExcisionTarget.indices'.
    #     """
    #     # Find unique gRNA sites within each feature
    #     # Use a sliding window to make a list of queries
    #     for feature in features: # This is the filtered list of features, so we want to find all of these
    #         
    #         feature_contig, feature_start, feature_end, feature_strand = features[feature]
    #         if (feature_end == None):
    #             feature_end = len(contigs[feature_contig])
    #         
    #         # Make sure the contig the feature is on is present in the FASTA
    #         if feature_contig in contigs:
    #             # Find a site within this feature that will serve as a unique gRNA
    #             
    #             # for each orientation:
    #             # if (args.strands in ['+', 'both']):
    #             #  targets.extend...
    #             # if (args.strands in ['-', 'both']):
    #             #  targets.extend...
    #             
    #             # Code to enable automatic feature expansion
    #             # expansion_size = 10 # number of nt to expand features in 5' and 3' direction
    #             # minimum_spacers_per_feature = 100
    #             # found_spacers = 0
    #             # while ((found_spacers < minimum_spacers_per_feature) and ((feature_start, feature_end) != (0, contig_length)):
    #             #     for orientation in ['+', '-']:
    #             #         for motif in OnTargetMotif.motifs:
    #             #             for m in matches:
    #             #                 for target in filtered_targets:
    #             #                     cls() # need to add the novel (feature_start, feature_end) so it can be stored in this object
    #             #                     found_spacers += 1.0/len(filtered_targets) # Penalize for targets at the exact same location (IUPAC disambiguation)
    #             #     # Expand feature before next iteration
    #             #     feature_start = max(0, feature_start-expansion_size)
    #             #     feature_end = min(contig_length, feature_end + expansion_size)
    #             
    #             # Search both the '+' and '-' strands
    #             for orientation in ['+', '-']:
    #                 if (orientation == '+'):
    #                     sequence = contigs[feature_contig][feature_start:feature_end]
    #                 else:
    #                     sequence = nucleotides.rc(contigs[feature_contig][feature_start:feature_end])
    #                 
    #                 #for i in range(len(args.parsed_motifs)):
    #                 for motif in OnTargetMotif.motifs:
    #                     #spacers, pams, side = args.parsed_motifs[i]
    #                     spacers, pams, side = motif.parsed_list
    #                     #compiled_regex = args.compiled_motifs[i]
    #                     #matches = nucleotides.motif_search(sequence, spacers, pams, side)
    #                     matches = nucleotides.motif_search2(sequence, side, motif.compiled_regex)
    #                     for seq, start, end, spacer, pam in matches:
    #                         if (orientation == '+'):
    #                             real_start = feature_start + start
    #                             real_end = feature_start + end
    #                         else:
    #                             real_start = len(sequence) - end + feature_start
    #                             real_end = len(sequence) - start + feature_start
    #                         t_upstream = sequence[start-10:start]
    #                         t_downstream = sequence[end:end+10]
    #                         filtered_targets = target_filter(seq, spacer, pam, t_upstream, t_downstream, args)
    #                         for filt_seq, filt_spacer, filt_pam in filtered_targets:
    #                             cls(feature, feature_contig, orientation, real_start, real_end, t_upstream, t_downstream, filt_seq, side, filt_spacer, filt_pam, motif.motif_string, motif.parsed_list)
    #                             # Maybe add to ReversionDonor here
    #         else:
    #             cls.logger.info("The contig '{}' for feature '{}' is not in the input FASTA.".format(feature_contig, feature))
    
    #def get_donors(self):
    #    return [ReversionDonor.indices[x] for x in self.get_contigs()]

class ReversionTarget(Target):
    prefix = 'reTarget'
    # Non-redundant dicts
    sequences = {} # key=nucleotide sequence, value=ReversionTarget object
    indices = {} # key=reTarget-234, value=ReversionTarget object
    #equivalents = {} # key=Target.name, value=set(Target, Target, ...)
    equivalents = {} # key=sequence, value=set(sequence, sequence, ...)
    
    logger = logger.getChild('ReversionTarget')
    
    @classmethod
    def create_target_objects(cls):
        """
        Extracts the spacer element from each ExcisionDonor object, and uses it
        to create the associated ReversionTarget objects.
        
        The constructor call adds the valid gRNA sites to 'ReversionTarget.sequences' and 'ReversionTarget.indices'.
        """
        from .donors import ExcisionDonor
        
        for dDNA, obj in ExcisionDonor.sequences.items():
            # TODO: WHy did I do this? Joining the feature names like this makes storing a 'location' more difficult
            #       It is causing a bug in April 2020.
            #       Apparently, it is used in 'get_reTarget_homologs()' and other functions (extensively) in '_subroutine_generate_all.py'
            #obj_features = ','.join([x[0] for x in list(obj.locations)])
            obj_features = ','.join(sorted(set([x[0] for x in obj.locations])))
            
            # Populate the ReversionTarget sequences indices dicts
            for t in obj.spacers: # [(orientation, start, end, filt_seq, side, filt_spacer, filt_pam), ...]
                #final_targets.append(tuple([obj.name, obj_feature, obj_contig] + list(t)))
                # t = (orientation, start, end, t_upstream, t_downstream, filt_seq, side, filt_spacer, filt_pam, args.motifs[i], args.parsed_motifs[i])
                ReversionTarget(obj_features, obj.name, t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], t[10])
    
    def get_donors(self):
        from .donors import ExcisionDonor
        return [ExcisionDonor.indices[x] for x in self.get_contigs()]
