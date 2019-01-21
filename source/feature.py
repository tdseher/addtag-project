#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/feature.py

# Import standard packages
import logging

# Import non-standard packages
import regex

# Import included AddTag-specific modules
#from .__init__ import Target
#from . import targets

class Feature(object):
    features = {}
    excluded_features = {}
    # Origin flags
    NONE=0
    INPUT=1
    DERIVED=2
    def __init__(self, contig, start, end, strand, name=None, attributes=None, source=None, feature_type=None, score=None, frame=None, origin=NONE, sep=';', expand_parent=None):
        self.contig = contig
        self.source = source
        self.feature_type = feature_type
        self.start = start
        self.end = end
        if (score != '.'):
            self.score = score
        else:
            self.score = None
        self.strand = strand
        if (frame != '.'):
            self.frame = frame
        else:
            self.frame = None
        self.attributes = {}
        if isinstance(attributes, dict):
            self.attributes = attributes
        elif isinstance(attributes, str): # "ID=12;Parent=nope"
            alist = regex.split(sep+'\s*', attributes) # ['ID=12', 'Parent=nope']
            self.attributes = dict([regex.split('\s*=\s*', x) for x in alist]) # {'ID': '12', 'Parent': 'nope'}
        self.name = name
        self.origin = origin
        self.expand_parent = expand_parent
    
    def get_expand_parent(self):
        if (self.expand_parent == None):
            return self
        else:
            return self.expand_parent.get_expand_parent()
    
    @classmethod
    def get_gene_from_feature(cls, feature_name, feature2gene):
        parent = cls.features[feature_name].get_expand_parent().name
        
        return feature2gene[parent]
    
    @classmethod
    def parse_gff_line(cls, text):
        obj = None
        
        line = text.rstrip()
        if not line.startswith('#'):
            sline = line.split('\t')
            if (len(sline) > 6):
                contig = sline[0]
                source = sline[1]
                feature_type = sline[2]
                start = int(sline[3])-1 # Should always be true: start <= end 
                end = int(sline[4]) # Should always be true: start <= end
                score = sline[5]
                strand = sline[6]
                frame = None
                attributes = None
                if (len(sline) > 7):
                    frame = sline[7]
                if (len(sline) > 8):
                    attributes = sline[8]
                obj = cls(contig, start, end, strand, source=source, feature_type=feature_type, score=score, frame=frame, attributes=attributes, origin=Feature.INPUT)
        return obj
    
    @classmethod
    def load_gff_file(cls, filename, feature_types, excluded_feature_types, selected_features, tag):
        """
        Load General Feature Format (GFF) file into dict()
        One line per feature, each containing 9 columns of data, plus optional
        track definition lines.
        
        Converts positions to 0-based index.
        """
        # Fields must be tab-separated. Also, all but the final field in each
        # feature line must contain a value; "empty" columns should be denoted
        # with a '.'
        #  1) seqid - name of the chromosome or scaffold; chromosome names can
        #     be given with or without the 'chr' prefix. Important note: the
        #     seqname must be a standard chromosome name or an identifier such as
        #     a scaffold ID, without any additional content such as species or
        #     assembly.
        #  2) source - name of the program that generated this feature, or the
        #     data source (database or project name)
        #  3) feature - feature type name, e.g. Gene, Variation, Similarity
        #  4) start - Start position of the feature, with sequence numbering
        #     starting at 1.
        #  5) end - End position of the feature, with sequence numbering
        #     starting at 1.
        #  6) score - A floating point value. As in earlier versions of the format,
        #     the semantics of the score are ill-defined. It is strongly
        #     recommended that E-values be used for sequence similarity features,
        #     and that P-values be used for ab initio gene prediction features.
        #     If there is no score, put a '.' (a period) in this field.
        #  7) strand - defined as '+' (forward), '-' (reverse), '.' (unstranded),
        #     '?' (relevant, but unknown).
        #  8) frame - for CDS features, '0', '1' or '2'. '0' indicates that the first base of
        #     the feature is the first base of a codon, '1' that the second base
        #     is the first base of a codon, and so on. Other features can use '.'.
        #  9) attribute - A semicolon-separated list of tag-value pairs, providing
        #     additional information about each feature. A list of feature
        #     attributes in the format tag=value. Multiple tag=value pairs are
        #     separated by semicolons. URL escaping rules are used for tags or
        #     values containing the following characters: ",=;". Spaces are allowed
        #     in this field, but tabs must be replaced with the %09 URL escape.
        #     This field is not required.
        #     Column 9 tags have predefined meanings:
        #       ID - Indicates the unique identifier of the feature. IDs must be
        #            unique within the scope of the GFF file.
        #       Name - Display name for the feature. This is the name to be displayed to the user. Unlike IDs, there is no requirement that the Name be unique within the file.
        #       Alias - A secondary name for the feature. It is suggested that this tag be used whenever a secondary identifier for the feature is needed, such as locus names and accession numbers. Unlike ID, there is no requirement that Alias be unique within the file.
        #       Parent - Indicates the parent of the feature. A parent ID can be used to group exons into transcripts, transcripts into genes, and so forth. A feature may have multiple parents. Parent can *only* be used to indicate a partof relationship.
        #       Target - Indicates the target of a nucleotide-to-nucleotide or protein-to-nucleotide alignment. The format of the value is "target_id start end [strand]", where strand is optional and may be "+" or "-". If the target_id contains spaces, they must be escaped as hex escape %20.
        #       Gap - The alignment of the feature to the target if the two are not collinear (e.g. contain gaps). The alignment format is taken from the CIGAR format described in the Exonerate documentation. http://cvsweb.sanger.ac.uk/cgi-bin/cvsweb.cgi/exonerate?cvsroot=Ensembl). See the GFF3 specification for more information.
        #       Derives_from - Used to disambiguate the relationship between one feature and another when the relationship is a temporal one rather than a purely structural "part of" one. This is needed for polycistronic genes. See the GFF3 specification for more information.
        #       Note - A free text note.
        #       Dbxref - A database cross reference. See the GFF3 specification for more information.
        #       Ontology_term - A cross reference to an ontology term. See the GFF3 specification for more information.
        #     Multiple attributes of the same type are indicated by separating the
        #     values with the comma "," character, as in: 'Parent=AF2312,AB2812,abc-3'
        #     Note that attribute names are case sensitive. "Parent" is not the
        #     same as "parent". All attributes that begin with an uppercase letter
        #     are reserved for later use. Attributes that begin with a lowercase
        #     letter can be used freely by applications. You can stash any
        #     semi-structured data into the database by using one or more
        #     unreserved (lowercase) tags.
        
        with open(filename, 'r') as flo:
            for line in flo:
                line = line.rstrip()
                obj = cls.parse_gff_line(line)
                if obj:
                    if ((obj.feature_type in feature_types) or ('all' in feature_types)):
                        obj.name = obj.attributes[tag]
                        
                        # Reduce the total number of features to just the ones indicated in the selection
                        if selected_features:
                            if (obj.name in selected_features):
                                cls.features[obj.name] = obj
                        else:
                            cls.features[obj.name] = obj
                    elif ((obj.feature_type in excluded_feature_types) or ('all' in excluded_feature_types)):
                        obj.name = obj.attributes[tag]
                        cls.excluded_features[obj.name] = obj
        
        #logger.info('GFF file parsed: {!r}'.format(filename))
    
    @classmethod
    def assert_features(cls, selection, contigs):
        # Current implementation limitations
        # GFF file only has 'Gene=' tag on ONE of the homologs, and not the other
        # User will have to specify the feature twice if he wants to target both homologs
        
        # Require at least one in order to filter
        # Feature has this format: key=gene/tag, value = (contig, start(bp), end(bp), strand)
        if selection:
            for s in selection:
                f = cls.features.get(s)
                if not f:
                    raise Exception("Feature specified as '--selection "+s+"' does not exist in input '--gff'.")
        
        for feature_name, f in cls.features.items():
            if not f.contig in contigs:
                raise Exception("Feature '{}' lies on contig '{}' which does not exist in input '--fasta'.".format(feature_name, f.contig))
        
        if (len(cls.features) == 0):
            raise Exception("Input '--gff' and '--selection' combination returns no valid features.")
    
    @classmethod
    def expand_all_features(cls, args, contigs):
        checked_keys = []
        # Some make-shift code to get this working when the 'Feature.features' dict changes size during the loop process
        while(len(checked_keys) < len(Feature.features)):
            for k, v in Feature.features.items():
                if (k not in checked_keys):
                    feature_name = k
                    f = v
                    break
            checked_keys.append(feature_name)
            
            contig_sequence = contigs[f.contig]
            if (f.origin == Feature.INPUT):
                #ef = f.expand_feature(args, contig_sequence)
                logging.info('Expanding feature: {}'.format(feature_name))
                f.expand_feature(args, contig_sequence)
    
    def expand_feature(self, args, contig_sequence):
        """
        New method of expanding a feature.
        Will generate all possible expanded features, which will be
        evaluated later to determine which is the best.
        """
        # Import included AddTag-specific modules
        from . import targets
        
        # Option 1 (Can't do because target evaluation happens later)
        #  keep expanding the feature until a target of minimum quality is found
        #  or until the maximum size 'args.feature_expansion_lengths[1]' is reached
        
        # Option 2 (Can't do because target evaluation happens later)
        #  scan the widest region for targets
        #  evaluate target scores/weights
        #  pick the best target
        #  center/justify the FEATURE and TARGET
        #  create the derived feature object
        
        # Option 3
        #  create all derived features within 'args.feature_expansion_lengths' limits
        #  these are centered/justified already around each target/feature (as specified by user command line options)
        #  later:
        #   evaluate the target scores/weights
        #   pick the feature with the best target
        
        contig_length = len(contig_sequence)
        
        # First we identify all the features from -max to +max
        feature_start = self.start
        feature_end = self.end or contig_length # if (self.end == None), then it will be len(contig_sequence)
        
        # Find the maximum distance the feature can be expanded both up- and down-stream:
        max_upstream_coord = 0
        max_downstream_coord = contig_length
        for exf_name, exf_obj in Feature.excluded_features.items():
            if (self.contig == exf_obj.contig):
                if (exf_obj.end < self.start):
                    max_upstream_coord = max(max_upstream_coord, exf_obj.end)
                if (exf_obj.start > self.end):
                    max_downstream_coord = min(max_downstream_coord, exf_obj.start)
                if (Feature.overlap_coverage(self.start, self.end, exf_obj.start, exf_obj.end) > 0):
                    logging.info("WARNING: Selected feature '{}' overlaps with excluded feature '{}'".format(self.name, exf_obj.name))
        
        logging.info('    feature: {}'.format(self.name))
        logging.info('     bounds: {}..{}'.format(self.start, self.end))
        logging.info('     limits: {}..{}'.format(max_upstream_coord, max_downstream_coord))
        
        
        
        #feature_sequence = contig_sequence[feature_start:feature_end]
        # The '10' is for upstream/downstream adjacent sequences
        #ex_feature_start = max(0, feature_start-(10+args.feature_expansion_lengths[1]))
        #ex_feature_end = min(feature_end+(10+args.feature_expansion_lengths[1]), contig_length)
        ex_feature_start = max(0, feature_start-args.feature_expansion_lengths[1])
        ex_feature_end = min(feature_end+args.feature_expansion_lengths[1], contig_length)
        
        #max_feature_sequence = contig_sequence[ex_feature_start:ex_feature_end]
        
        # We scan for targets only once
        targetsl = targets.Target.get_targets(args, contig_sequence, start=ex_feature_start, end=ex_feature_end) # Does both orientations (+/-)
        logging.info('max targets: {}'.format(len(targetsl)))
        
        # The relative coordinates of FEATURE within EX_FEATURE:
        #relative_feature_start = feature_start - ex_feature_start
        #relative_feature_end = feature_end - ex_feature_start
        
        derived_sets = set()
        
        # We generate a derived feature for each identified target
        for t in targetsl:
            # t = (orientation, start, end, upstream, downstream, filt_seq, side, filt_spacer, filt_pam, mymotif.motif_string, tuple([tuple(x) if isinstance(x, list) else x for x in mymotif.parsed_list])))
            
            target_start = t[1]
            target_end = t[2]
            
            #  center_feature: --------HHHH[...............FEATURE.........TARGET]HHHH------------------------
            #                  -----HHHH[..................FEATURE.........TARGET...]HHHH--------------------- pad=3
            #   center_target: -----------------------HHHH[FEATURE.........TARGET................]HHHH--------
            #                  --------------------HHHH[...FEATURE.........TARGET...................]HHHH----- pad=3
            #     center_both: -----------------------HHHH[FEATURE.........TARGET]HHHH------------------------
            #                  --------------------HHHH[...FEATURE.........TARGET...]HHHH--------------------- pad=3
            # justify_feature: -----------------------HHHH[FEATURE.........TARGET]HHHH------------------------
            #                  -----------------------HHHH[FEATURE.........TARGET...]HHHH--------------------- pad=3
            #  justify_target: -----------------------HHHH[FEATURE.........TARGET]HHHH------------------------
            #                  --------------------HHHH[...FEATURE.........TARGET]HHHH------------------------ pad=3
            
            # max(end2-start1, end1-start2) # Full distance including X and Y = 19
            # max(start2-end1, start1-end2) # Distance in between between X and Y = 3
            #     1 XXXXXXXXXX
            #     2              YYYYYY
            #    
            #     1          XXXXXXXXXX
            #     2 YYYYYY
            
            if (args.feature_expansion_method == 'center_feature'):
                # ----[........FFF....TTTT]----  Terminal 3'
                # ----[TTTT....FFF........]----  Terminal 5'
                # ----------[..fffTT]----------  Overlapping 3'
                # ----------[TTfff..]----------  Overlapping 5'
                # ----------[.TfffTT]----------  Complete overlap 3' longer
                # ----------[TTfffT.]----------  Complete overlap 5' longer
                full_dist = max(target_end-feature_start, feature_end-target_start)
                derived_start = feature_end - full_dist - args.feature_expansion_pad
                derived_end = feature_start + full_dist + args.feature_expansion_pad
            
            elif (args.feature_expansion_method == 'center_target'):
                # ----[....TTTT.FFF]----
                # ----[FFF.TTTT....]----
                full_dist = max(target_end-feature_start, feature_end-target_start)
                derived_start = target_end - full_dist - args.feature_expansion_pad
                derived_end = target_start + full_dist + args.feature_expansion_pad
            
            elif (args.feature_expansion_method == 'center_both'):
                # ----[TTTT...FFF]----
                # ----[FFF...TTTT]----
                derived_start = min(feature_start, target_start) - args.feature_expansion_pad
                derived_end = max(target_end, feature_end) + args.feature_expansion_pad
            
            elif (args.feature_expansion_method == 'justify_feature'):
                # ----[FFF....TTTTxxx]----
                # ----[xxxTTTT....FFF]----
                derived_start = min(feature_start, target_start)
                derived_end = max(target_end, feature_end)
                if (derived_start == feature_start):
                    derived_end += args.feature_expansion_pad
                else:
                    derived_start -= args.feature_expansion_pad
                
            elif (args.feature_expansion_method == 'justify_target'):
                # ----[TTTT....FFFxxx]----
                # ----[xxxFFF....TTTT]----
                derived_start = min(feature_start, target_start)
                derived_end = max(target_end, feature_end)
                if (derived_start == target_start):
                    derived_end += args.feature_expansion_pad
                else:
                    derived_start -= args.feature_expansion_pad
            
            # If derived feature is too small, then we automatically pad it
            # so it reaches the minimum length
            if (derived_end - derived_start < args.feature_expansion_lengths[0]):
                size_difference = args.feature_expansion_lengths[0] - (derived_end - derived_start)
                if (args.feature_expansion_method in ['center_feature', 'center_target', 'center_both']):
                    new_pad = (size_difference+1)//2
                    derived_start -= new_pad
                    derived_end += new_pad
                    logging.info("Added {} nt of padding bases to either side of DERIVED FEATURE".format(new_pad))
                
                elif (args.feature_expansion_method == 'justify_feature'):
                    if (derived_start == feature_start):
                        derived_end += size_difference
                    else:
                        derived_start -= size_difference
                    logging.info("Added {} nt of padding bases to one side of DERIVED FEATURE".format(size_difference))
                
                elif (args.feature_expansion_method == 'justify_target'):
                    if (derived_start == target_start):
                        derived_end += size_difference
                    else:
                        derived_start += size_difference
                    logging.info("Added {} nt of padding bases to one side of DERIVED FEATURE".format(size_difference))
            
            logging.info("derived_start = {}".format(derived_start))
            logging.info("  derived_end = {}".format(derived_end))
            logging.info("derived_end - derived_start = {}".format(derived_end - derived_start))
            derived_sets.add((derived_start, derived_end))
        
        # Create a counter for the derived feature
        count = 0
        
        for derived_start, derived_end in derived_sets:
            # If derived feature length is within the desired size range, then create the derived feature
            if (args.feature_expansion_lengths[0] <= derived_end - derived_start <= args.feature_expansion_lengths[1]):
                if (max_upstream_coord <= derived_start <= derived_end <= max_downstream_coord):
                    new_name = self.name + '_derived-' + str(count)
                    new_attributes = self.attributes.copy()
                    new_attributes[args.tag] = new_name
                    new_feature = Feature(self.contig, derived_start, derived_end, self.strand, name=new_name, source=self.source, feature_type=self.feature_type, score=self.score, frame=self.frame, attributes=new_attributes, origin=Feature.DERIVED, expand_parent=self)
                    Feature.features[new_name] = new_feature
                    count += 1
                    logging.info("DERIVED FEATURE '{}' created.".format(new_name))
                else:
                    logging.info("DERIVED FEATURE would overlap with excluded feature.")
        # END expand_feature()
    
    def previous_expand_feature(self, args, contig_sequence, expansion_size=20, minimum_targets_per_feature=5):
        # Import included AddTag-specific modules
        from . import targets
        
        feature_start = self.start
        feature_end = self.end
        
        contig_length = len(contig_sequence)
        
        if (feature_end == None):
            feature_end = len(contig_sequence)
        
        # Find the maximum distance the feature can be expanded both up- and down-stream:
        max_upstream_coord = 0
        max_downstream_coord = contig_length
        for exf_name, exf_obj in Feature.excluded_features.items():
            if (self.contig == exf_obj.contig):
                if (exf_obj.end < self.start):
                    max_upstream_coord = max(max_upstream_coord, exf_obj.end)
                if (exf_obj.start > self.end):
                    max_downstream_coord = min(max_downstream_coord, exf_obj.start)
                if (Feature.overlap_coverage(self.start, self.end, exf_obj.start, exf_obj.end) > 0):
                    logging.info("WARNING: Selected feature '{}' overlaps with excluded feature '{}'".format(self.name, exf_obj.name))
        
        logging.info('feature: {}'.format(self.name))
        logging.info( 'bounds: {}..{}'.format(self.start, self.end))
        logging.info(' limits: {}..{}'.format(max_upstream_coord, max_downstream_coord))
        
        # Gradually expand feature size until the minimum number of targets is found
        feature_sequence = contig_sequence[feature_start:feature_end]
        targetsl = targets.Target.get_targets(args, feature_sequence) # Does both orientations (+/-)
        logging.info('targets: {}'.format(len(targetsl)))
        count = 0
        while ((len(targetsl) < minimum_targets_per_feature) and (((feature_start, feature_end) != (0, contig_length)) or ((feature_start, feature_end) != (max_upstream_coord, max_downstream_coord)))):
            feature_start = max(0, feature_start-expansion_size, max_upstream_coord)
            feature_end = min(contig_length, feature_end + expansion_size, max_downstream_coord)
            
            feature_sequence = contig_sequence[feature_start:feature_end]
            targetsl = targets.Target.get_targets(args, feature_sequence) # Does both orientations (+/-)
            logging.info('targets: {}'.format(len(targetsl)))
            count += 1
        
        # Save the new feature as a Feature object
        if (count > 0):
            new_name = self.name + '_derived'
            new_attributes = self.attributes.copy()
            new_attributes[args.tag] = new_name
            new_feature = Feature(self.contig, feature_start, feature_end, self.strand, name=new_name, source=self.source, feature_type=self.feature_type, score=self.score, frame=self.frame, attributes=new_attributes, origin=Feature.DERIVED, expand_parent=self)
            Feature.features[new_name] = new_feature
            return new_feature
        else:
            return self
    
#    def overlap_distance(self, start1, end1, start2, end2):
#        coverage = self.overlap_coverage(start1, end1, start2, end2)
#        if (coverage > 0):
#            return -coverage
#        else:
#            #return min(abs(start2-end1), abs(start1-end2))
#            return max(start2-end1, start1-end2)
#            # 1 ..........
#            # 2              XXXXXX
#            
#            # 1          ........
#            # 2 XXXXXX
    
#    def overlap_percent_coverage(self, start1, end1, start2, end2):
#        coverage = self.overlap_coverage(start1, end1, start2, end2)
#        return float(coverage) / max(abs(end1 - start1), abs(end2 - start2))
    
#    def percent_full_length_coverage(self, start1, end1, len1, start2, end2, len2):
#        coverage = self.overlap_coverage(start1, end1, start2, end2)
#        return float(coverage) / max(len1, len2)
    
    @staticmethod
    def overlap_coverage(start1, end1, start2, end2, index_base=0):
        coverage = 0
        
        # Expects 0-based, left-inclusive, right-exclusive indexing
        if (start2 <= start1 <= end1 <= end2):
            # 1      .............
            # 2   XXXXXXXXXXXXXXXXXX
            coverage = abs(end1 - start1)
        elif (start1 <= start2 <= end2 <= end1):
            # 1      .....................
            # 2            XXXXXXXXXX
            coverage = abs(end2 - start2)
        elif (start1 <= start2 <= end1 <= end2):
            # 1      ........
            # 2        XXXXXXXX
            coverage = abs(end1 - start2)
        elif (start2 <= start1 <= end2 <= end1):
            # 1      .............
            # 2  XXXXXXXXXX
            coverage = abs(end2 - start1)
        
        # if 1-based, left-inclusive, right-inclusive indexing,
        # then give +1 to coverage
        if index_base == 1:
            coverage += 1
        
        return coverage
    
    @classmethod
    def get_overlapping_features(cls, contig, start, end):
        """Returns list (sorted) of features that overlap with arguments"""
        the_features = set()
        for feature_name, f in cls.features.items():
            if (contig == f.contig):
                if (cls.overlap_coverage(start, end, f.start, f.end) > 0):
                    the_features.add(feature_name)
        
        return sorted(the_features)