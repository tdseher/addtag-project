#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/feature.py

# Import standard packages
import logging
import itertools
import os
import sys
#from collections import Counter
from itertools import permutations, cycle
from copy import copy, deepcopy
from collections import namedtuple

# Import non-standard packages
import regex

# Import included AddTag-specific modules
#from .__init__ import Target
#from . import targets
from . import nucleotides
from . import utils

logger = logging.getLogger(__name__)

class Feature(object):
    features = {}
    excluded_features = {}
    # Origin flags
    NONE=0
    INPUT=1
    DERIVED=2
    
    logger = logger.getChild('Feature')
    
    def __init__(self, contig, start, end, strand, name=None, attributes=None, source=None, feature_type=None, score=None, frame=None, origin=NONE, sep=';', parent=None, gene=None, homologs=None):
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
            alist = regex.split(sep+r'\s*', attributes) # ['ID=12', 'Parent=nope']
            self.attributes = dict([regex.split(r'\s*=\s*', x) for x in alist]) # {'ID': '12', 'Parent': 'nope'}
        self.name = name
        self.origin = origin
        self.parent = parent
        
        if homologs:
            self.homologs = homologs
        else:
            self.homologs = []
        self.gene = gene
    
    def get_homologs(self):
        return self.homologs
    
    def get_gene(self):
        return self.gene
    
    def get_parent(self):
        if ((self.parent == None) or (self.parent == self)):
            return self
        else:
            return self.parent.get_parent()
    
    @classmethod
    def get_gene_from_feature(cls, feature_name, feature2gene):
        parent = cls.features[feature_name].get_parent().name
        
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
    def group_features_by_gene(cls, feature2gene):
        '''
        Returns a dict of Feature groups, key=gene name, value = list of Features
        :param feature2gene: 
        :return: 
        '''
        groups = {}
        for fname, f in cls.features.items():
            pname = f.get_parent().name
            gene = feature2gene[pname]
            groups.setdefault(gene, []).append(f)
        return groups
    
    @classmethod
    def assign_homologs(cls, homologs):
        """
        Give all features the 'self.homologs' attribute
        """
        # 'self.homologs' holds a list of sets (why?), and an example would look like this:
        #    [
        #      {
        #        Feature(name=C3_04520C_A, gene=ADE2, location=Ca22chr3A_C_albicans_SC5314:-:943415..945122),
        #        Feature(name=C3_04520C_B, gene=ADE2, location=Ca22chr3B_C_albicans_SC5314:-:943393..945100)
        #      },
        #      ...
        #    ]
        for feature_name, f in cls.features.items():
            if homologs:
                hname_set = homologs[f.get_parent().name]
                h_set = set() # This probably doesn't need to be a set!
                for fname, f2 in cls.features.items():
                    if fname in hname_set:
                        h_set.add(f2)
                if h_set not in f.homologs:
                    f.homologs.append(h_set)
    
    @classmethod
    def assign_gene(cls, feature2gene):
        """
        Give all features the 'self.gene' attribute
        """
        for feature_name, f in cls.features.items():
            if feature2gene:
                f.gene = feature2gene[f.get_parent().name]
    
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
    
    # TODO: Delete previous versions of 'Feature.new_new_expand_all_features()'
    
    @classmethod
    def new_new_expand_all_features(cls, args, contigs, h_groups):
        '''
        New method makes Features expand by homology group (instead of in isolation).
        :param args: Argparse namespace object
        :param contigs: Dict with key=header, value=DNA sequence
        :param h_groups: Output of 'Feature.group_features_by_gene(feature2gene)' (dict with key=gene, value=[Feature, Feature, ...]
        :return: None
        '''
        # This function does not determine if the Features should be expanded--it ASSUMES that they should be expanded.
        cls.logger.info('Starting newer Feauture expansion method')
        
        for gene, feature_list in h_groups.items():
            cls.logger.info('Working on gene: {}'.format(gene))
            
            # Find equivalent Targets, respecting 'args.target_distance'
            cls.logger.info('Finding equivalent Targets...')
            equivalents = cls.expand_for_targets(args, contigs, feature_list)
            
            # For each equivalent, we create expanded Features respecting 'args.feature_expansion_format'
            cls.logger.info('Formatting Feature+Target equivalents...')
            bounds = cls.expand_for_format(args, contigs, feature_list, equivalents)
            
            # For each formatted equivalent, we expand to find flanking homology regions with 
            # desired polymorphism level, respecting 'args.homology_distance'
            cls.logger.info('Finding flanking homology arms...')
            bounds = cls.expand_for_homology(args, contigs, feature_list, equivalents, bounds)
            
            # TODO: Need to add code that checks that the homology arm sequences don't also exist within the
            #       Feature that will be deleted, as this might cause internal recombination, and thus,
            #       the Feature will not be completely removed.
            
            # We filter by size (needs to be within '--feature_expansion_lengths MIN MAX')
            equivalents, bounds = cls.expand_size_filter(args, equivalents, bounds)
            
            # Finally, we create 'derived' Feature objects for each
            cls.logger.info("Creating derived 'Feature' objects...")
            cls.create_derived_features(args, contigs, feature_list, equivalents, bounds)
    
    @classmethod
    def new_expand_all_features(cls, args, contigs, h_groups):
        cls.logger.info('Starting new Feauture expansion method')
        
        from . import targets
        
        # for each gene/homology group:
        for gene, feature_list in h_groups.items():
            cls.logger.info('Working on gene: {}'.format(gene))
            # Find all Targets within all homologous Features
            # Create a matrix that stores which Targets are observed within which Features
            #       F1  F2  F3 ... FN
            #   T1  X   X   -      -
            #   T2  X   X   X      X
            targets_dict = {} # key = target sequence, value = list of Features containing this target sequence
            for f in feature_list:
                
                targets_list = targets.Target.get_targets(args, contigs[f.contig], start=f.start, end=f.end) # Does both orientations (+/-)
                
                for t in targets_list:
                    # t = (orientation, start, end, upstream, downstream, filt_seq, side, filt_spacer, filt_pam, mymotif.motif_string, tuple([tuple(x) if isinstance(x, list) else x for x in mymotif.parsed_list])))
                    #target_start = t[1]
                    #target_end = t[2]
                    t_seq = t[5]
                    # TODO: Somehow incorporate 'args.target_distance' to allow for fuzzy matches, thus giving more precise final output
                    x = targets_dict.setdefault(t_seq, [])
                    if f not in x:
                        x.append(f)
            
            # Temp PRINT to log
            cls.logger.info('Printing Feature/Targets matrix...')
            for k, v in targets_dict.items():
                cls.logger.info('{} {}'.format(k, v))
            
            # When matching derived Features (exFeatures) with one another, there are 3 points of homology.
            #  1) the Feature as defined in GFF input. This is assumed to be 100% identical, and will not be aligned
            #  2) the Target. If a Target (within specified 'args.target_distance') exists on all homologous Features,
            #     then it will be considered homologous
            #  3) the flanking homology regions outside the full exFeature.
            
            cls.logger.info('--target_specificity={}'.format(args.target_specificity))
            if (args.target_specificity == 'all'): # multi-allelic
                shared_targets_list = []
                # Do some Targets align against all Features?
                for t_seq, tf_list in targets_dict.items():
                    if (len(tf_list) == len(feature_list)):
                        shared_targets_list.append(t_seq)
                
                # If not, then expand all Features the same amount so they have the same Target (using args.feature_expansion_format)
                if (len(shared_targets_list) == 0):
                    pass
            
            elif (args.target_specificity == 'exclusive'): # allele-specific
                unshared_targets_dict = {} # key=feature_name, value=list of target_seqs
                # Does each Feature have a unique Target?
                for t_seq, tf_list in targets_dict.items():
                    if (len(tf_list) == 1):
                        unshared_targets_dict.setdefault(tf_list[0], []).append(t_seq)
                
                # If not, then expand all Features the same amount until each has a unique Target (using args.feature_expansion_format)
                if (len(unshared_targets_dict) != len(feature_list)):
                    pass
            
            elif (args.target_specificity == 'any'): # allele-agnostic
                feature_target_counts = {x.name: 0 for x in feature_list} # key=feature_name, value=count of Targets for that Feature
                # Does each Features have at least one Target?
                for t_seq, tf_list in targets_dict.items():
                    for tf_list_feature in tf_list:
                        feature_target_counts[tf_list_feature] += 1
                
                # If not, then expand all Features the same amount until all Features have a Target (using args.feature_expansion_format)
                if (len([k for k, v in feature_target_counts.items() if v > 0]) != len(feature_list)):
                    pass
            
            
            cls.logger.info('--donor_specificity={}'.format(args.donor_specificity))
            if (args.donor_specificity == 'all'):
                # All homology arms must be within specified distance (args.homology_distance)
                # Do quick check -- if homology arms are within edit distance (error threshold), then no further expansion needed
                # if quick check fails, then do MSA, and identify the closest US/DS homology arms that are within specified distance for all Features
                # Set each start/end for all Features in this homology group to the newly-found expanded 
                pass
            elif (args.donor_specificity == 'exclusive'):
                # All homology arms must be different
                # Do quick check
                # If quick check fails, then do MSA, and find closest US/DS regions with max polymorphism
                # set new start/end locations for all Features in homology group
                pass
            elif (args.donor_specificity == 'any'):
                # Do nothing
                pass
        pass
    
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
                cls.logger.info('Expanding feature: {}'.format(feature_name))
                f.expand_feature(args, contig_sequence)
    
    def calc_homology_similarity(self, args, contigs, f1, f2):
        """
        Calculate the similarity of two features. Takes into account:
         * Feature length (not yet)
         * Similarity of US homology region
         * Similarity of DS homology region
        
        Return a number between 0 and 1, with 1 being similar, and 0 being dissimilar
        """
        
        # Eventually, this should be calculated from '--excise_upstream_homology MIN MAX' and '--excise_downstream_homology MIN MAX'
        # The homology length is a function of '--excise_donor_lengths 100 100' and the insert size
        # We just divide the max length by two (most-extreme case is 'mintag' with 0 nt insert)
        max_homology_length = (max(args.excise_donor_lengths)+1)//2
        
        # We get the up/downstream homology regions up to this 'max_homology_length' for 'f1'
        f = f1
        us_start, us_end = f.start-max_homology_length, f.start
        f1_upstream = contigs[f.contig][us_start:us_end]
        
        ds_start, ds_end = f.end, f.end+max_homology_length
        f1_downstream = contigs[f.contig][ds_start:ds_end]
        
        # Get the homology regions for 'f2'
        f = f2
        us_start, us_end = f.start-max_homology_length, f.start
        f2_upstream = contigs[f.contig][us_start:us_end]
        
        ds_start, ds_end = f.end, f.end+max_homology_length
        f2_downstream = contigs[f.contig][ds_start:ds_end]
        
        
        # Calculate the similarity
        # Will treat any ambiguous characters as errors
        us_errors = sum(nucleotides.count_errors(f1_upstream, f2_upstream))
        ds_errors = sum(nucleotides.count_errors(f1_downstream, f2_downstream))
        
        us_len = max(len(f1_upstream), len(f2_upstream))
        ds_len = max(len(f1_downstream), len(f2_downstream))
        
        sim = (us_len+ds_len-us_errors-ds_errors)/(us_len+ds_len)
        
        return sim
    
    @classmethod
    def calc_homology_errors(cls, args, contigs, f1, f2):
        """
        Calculate the number of errors in alignment of homology regions of the two features.
        
        Returns an int.
        """
        # TODO: Modify this function so it accepts a list of features instead of only 2
        #       (See the function that calls this one 'match_features_by_homology()' for implementation.)
        
        # Eventually, this should be calculated from '--excise_upstream_homology MIN MAX' and '--excise_downstream_homology MIN MAX'
        # The homology length is a function of '--excise_donor_lengths 100 100' and the insert size
        # We just divide the max length by two (most-extreme case is 'mintag' with 0 nt insert)
        max_homology_length = (max(args.excise_donor_lengths)+1)//2
        
        # We get the up/downstream homology regions up to this 'max_homology_length' for 'f1'
        f = f1
        us_start, us_end = f.start-max_homology_length, f.start
        f1_upstream = contigs[f.contig][us_start:us_end]
        
        ds_start, ds_end = f.end, f.end+max_homology_length
        f1_downstream = contigs[f.contig][ds_start:ds_end]
        
        # Get the homology regions for 'f2'
        f = f2
        us_start, us_end = f.start-max_homology_length, f.start
        f2_upstream = contigs[f.contig][us_start:us_end]
        
        ds_start, ds_end = f.end, f.end+max_homology_length
        f2_downstream = contigs[f.contig][ds_start:ds_end]
        
        
        # Calculate the similarity
        # Will treat any ambiguous characters as errors
        us_errors = sum(nucleotides.count_errors(f1_upstream, f2_upstream))
        ds_errors = sum(nucleotides.count_errors(f1_downstream, f2_downstream))
        
        return us_errors, ds_errors
    
    @classmethod
    def match_features_by_homology(cls, args, contigs):
        """
        Group all features (including derived ones) into homologs,
        taking feature length and us/ds homology into account.
        
        Use this function after creating derived features. 
        This populates the 'feature.homologs' attribute and returns a final list of 'features_to_keep'
        """
        features_to_keep = set()
        if (args.donor_specificity == 'all'): # Multi-allelic
            
            # Create a dict where key=gene, value=set of all feature parents
            parents_per_gene = {}
            for fname, feature in cls.features.items():
                gene = feature.get_gene()
                parent = feature.get_parent()
                parents_per_gene.setdefault(gene, set()).add(parent)
            
            cls.logger.info('parents_per_gene:')
            for k, v in parents_per_gene.items():
                cls.logger.info(' {} {}'.format(k, v)) # BRG1 {'C1_05140W_B', 'C1_05140W_C', 'C1_05140W_A'}
            
            # We will need to do this for each gene
            for G, parents in parents_per_gene.items():
                # Create a dict where key=parent, value=feature
                vals = {}
                for fname, feature in cls.features.items():
                    gene = feature.get_gene()
                    parent = feature.get_parent()
                    if (gene == G):
                        vals.setdefault(parent, list()).append(feature)
                order = sorted(vals, key=lambda x: x.name)
                cls.logger.info('order = {}'.format(order)) # order = ['C1_05140W_A', 'C1_05140W_B', 'C1_05140W_C']
                
                # Separate features into lists based on their parent
                # All derived features from the same parent will be in a list together
                odat = [vals[x] for x in order]
                cls.logger.info('odat = {}'.format(odat)) # odat = [[Feature(C1_05140W_A_derived-0), Feature(C1_05140W_A_derived-1), Feature(C1_05140W_A_derived-2)], [Feature(C1_05140W_B_derived-0), Feature(C1_05140W_B_derived-1)], [Feature(C1_05140W_C_derived-0)]]
                
                # Perform pairwise comparisons for every combination of input homologies
                similarity = {}
                c_similarities = {}
                couplings = list(itertools.product(*odat))
                for c1 in couplings:
                    cls.logger.info(c1)
                    sims = []
                    comparisons = list(itertools.combinations(c1, 2))
                    for c2 in comparisons:
                        cls.logger.info('  {}'.format(c2))
                        #sims.append(similarity.setdefault(c2, self.calc_homology_similarity(args, contigs, c2[0], c2[1])))
                        sim = similarity.get(c2, None)
                        if sim:
                            sims.append(sim)
                        else:
                            sims.append(similarity.setdefault(c2, cls.calc_homology_errors(args, contigs, c2[0], c2[1])))
                        
                    c_similarities[c1] = sims
                    # (Feature(C1_05140W_A_derived-0), Feature(C1_05140W_B_derived-0), Feature(C1_05140W_C_derived-0))
                    #   (Feature(C1_05140W_A_derived-0), Feature(C1_05140W_B_derived-0))
                    #   (Feature(C1_05140W_A_derived-0), Feature(C1_05140W_C_derived-0))
                    #   (Feature(C1_05140W_B_derived-0), Feature(C1_05140W_C_derived-0))
                    
                
                # Test every possible input homology group to see if they are similar enough
                cls.logger.info('c_similarities:')
                for k, v in c_similarities.items():
                    #cls.logger.info(' ', k, v, all(x > 0.95 for x in v)) # (Feature(C1_05140W_A_derived-0), Feature(C1_05140W_B_derived-0), Feature(C1_05140W_C_derived-0)) [0.99, 0.98, 0.97] True
                    verdict = all(((x[0] <= args.homology_distance[2]) and (x[1] <= args.homology_distance[2])) for x in v)
                    cls.logger.info(' {} {} {}'.format(k, v, verdict))
                    if verdict:
                        for f in k:
                            features_to_keep.add(f)
                            # Add this new set of homologs to the 'self.homologs' attribute
                            # (each feature can thus have multiple sets of homologs)
                            # Ideally, it should resemble the following:
                            #   self.homologs = [[Feature(), Feature(), Feature()], [Feature(), Feature(), Feature()]]
                            hset = set(k)
                            if hset not in f.homologs:
                                f.homologs.append(hset)
                            
        
        elif (args.donor_specificity == 'exclusive'): # Uni-alleleic
            # TODO: Finish writing Feature match function for exclusive/uni-allelic case
            # Either the length of the insert should be diagnostically different
            # Or the homology regions should have polymorphisms (maximize polymorphisms)
            # or both
            
            
            # Placeholder code
            for fname, f in cls.features.items():
                features_to_keep.add(f)
            
        elif (args.donor_specificity == 'any'): # Allele-agnostic
            # TODO: Finish writing Feature match function for any/allele-agnostic case
            # Placeholder code
            for fname, f in cls.features.items():
                features_to_keep.add(f)
    
        return list(features_to_keep)
    
    @classmethod
    def filter_features(cls, feature_names):
        """
        Reduce the total number of features to just the ones indicated in the selection
        """
        
        # Require at least one in order to filter
        # Feature has this format: key=gene/tag, value = (contig, start(bp), end(bp), strand)
        if feature_names:
            new_features = {}
            for fn in feature_names:
                f = cls.features.get(fn)
                if f:
                    new_features[fn] = f
                else:
                    raise Exception("Feature specified as '"+fn+"' does not exist in 'Feature.features'.")
                    #sys.exit(1)
            if (len(new_features) == 0):
                raise Exception("No features to be filtered are present in 'Feature.features'.")
                #sys.exit(1)
            cls.features = new_features
    
    def expand_feature_for_homology(self, args, contigs):
        """
        Expand the feature to ensure its upstream/downstream flanking sequences
        have no polymorphisms
        
        Assumes each 'Feature.homologs' list is populated
        """
        
        from . import aligners
        
        for side in ('upstream', 'downstream'):
            dist = 0
            found = False
            while (found == False):
                dist += 100
                
                # Get upstream and downstream regions of each homolog
                hom_list = []
                
                for f in self.homologs:
                    hom_list.append(None)
                    
                folder = os.path.join(args.temp_folder, 'addtag', os.path.basename(args.folder))
                query_filename = os.path.join(folder, side+'.fasta')
                
                # Store upstream and downstream regions into a multi-fasta file
                with open(query_filename, 'w') as flo:
                    for i, seq in enumerate(hom_list):
                        print('>{}\n{}'.format(i, seq), file=flo)
                
                
                # Get the MSA aligner object
                for a in aligners.ms_aligners:
                    if (a.name == 'mafft'):
                        aligner = a
                        break
                
                # Align sequences using MSA
                aln_filename = aligner.align(us_query_filename, side+'.aln', folder, args.processors)
                
                # Read the MSA
                # Extract the longest region with perfect homology
                # check to see if region is long enough
                if (hom_region_len >= min_hom_len):
                    found = True
        
    
    def expand_homologous_feature_bad(self, args, contigs, homologs):
        """
        Expand the feature to ensure its upstream/downstream flanking sequences
        have no polymorphisms
        """
        
        # Should run MSA on upstream/downstream regions for all homologous features
        # From this, pick the shared us/ds homology regions
        # and choose different sizes to expand each feature depending on this MSA
        
        
        #min_homology_length = 50 # Eventually, this should be calculated from '--excise_upstream_homology MIN MAX' and '--excise_downstream_homology MIN MAX'
        # The homology length is a function of '--excise_donor_lengths 100 100' and the insert size
        # We just divide the max length by two (most-extreme case is 'mintag' with 0 nt insert)
        min_homology_length = (max(args.excise_donor_lengths)+1)//2
        
        # Get the IDs of the features homologous to self
        try:
            homolog_set = homologs[self.name]
        except KeyError:
            # This feature name does not have a homolog in the 'homologs' dict (probably because it is a derived feature)
            homolog_set = set()
        
        # Get the objects associated with each of the homologous feature names
        homolog_features = []
        
        for fname, f in Feature.features.items():
            if fname in homolog_set: # <== temporary code: needs to be fixed
                if (self != f):
                    if (f.origin == Feature.INPUT):
                        homolog_features.append(f)
        
        lcs_size = 0
        lcs_start, lcs_end
        n = 20
        should_end = False
        while (not should_end and (lcs_size < min_homology_length)):
            
            if (self.start - n < 0):
                should_end = True
            else:
                ref_us = contigs[self.contig][self.start-n:self.start]
                
                for f in homolog_features:
                    if (f.start - n < 0):
                        should_end = True
                    else:
                        f_us = contigs[f.contig][f.start-n:f.start]
                        m = nucleotides.lcs(ref_us, f_us)
                lcs_size = m.size
                lcs_start = self.start-n + m.a
                lcs_end = self.start-n + m.a + m.size
            
            n += 1
        
            
        # We search the upstream region
        n = 1 # Start at 20 nt upstream of feature start
        should_end = False
        while(should_end == False):
            if (self.start-n < 0):
                should_end = True
                n = self.start # Limit n
                break
            else:
                ref_us_region = contigs[self.contig][self.start-n:self.start]
                lcs_lengths = []
                for f in homolog_features:
                    if (f.start - n < 0):
                        should_end = True
                        n = f.start # Limit n
                    else:
                        f_us_region = contigs[f.contig][f.start-n:f.start]
                        m = nucleotides.lcs(ref_us_region, f_us_region)
                        lcs_lengths.append(m.size)
                if (min(lcs_lengths) > min_homology_length):
                    # Set n to the precise distance to make min(lcs_lengths) equal min_homology_length
                    n = n-(min(lcs_lengths) - min_homology_length)
                    should_end = True
                    break
            n += 10
    
    @classmethod
    def quick_msa(cls):
        # Ideas
        
        ### Method 1 - no expansion ###
        # position_lists = []
        # For each pairwise combination
        #   match_positions = []
        #   gap_positions = []
        #   error_positions = []
        #   
        #   use regex to compare just the homology arms (not the full expansion)
        #   Find the M/G/E, and note its position
        #   position_lists.append(match_positions, gap_positions, error_positions)
        #    
        # for each position with M/G/E
        #   if G:
        #     g_count += 1
        #     e_count += 1
        #   if M:
        #     m_count += 1
        #     e_count += 1
        
        ### Method 2 - with anchored expansion ###
        # For each pairwise combination
        #   do a regex search with specified M/G/E of first homology region vs expanded second
        #   note any match locations
        # 
        # if 'all'
        #   make list of match locations shared by all pairwise combinations
        #   if upstream
        #     take right-most match location shared
        #   if downstream
        #     take left-most match location shared
        # if 'exclusive'
        #   for each feature
        #     make list of locations without matches
        #     if upstream
        #       take right-most exclusive location
        #     if downstream
        #       take left-most exclusive location
        # if 'any'
        #   don't do any matching, and just take derived_start, derived_end coordinates as-is
        # 
        # An example:
        #   a  TTTTTTT
        #   b  GTATAGG
        #   c  CGACTAT
        #   af TTTTTTTCAATCAGCTAAGCCGATACGA
        #   bf GTATAGGCCTTTTTTTAGCCTAGAGCTA
        #   cf CGACTATTTTTTTGCATAAGACAAACCA
        #   
        #   a           TTTTTTT
        #   bf GTATAGGCCTTTTTTTAGCCTAGAGCTA
        #   
        #   a        TTTTTTT
        #   cf CGACTATTTTTTTGCATAAGACAAACCA
        #   
        #   b                      GTATAGG
        #   af TTTTTTTCAATCAGCTAAGCCGATACGA
        #   
        #   b               GTATAGG
        #   cf CGACTATTTTTTTGCATAAGACAAACCA
        #   
        #   c              CGACTAT
        #   af TTTTTTTCAATCAGCTAAGCCGATACGA
        #   
        #   c                  CGACTAT
        #   bf GTATAGGCCTTTTTTTAGCCTAGAGCTA
        #   
        #   In this example, the 'a' homology arm is best, so this is the one that will be used
        pass
    
    @classmethod
    def create_derived_features(cls, args, contigs, feature_list, equivalents_list, bounds_list):
        
        cls.logger.info("Running function 'create_derived_features()'...")
        
        # Non-redundant list of features to add (same as a 'bounds' tuple)
        #bounds_set = set() # (fi, start, end)
        bounds_dict = {} # key=(fi, start, end), value=[i, i, ...] list of equivalents group indices
        targets_dict = {}
        new_features_dict = {}
        
        for i, (el, bl) in enumerate(zip(equivalents_list, bounds_list)):
            # We make non-redundant list of features to add
            #for bounds in bl:
            #    #bounds_set.add(bounds)
            #    bounds_dict.setdefault(bounds, []).append(i) # For every bounds, we make their homologs list
            
            #for (e_fi, e_t), (d_fi, d_start, d_end) in zip(el, bl):
            for (fi, t), bounds in zip(el, bl):
                # For every bounds, we make their homologs list
                bounds_dict.setdefault(bounds, []).append(i)
                
                # For every bounds, we make their Targets/spacers list
                targets_dict.setdefault(bounds, []).append(t)
        
        # For every bounds, we create the 'Feature' object, starting from smallest to largest
        #for i, (d_fi, d_start, d_end) in enumerate(sorted(bounds_set, key=lambda x: x[2]-x[1])):
        
        for i, (bounds, eqi_list) in enumerate(bounds_dict.items()):
            # TODO: Ideally, every 'Feature.name' in an equivalence group would have the same '_derived-N' suffix
            #       However, the program is currently designed so only one Feature object per sequence can exist.
            #       If there are multiple Feature objects with the same sequence, then there will likely be a problem
            #       So as a work-around, the names will look like this: '_derived-1/3/5/7'
            #       ('/' because ',' is already taken) with a comma-separated list of homology groups
            d_fi, d_start, d_end = bounds
            f = feature_list[d_fi]
            new_name = '{}_derived-{}'.format(f.name, '/'.join(map(str, eqi_list)))
            new_attributes = f.attributes.copy()
            new_attributes[args.tag] = new_name
            new_feature = Feature(f.contig, d_start, d_end, f.strand, name=new_name, source=f.source, feature_type=f.feature_type, score=f.score, frame=f.frame, attributes=new_attributes, origin=Feature.DERIVED, parent=f, gene=f.gene)
            Feature.features[new_name] = new_feature
            new_features_dict[bounds] = new_feature
            
            cls.logger.info("NEW DERIVED FEATURE created ({}): {}".format(i, new_feature))
        cls.logger.info('') # Blank line
        
        # Now that the derived Features were created, we add their homolog groups to the 'Feature.homologs' attributes
        # Feature.homologs = [set(Feature, Feature, ...), set(Feature, Feature, ...)]
        for bl in bounds_list:
            s = set(new_features_dict[bounds] for bounds in bl)
            
            for f in s:
                if s not in f.homologs:
                    f.homologs.append(s)
        
        # Log how homologs appear after adding them
        cls.logger.info('Added homologs to derived features: (d_fi, d_start, d_end), Feature.name, Feature.homologs')
        for bounds, f in new_features_dict.items():
            cls.logger.info(' {}, {}, {}'.format(bounds, f.name, f.homologs))
        
        # Get index order from smallest average Feature to largest
        #for i, (el, bl) in enumerate(sorted(zip(equivalents_list, bounds_list), key=lambda x: sum(y[2]-y[1] for y in x[1])/len(x[1]))):
        #    cls.logger.info('{}: avg_len={}'.format(i, sum(y[2]-y[1] for y in bl)/len(bl)))
        #    for e in el:
        #        cls.logger.info('  {}'.format(e))
        #    for b in bl:
        #        cls.logger.info('  {}'.format(b))
            
            # Several bound entries (d_fi, d_start, d_end) will be identical even though they all have different Target (equivalents)
            # These all should only make a single Feature, with several Targets linked /homologs added?
            # Later on, Target objects are created in the '_subroutine_generate_all.py' file by the following commands:
            #   ExcisionTarget.search_all_features(args, contig_sequences)
            #     and
            #   ReversionTarget.get_targets()
            # And in the 'donors.py' file by the following commands in the 'ExcisionDonor' class:
            #   cls(f.name, f.contig, orientation, dDNA, (start1, end1), mAT, (start2, end2), spacer=t)
            #   cls(f.name, f.contig, orientation, dDNA, (us_start, us_end), addtag, (ds_start, ds_end), spacer=t)
            #   cls(f.name, f.contig, orientation, dDNA, (us_start, us_end), unitag, (ds_start, ds_end), spacer=t)
            # No Target objects are created by the 'ReversionDonor' class
            
            # An alternative would be to create a single Target now, and add it to the 'Donor.spacers' set for each 
            # derived Feature this Equivalent should target.
            
            # One problem is that Target objects are sequence-specific, and do not respect equivalence (yet)
            # Should each Target's equivalence be stored in the 'Target' class?
            #   option 1:
            #     Target.equivalents = {'ACAGCT': ['CCAGCT', 'GCAGCT', ...]} # These are the 'key's to Target.sequences
            #   option 2:
            #     Target.equivalents = {'exTarget-12': {Target, Target, ...}} # key=Target.name, value=set(Target, Target, ...)
            # Option 2 won't work, because at this stage, Target objects don't exist, so they don't have 'Target.name'
            # So it NEEDS to be option 1.
            
            
            
            
        
        # for derived_start, derived_end in sorted(derived_sets, key=lambda x: x[1]-x[0]): # in order from smallest to largest
        #     # If derived feature length is within the desired size range, then create the derived feature
        #     if (args.feature_expansion_lengths[0] <= derived_end - derived_start <= args.feature_expansion_lengths[1]):
        #         if (max_upstream_coord <= derived_start <= derived_end <= max_downstream_coord):
        #             new_name = self.name + '_derived-' + str(count)
        #             new_attributes = self.attributes.copy()
        #             new_attributes[args.tag] = new_name
        #             new_feature = Feature(self.contig, derived_start, derived_end, self.strand, name=new_name, source=self.source, feature_type=self.feature_type, score=self.score, frame=self.frame, attributes=new_attributes, origin=Feature.DERIVED, parent=self, gene=self.gene)
        #             Feature.features[new_name] = new_feature
        #             count += 1
        #             self.logger.info("DERIVED FEATURE '{}' created.".format(new_name))
        #         else:
        #             self.logger.info("DERIVED FEATURE would overlap with excluded feature.")
        
    
    @classmethod
    def expand_size_filter(cls, args, equivalents_list, bounds_list):
        new_equivalents_list = []
        new_bounds_list = []
        
        # If any Feature in the group fails, then none of the equivalents are added as derived Features
        for el, bl in zip(equivalents_list, bounds_list):
            if all([args.feature_expansion_lengths[0] <= d_end-d_start <= args.feature_expansion_lengths[1] for d_fi, d_start, d_end in bl]):
                new_equivalents_list.append(el)
                new_bounds_list.append(bl)
        
        return new_equivalents_list, new_bounds_list
    
    @classmethod
    def expand_for_homology(cls, args, contigs, feature_list, equivalents_list, bounds_list):
        # Define convenience variables
        mismatches, gaps, errors = args.homology_distance
        
        # Import the aligners
        from . import aligners
        
        # Get the MSA aligner object
        aligner = None
        for a in aligners.ms_aligners:
            if (a.name == 'mafft'):
                aligner = a
                break
        
        # Calculate the length of the homology arm to search for
        length = (max(args.excise_donor_lengths)+1)//2
        
        # Make a 'dict' for storing already-computed homologs for this gene (feature group)
        #   key=((fi, start, end), (fi, start, end), ...) for input bounds to the MSA
        #   value=records (as parsed)
        # TODO: make program skip the records analysis as well. Thus:
        #       value=(d_start, d_end) with MSA result, or (start, end) if no MSA performed 
        msa_dict = {}
        
        # 'args.feature_expansion_lengths[1]' is the max distance up/down-stream to expand the feature
        # This is the length of the homology region to search within
        
        # Define an empty list to populate with results
        # Should have the following format (where the lists are equivalents, 
        # and the tuples are Features within the equivalent):
        #   new_bounds_list = [
        #     [(feature_index, derived_start, derived_end), (feature_index, derived_start, derived_end), ...],
        #     [(feature_index, derived_start, derived_end), (feature_index, derived_start, derived_end), ...],
        #     ...
        #   ]
        new_bounds_list = []
        
        # If/elif/elif tatement for selecting between 'all', 'exclusive', and 'any' dDNA homologies based on input parameters
        if (args.donor_specificity == 'any'): # Allele-agnostic
            # For allele-agnostic ('any') dDNA, then no MSA is performed
            
            # Thus, we just return the 'bounds_list' with just -/+ length to US/DS
            # for blist in bounds_list:
            #     new_bounds = []
            #     for fi, d_start, d_end in blist:
            #         new_bounds.append((fi, max(0, d_start-length), min(d_end+length, len(contigs[feature_list[fi].contig]))))
            #     new_bounds_list.append(new_bounds)
            #     
            # return new_bounds_list
            
            return bounds_list
        
        else: # 'all' or 'exclusive'
            for eqi, (elist, blist) in enumerate(zip(equivalents_list, bounds_list)):
                gene = feature_list[0].get_gene() # They should all be the same gene
                
                # 'good_windows' stores closest determined flanking homology arms
                # good_windows = [
                #   (us_consensus_seq, [WindowRecord(), WindowRecord(), ...]),
                #   (ds_consensus_seq, [WindowRecord(), WindowRecord(), ...]),
                # ]
                good_windows = [None, None] # [US, DS]
                
                for side_i, side in enumerate(['us', 'ds']):
                    logging.info('Working on MSA for side: {}'.format(side))
                    msa_sequences = []
                    msa_names = []
                    msa_contigs = []
                    msa_starts = []
                    msa_ends = []
                    msa_feature_indices = []
                    
                    if side in ['us', 'US', 'left', "5'", 'upstream']:
                        # Find the US region of all homologous features
                    
                        for (fi, t), (fi2, d_start, d_end) in zip(elist, blist):
                            f = feature_list[fi]
                            
                            us_start = max(0, d_start-args.feature_expansion_lengths[1]-length)
                            us_end = min(d_start, len(contigs[f.contig]))
                            
                            msa_sequences.append(contigs[f.contig][us_start:us_end])
                            msa_names.append('{}:EQ{}:{}:{}:{}..{}'.format(gene, eqi, f.name, f.contig, us_start, us_end))
                            msa_contigs.append(f.contig)
                            msa_starts.append(us_start)
                            msa_ends.append(us_end)
                            msa_feature_indices.append(fi)
                    
                    elif side in ['ds', 'DS', 'right', "3'", 'downstream']:
                        # Find the DS region of all homologous features
                        
                        for (fi, t), (fi2, d_start, d_end) in zip(elist, blist):
                            f = feature_list[fi]
                            
                            us_start = max(0, min(d_end, len(contigs[f.contig])))
                            us_end = min(d_end+args.feature_expansion_lengths[1]+length, len(contigs[f.contig]))
                            
                            msa_sequences.append(contigs[f.contig][us_start:us_end])
                            msa_names.append('{}:EQ{}:{}:{}:{}..{}'.format(gene, eqi, f.name, f.contig, us_start, us_end))
                            msa_contigs.append(f.contig)
                            msa_starts.append(us_start)
                            msa_ends.append(us_end)
                            msa_feature_indices.append(fi)
                    
                    # TODO: Make a quick-align function to try first, BEFORE doing a third-party MSA
                    #       see 'quick-msa()' for my ideas.
                    #       If the homologous flanking regions pass 'quick-msa()', then there is no
                    #       reason to do a potentially expensive/erroneous third-party MSA
                    #       Perhaps, the 'quick-msa()' will only naively calculate whether-or-not expansion
                    #       is necessary, and not attempt to do it at all.
                    
                    # Check to see if the MSA was already performed.
                    # If so, then we skip the redundant calculation, and just retrieve the previous result.
                    msa_key = tuple(zip(msa_feature_indices, msa_starts, msa_ends))
                    records = msa_dict.get(msa_key) # Defaults to None if key does not exist in dict
                    
                    if not records:
                        # MSA FASTA files need to be in their own sub-folder because there are so many
                        msa_folder = os.path.join(args.folder, 'msa')
                        os.makedirs(msa_folder, exist_ok=True)
                        
                        # Write sequences to FASTA file
                        msa_input_path = utils.write_merged_fasta((msa_names, msa_sequences), os.path.join(msa_folder, 'msa-{}-EQ{}-{}-input.fasta'.format(utils.slugify(gene), eqi, side)))
                        
                        # Perform MSA
                        msa_output_path = aligner.align(msa_input_path, 'msa-{}-EQ{}-{}-output.fasta'.format(utils.slugify(gene), eqi, side), msa_folder, args.processors)
                        
                        # Parse MSA output to get records
                        # We assume the ordering of sequences in 'records' mirrors the input sequence ordering
                        records = list(aligner.load(msa_output_path))
                        
                        # Store records in 'msa_dict'
                        msa_dict[msa_key] = records
                    
                    # Print records to log
                    cls.logger.info('records = [')
                    for r in records:
                        cls.logger.info('  {}'.format(r))
                    cls.logger.info(']')
                    
                    # Calculate conservation of MSA
                    conservation_list, consensus_list, counts_list, freqs_list = cls.calculate_msa_conservation(records, ambiguities=False)
                    
                    # Log the MSA
                    for line in cls.format_msa(records, conservation_list, consensus_list, counts_list):
                        cls.logger.info(line)
                    
                    # TODO: Put the 'window_iterator' code (below) into its own separate function
                    
                    if side in ['us', 'US', 'left', "5'", 'upstream']:
                        # We find right-most acceptable homology region
                        #window_iterator = range(len(records[0].sequence)-length, -1, -1)
                        window_iterator = range(len(records[0].sequence)-length+1)[::-1]
                    elif side in ['ds', 'DS', 'right', "3'", 'downstream']:
                        # We find left-most acceptable homology region
                        window_iterator = range(len(records[0].sequence)-length+1)
                    
                    # Identify all regions of acceptable 'length' that fulfill stringency requirements
                    # An alternative method:
                    #   Extract the longest region with acceptable homology from the MSA, 
                    #   then check to see if the region is long enough (>length)
                    
                    for w in window_iterator: # w is the position in the alignment of the candidate window
                        
                        all_records_passed = False
                        
                        mismatch_counts = [0] * len(records)
                        gap_counts = [0] * len(records)
                        error_counts = [0] * len(records)
                        
                        if (args.donor_specificity == 'all'): # Multi-allelic
                            ###### Start multi-allelic ('all') ######
                            mismatch_count = 0
                            gap_count = 0
                            error_count = 0
                            
                            for i in range(w, w+length):
                                nt_list = [r.sequence[i] for r in records]
                                
                                # Handle mismatches (treats any level of conservation as equivalent)
                                # TODO: Make it so conservation level within mismatch positions actually matters
                                if (len(nt_list)-counts_list[i]-nt_list.count('-') > 0):
                                    error_count += 1
                                    mismatch_count += 1
                                
                                # Handle gaps
                                if ('-' in nt_list):
                                    error_count += 1
                                    gap_count += 1
                                
                                if (mismatch_count <= mismatches):
                                    mismatch_check_passed = True
                                else:
                                    mismatch_check_passed = False
                                
                                if (gap_count <= gaps):
                                    gap_check_passed = True
                                else:
                                    gap_check_passed = False
                                
                                if (error_count <= errors):
                                    error_check_passed = True
                                else:
                                    error_check_passed = False
                                
                                if all([mismatch_check_passed, gap_check_passed, error_check_passed]):
                                    all_records_passed = True
                                else:
                                    all_records_passed = False
                                    break
                            
                            mismatch_counts = [mismatch_count for m in mismatch_counts]
                            gap_counts = [gap_count for m in gap_counts]
                            error_counts = [error_count for m in error_counts]
                            
                            ###### End multi-allelic ('all') ######
                        
                        elif (args.donor_specificity == 'exclusive'): # Uni-alleleic/allele-specific
                            ###### Start allele-specific ('exclusive') ######
                            
                            # For allele-specific ('exclusive') dDNA, each Feature within the homology group(gene) should have more than the number of permitted mismatches/gaps/errors
                            #
                            # If a residue is unique (not shared) amongst any other record at a position, it will be
                            # recorded as '1' or '2'. If the residue is redundant, then it will be a '0'
                            # Will turn this:
                            #   [
                            #     'GTTGGCTCGAATYGCGCGCTTA-----TTACGCTAAAATAGATCGCGCTAAGGCTATAATAGG-AGCTCCCGGAGCGGCCCGCGAA',
                            #     '----GTTGGCTTCGCGCGCAAA-----TTACGCTAAAATAACTCGCGCTAAGGCTATAATAGGAAGCTCCCGGAGCGGCC-GCGAA',
                            #     'ATACACTCGAATCGCAAACTTAGGATTTTACGCTAAAATAGATCGGGCTAAGGCTATACTAGG-AGCTC---GAGCGGCCCGCGAA',
                            #     'ATACACTGGAATC-CAAACTTAGGATT---CGCTACAATAGATCGGGCTAAGGCTATACTAGGCAGCTC---GAGCGGCCAGCGAA',
                            #   ]
                            # Into this:
                            #   [
                            #     '10110000000010000000000000000000000000000000000000000000000000000000000000000000000000',
                            #     '22220100011000000001100000000000000000001100000000000000000000010000000000000000200000',
                            #     '00000000000000000000000000000000000000000000000000000000000000000000000000000000000000',
                            #     '00000000000002000000000000022200000100000000000000000000000000010000000000000000100000',
                            #   ]
                            IDENTICAL = 0b0
                            MISMATCH = 0b1
                            GAP = 0b10
                            uniques = []
                            #for i in range(len(records[0].sequence)):
                            for i in range(w, w+length):
                                ulist = []
                                nt_list = [r.sequence[i] for r in records]
                                for nti1, nt1 in enumerate(nt_list):
                                    u = IDENTICAL
                                    for nti2, nt2 in enumerate(nt_list):
                                        if (nti1 != nti2):
                                            if (nt1 == nt2):
                                                u = IDENTICAL # If the residue is duplicated, then it is not unique
                                                break
                                            if (nt1 == '-'):
                                                u |= GAP
                                            else:
                                                u |= MISMATCH
                                    ulist.append(u)
                                uniques.append(ulist)
                            unique_list = [[u[i] for u in uniques] for i in range(len(records))]
                            useq_list = [''.join(str(y) for y in x) for x in unique_list]
                            
                            # Print useq strings to log
                            cls.logger.info('uniques = [')
                            for useq in useq_list:
                                cls.logger.info('  {}'.format(useq))
                            cls.logger.info(']')
                            
                            # We count up the non-zero entries. For each record's uniqueness. If it is > the specified
                            # M/G/E, then we call that window "unique" for that record only.
                            unique_windows = [] # One element per record. True if unique, otherwise False
                            for ui, ulist in enumerate(unique_list):
                                #u_mismatch_count = 0
                                #u_gap_count = 0
                                #u_error_count = 0
                                for u in ulist:
                                    if (u & MISMATCH):
                                        #u_mismatch_count += 1
                                        #u_error_count += 1
                                        mismatch_counts[ui] += 1
                                        error_counts[ui] += 1
                                    if (u & GAP):
                                        #u_gap_count += 1
                                        #u_error_count += 1
                                        gap_counts[ui] += 1
                                        error_counts[ui] += 1
                                
                                #if ((u_mismatch_count > mismatches) or (u_gap_count > gaps) or (u_error_count > errors)):
                                if ((mismatch_counts[ui] > mismatches) or (gap_counts[ui] > gaps) or (error_counts[ui] > errors)):
                                    unique_windows.append(True)
                                else:
                                    unique_windows.append(False)
                            
                            # If all records have greater than the threshold of permissable errors with each other,
                            # Then the window is allele-specific for all of them
                            if all(unique_windows):
                                all_records_passed = True
                            
                            ###### End allele-specific ('exclusive') ######
                        
                        # If 'all_records_passed', then this window of the MSA should be a candidate for a flanking homology region
                        if all_records_passed:
                            # TODO: The window at (w, w+length) will be the MAXIMUM homology length. Thus, if there is
                            #       a position with the alignment ['A', '-', '-', '-'], then it will be 'A' in the window,
                            #       where the majority/consensus would list it as '-' (to be removed)
                            #       Thus, the window must be filtered to remove majority '-' characters
                            #       BEFORE the window of size 'length' is searched
                            
                            WindowRecord = namedtuple('WindowRecord', ['msa_start', 'msa_end', 'msa_sequence', 'contig', 'start', 'end', 'sequence', 'feature_index', 'feature', 'mismatches', 'gaps', 'errors'])
                            
                            window_list = []
                            
                            # We need to convert the window's MSA coordinates into genomic coordinates
                            for ri, r in enumerate(records):
                                window_start = msa_starts[ri] + w - r.sequence[:w].count('-')
                                window_end = msa_starts[ri] + w + length - r.sequence[:w+length].count('-')
                                wr = WindowRecord(
                                    msa_start=w,
                                    msa_end=w+length,
                                    msa_sequence=r.sequence[w:w+length],
                                    #msa_consensus=''.join(consensus_list[w:w+length]),
                                    contig=msa_contigs[ri],
                                    start=window_start,
                                    end=window_end,
                                    sequence=contigs[msa_contigs[ri]][window_start:window_end],
                                    feature_index=msa_feature_indices[ri],
                                    feature=feature_list[msa_feature_indices[ri]].name,
                                    mismatches=mismatch_counts[ri],
                                    gaps=gap_counts[ri],
                                    errors=error_counts[ri],
                                )
                                window_list.append(wr)
                            
                            # Print to log
                            cls.logger.info('window_list = [')
                            for wr in window_list:
                                cls.logger.info('  {}'.format(wr))
                            cls.logger.info(']')
                            cls.logger.info('') # Blank line
                            
                            # TODO: What if the alignment has gaps, so the actual window is less than 'length' (50 nt)?
                            #       nevermind... the 'consensus' is used!
                            
                            # We add all this info to our master list of good windows
                            good_windows[side_i] = (''.join(consensus_list[w:w+length]), window_list)
                            
                            # We terminate after finding the first good window
                            break
                
                # Show abbreviated contents of 'good_windows'
                cls.logger.info('good_windows = [')
                for gw in good_windows:
                    if gw:
                        cls.logger.info('  ({}, ...)'.format(gw[0]))
                    else:
                        cls.logger.info('  None')
                cls.logger.info(']')
                cls.logger.info('') # Blank line
                
                # Print info message if no good windows found for either US or DS region
                if None in good_windows:
                    cls.logger.info('No adequate Windows meeting divergence thresholds found!')
                    cls.logger.info("Falling back on '--donor_specificity=any' implementation for this equivalent set")
                
                # We convert the 'good_windows' into new bounds (start/end locations) for the 'new_bounds_list'
                d_fis = []
                d_starts = []
                d_ends = []
                
                if (good_windows[0] == None):
                    d_fis.append(bounds_list[eqi][0])
                    #d_starts.append(max(0, bounds_list[eqi][1]-length)) # When including homology arm
                    d_starts.append(bounds_list[eqi][1])
                else:
                    for wr in good_windows[0][1]: # US
                        d_fis.append(wr.feature_index)
                        #d_starts.append(wr.start) # When including homology arm
                        d_starts.append(wr.end)
                
                if (good_windows[1] == None):
                    #d_ends.append(min(bounds_list[eqi][2]+length, len(contigs[feature_list[bounds_list[eqi][0]].contig]))) # When including homology arm
                    d_ends.append(bounds_list[eqi][2])
                else:
                    for wr in good_windows[1][1]: # DS
                        #d_ends.append(wr.end) # When including homology arm
                        d_ends.append(wr.start)
                
                new_bounds = list(zip(d_fis, d_starts, d_ends))
                new_bounds_list.append(new_bounds)
        
        # Log old and new bounds to make sure they were expanded correctly
        cls.logger.info('OLD BOUNDS = [')
        for i, b in enumerate(bounds_list):
            cls.logger.info('  {}: {}'.format(i, b))
        cls.logger.info(']')
        cls.logger.info('NEW BOUNDS = [')
        for i, b in enumerate(new_bounds_list):
            cls.logger.info('  {}: {}'.format(i, b))
        cls.logger.info(']')
        cls.logger.info('')
        
        return new_bounds_list
    
    @classmethod
    def expand_for_format(cls, args, contigs, feature_list, equivalents_list):
        
        feature_bounds = [] # Same internal structure as 'equivalents_list'
        
        for eqi, be in enumerate(equivalents_list):
            fb = []
            for fi, t in be:
                # t = (orientation, start, end, upstream, downstream, filt_seq, side, filt_spacer, filt_pam, mymotif.motif_string, tuple([tuple(x) if isinstance(x, list) else x for x in mymotif.parsed_list])))
                
                f = feature_list[fi]
                
                feature_start = f.start
                feature_end = f.end or len(contigs[f.contig]) # if (f.end == None), then it will be len(contig_sequence)
                
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
                
                if (args.feature_expansion_format == 'center_feature'):
                    # ----[........FFF....TTTT]----  Terminal 3'
                    # ----[TTTT....FFF........]----  Terminal 5'
                    # ----------[..fffTT]----------  Overlapping 3'
                    # ----------[TTfff..]----------  Overlapping 5'
                    # ----------[.TfffTT]----------  Complete overlap 3' longer
                    # ----------[TTfffT.]----------  Complete overlap 5' longer
                    full_dist = max(target_end-feature_start, feature_end-target_start)
                    derived_start = feature_end - full_dist - args.feature_expansion_pad
                    derived_end = feature_start + full_dist + args.feature_expansion_pad
                
                elif (args.feature_expansion_format == 'center_target'):
                    # ----[....TTTT.FFF]----
                    # ----[FFF.TTTT....]----
                    full_dist = max(target_end-feature_start, feature_end-target_start)
                    derived_start = target_end - full_dist - args.feature_expansion_pad
                    derived_end = target_start + full_dist + args.feature_expansion_pad
                
                #elif (args.feature_expansion_format == 'center_both'):
                #    # ----[TTTT...FFF]----
                #    # ----[FFF...TTTT]----
                #    derived_start = min(feature_start, target_start) - args.feature_expansion_pad
                #    derived_end = max(target_end, feature_end) + args.feature_expansion_pad
                
                elif (args.feature_expansion_format == 'justify_feature'):
                    # ----[FFF....TTTTxxx]----
                    # ----[xxxTTTT....FFF]----
                    derived_start = min(feature_start, target_start)
                    derived_end = max(target_end, feature_end)
                    if (derived_start == feature_start):
                        derived_end += args.feature_expansion_pad
                    else:
                        derived_start -= args.feature_expansion_pad
                    
                elif (args.feature_expansion_format == 'justify_target'):
                    # ----[TTTT....FFFxxx]----
                    # ----[xxxFFF....TTTT]----
                    derived_start = min(feature_start, target_start)
                    derived_end = max(target_end, feature_end)
                    if (derived_start == target_start):
                        derived_end += args.feature_expansion_pad
                    else:
                        derived_start -= args.feature_expansion_pad
                
                # By default, 'args.feature_expansion_format' equals 'None', in which case, there should be no
                # extraneous formatting...
                else: # if (args.feature_expansion_format == 'center_both'):
                    # ----[TTTT...FFF]----
                    # ----[FFF...TTTT]----
                    derived_start = min(feature_start, target_start) - args.feature_expansion_pad
                    derived_end = max(target_end, feature_end) + args.feature_expansion_pad
                
                # If derived feature is too small, then we automatically pad it
                # so it reaches the minimum length
                if (derived_end - derived_start < args.feature_expansion_lengths[0]):
                    size_difference = args.feature_expansion_lengths[0] - (derived_end - derived_start)
                    if (args.feature_expansion_format in ['center_feature', 'center_target', 'center_both']):
                        new_pad = (size_difference+1)//2
                        derived_start -= new_pad
                        derived_end += new_pad
                        cls.logger.info("Added {} nt of padding bases to either side of DERIVED FEATURE".format(new_pad))
                    
                    elif (args.feature_expansion_format == 'justify_feature'):
                        if (derived_start == feature_start):
                            derived_end += size_difference
                        else:
                            derived_start -= size_difference
                        cls.logger.info("Added {} nt of padding bases to one side of DERIVED FEATURE".format(size_difference))
                    
                    elif (args.feature_expansion_format == 'justify_target'):
                        if (derived_start == target_start):
                            derived_end += size_difference
                        else:
                            derived_start += size_difference
                        cls.logger.info("Added {} nt of padding bases to one side of DERIVED FEATURE".format(size_difference))
                
                cls.logger.info("  eqi={}, fi={}, derived=({}, {}), len(derived)={}, target=({}, {}), feature={}".format(eqi, fi, derived_start, derived_end, derived_end-derived_start, target_start, target_end, f))
                
                fb.append((fi, derived_start, derived_end))
            
            feature_bounds.append(fb) # Will contain redundant bounds if a Target appears multiple times
        
        # Return bounds
        return feature_bounds
    
    @classmethod
    def expand_for_targets(cls, args, contigs, feature_list):
        '''
        Expand all input features in parallel, the same amounts, taking 'args.feature_expansion_format' into account
        Will expand to all appropriate sizes within 'args.feature_expansion_lengths'.
        :param args: argparse Namespace object
        :param contigs: dict with key=header value=sequence
        :param feature_list: List of all features that make up a homologous group
        :return: 
        '''
        
        # TODO: Check to see if this function is supposed to match allele-specific Target permutations as equivalence groups
        #       If it is supposed to, then make sure it is actually doing this (I forgot if it is).
        #       If it is not supposed to, then make sure it doesn't actually do this.
        #       (For conceptual details, see Jan 2020 presentation)
        
        from . import targets
        
        # 'all'       Expand all Features the same amount so they have the same Target (using args.feature_expansion_format)
        # 'exclusive' Expand all Features the same amount until each has a unique Target (using args.feature_expansion_format)
        # 'any'       Expand all Features the same amount until all Features have a Target (using args.feature_expansion_format)
        
        # Create a matrix that stores which Targets are observed within which Features
        #       F1  F2  F3 ... FN
        #   T1  X   X   -      -
        #   T2  X   X   X      X
        
        targets_list = []
        for f in feature_list:
            
            feature_start = f.start
            feature_end = f.end or len(contigs[f.contig]) # if (f.end == None), then it will be len(contig_sequence)
            
            ex_feature_start = max(0, feature_start-args.feature_expansion_lengths[1])
            ex_feature_end = min(feature_end+args.feature_expansion_lengths[1], len(contigs[f.contig]))
            
            # Identify all Targets within FULL region surrounding each homologous Feature
            t_list = targets.Target.get_targets(args, contigs[f.contig], start=ex_feature_start, end=ex_feature_end) # Does both orientations (+/-)
            targets_list.append(t_list)
        
        # We match all identified Targets to their putative homologs in other Features (to-be) (derived)
        # For instance, the homologous Features A and B below
        #   A  [F..T1..T2...T3a....T4...T3b]
        #   B  [F...T3...T5]
        # Would necessitate the following possible pairing and results
        #   A (T1 ) x B (T3) - different
        #   A (T1 ) x B (T5) - different
        #   A (T2 ) x B (T3) - different
        #   A (T2 ) x B (T5) - different
        #   A (T3a) x B (T3) - equivalent
        #   A (T3a) x B (T5) - different
        #   A (T4 ) x B (T3) - different
        #   A (T4 ) x B (T5) - different
        #   A (T3b) x B (T3) - equivalent
        #   A (T3b) x B (T5) - different
        # Thus there would be two candidates for derived Feature homologs:
        #   1) A  [F...........T3a]
        #      B  [F...T3]
        #   2) A  [F.......................T3b]
        #      B  [F...T3]
        
        
        
        # Need to find out which Targets are present (within 'args.target_distance') in each feature's t_list...
        
        # First, we iterate through ALL Targets, regardless of their Feature origin,
        # and we compare to all Targets
        #flat_targets = [t for tlist in targets_list for t in tlist]
        #flat_target_origins = [i for i, tlist in enumerate(targets_list) for t in tlist]
        #
        #i2t = {}
        #t2fi = {}
        #for i, (fi, t) in enumerate(zip(flat_target_origins, flat_targets)):
        #    i2t[i] = t
        #    t2fi[t] = fi
        
        
        t2fi = {}
        for fi, tlist in enumerate(targets_list):
            for ti, t in enumerate(tlist):
                t2fi[t] = fi
        
        
        new_equivalents = {}
        
        for fi1, tlist1 in enumerate(targets_list):
            for ti1, t1 in enumerate(tlist1):
                ne = new_equivalents.setdefault((fi1, ti1), dict())
                ne.setdefault(fi1, []).append(t1)
                for fi2, tlist2 in enumerate(targets_list):
                    if (fi1 > fi2):
                        for ti2, t2 in enumerate(tlist2):
                            #print(fi1, ti1, fi2, ti2)
                            
                            pattern = nucleotides.build_regex_pattern(
                                t1[5],
                                max_substitutions=args.target_distance[0],
                                max_insertions=args.target_distance[1],
                                max_deletions=args.target_distance[1],
                                max_errors=args.target_distance[2]
                            )
                            m = regex.match(pattern, t2[5], flags=regex.BESTMATCH|regex.IGNORECASE) # Code to allow fuzzy regex matching
                            #m = regex.match(t1[5], t2[5]) # Code for non-fuzzy regex matching
                            if m:
                                ne = new_equivalents.setdefault((fi1, ti1), dict())
                                ne.setdefault(fi2, []).append(t2)
                                
                                ne = new_equivalents.setdefault((fi2, ti2), dict())
                                ne.setdefault(fi1, []).append(t1)
                                
                                # We populate 'Target.equivalents'
                                #   key=seq, value=set(seq, seq, ...)
                                #   # Example: { 'ACA': {'ACA', 'CCA'}, 'CCA': {'ACA', 'CCA'}}
                                # This will only store equivalents if they are on different features (not the same feature)
                                s = targets.Target.equivalents.setdefault(t1[5], set())
                                s.add(t1[5])
                                s.add(t2[5])
                                
                                s = targets.Target.equivalents.setdefault(t2[5], set())
                                s.add(t1[5])
                                s.add(t2[5])
        
        # Make non-redundant list of equivalent Targets.
        # If a Feature would have several equivalent Targets, then each would give a new entry. For example:
        #    F1   [T1]
        #    F2   [T2a..T2b]
        #    F3   [T3a....T3b]
        # Would yield the following equivalents:
        #    [[T1, T2a, T3a], [T1, T2a, T3b], [T1, T2b, T3a], [T1, T2b, T3b]]
        
        better_equivalents = [] # List of list of Targets: [[T1, T2], [T3], [T4, T5, T6]]
        for (fi, ti), fi_target_dict in new_equivalents.items():
            be_list = [[]]
            for fi2, tlist2 in fi_target_dict.items():
                #if (len(tlist2) == 1):
                #    for be in be_list:
                #        be.append(tlist2[0])
                #elif (len(tlist2) > 1):
                be_list = [copy(be) for be in be_list for t in tlist2]
                
                for t, be in zip(cycle(tlist2), be_list):
                    be.append((fi2, t))
            
            for be in be_list:
                sbe = sorted(be)
                if sbe not in better_equivalents:
                    better_equivalents.append(sbe)
        
        spec_dict = {
            'all': 0b1, # multi-allelic
            'exclusive': 0b10, # allele-specific
            'any': 0b100, # allele-agnostic
        }
        
        specificities = []
        
        # Calculate specificities and write equivalents to log
        for i, be in enumerate(better_equivalents):
            specs = 0b0
            
            if (len(be) == len(feature_list)): # multi-allelic
                specs |= spec_dict['all']
            if (len(be) == 1): # allele-specific
                specs |= spec_dict['exclusive']
            specs |= spec_dict['any'] # allele-agnostic
            specificities.append(specs)
            
            cls.logger.info('  Equivalent: {}, ALL={}, EXCLUSIVE={}, ANY={}'.format(i, specs & spec_dict['all'], specs & spec_dict['exclusive'], specs & spec_dict['any']))
            for fi, t in be:
                cls.logger.info('    {}: {}'.format(fi, t))
        cls.logger.info('')
        
        return_list = []
        
        # We go through the 'better_equivalents' list to find ones that match the selected 'args.target_specificity'
        cls.logger.info("Filtering equivalents to match '--target_specificity'...")
        for i, (be, specs) in enumerate(zip(better_equivalents, specificities)):
            if (specs & spec_dict[args.target_specificity]):
                return_list.append(be)
                
                cls.logger.info('  Equivalent: {}'.format(i))
                for fi, t in be:
                    cls.logger.info('    {}: {}'.format(fi, t))
        cls.logger.info('')
        
        # Return equivalents list
        return return_list
    
    def expand_feature(self, args, contig_sequence):
        """
        Expand the feature to ensure it has a viable 'Target' sequence
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
                    self.logger.info("WARNING: Selected feature '{}' overlaps with excluded feature '{}'".format(self.name, exf_obj.name))
        
        self.logger.info('    feature: {}'.format(self.name))
        self.logger.info('     bounds: {}..{}'.format(self.start, self.end))
        self.logger.info('     limits: {}..{}'.format(max_upstream_coord, max_downstream_coord))
        
        
        
        #feature_sequence = contig_sequence[feature_start:feature_end]
        # The '10' is for upstream/downstream adjacent sequences
        #ex_feature_start = max(0, feature_start-(10+args.feature_expansion_lengths[1]))
        #ex_feature_end = min(feature_end+(10+args.feature_expansion_lengths[1]), contig_length)
        ex_feature_start = max(0, feature_start-args.feature_expansion_lengths[1])
        ex_feature_end = min(feature_end+args.feature_expansion_lengths[1], contig_length)
        
        #max_feature_sequence = contig_sequence[ex_feature_start:ex_feature_end]
        
        # We scan for targets only once
        targetsl = targets.Target.get_targets(args, contig_sequence, start=ex_feature_start, end=ex_feature_end) # Does both orientations (+/-)
        self.logger.info('max targets: {}'.format(len(targetsl)))
        
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
            
            if (args.feature_expansion_format == 'center_feature'):
                # ----[........FFF....TTTT]----  Terminal 3'
                # ----[TTTT....FFF........]----  Terminal 5'
                # ----------[..fffTT]----------  Overlapping 3'
                # ----------[TTfff..]----------  Overlapping 5'
                # ----------[.TfffTT]----------  Complete overlap 3' longer
                # ----------[TTfffT.]----------  Complete overlap 5' longer
                full_dist = max(target_end-feature_start, feature_end-target_start)
                derived_start = feature_end - full_dist - args.feature_expansion_pad
                derived_end = feature_start + full_dist + args.feature_expansion_pad
            
            elif (args.feature_expansion_format == 'center_target'):
                # ----[....TTTT.FFF]----
                # ----[FFF.TTTT....]----
                full_dist = max(target_end-feature_start, feature_end-target_start)
                derived_start = target_end - full_dist - args.feature_expansion_pad
                derived_end = target_start + full_dist + args.feature_expansion_pad
            
            elif (args.feature_expansion_format == 'center_both'):
                # ----[TTTT...FFF]----
                # ----[FFF...TTTT]----
                derived_start = min(feature_start, target_start) - args.feature_expansion_pad
                derived_end = max(target_end, feature_end) + args.feature_expansion_pad
            
            elif (args.feature_expansion_format == 'justify_feature'):
                # ----[FFF....TTTTxxx]----
                # ----[xxxTTTT....FFF]----
                derived_start = min(feature_start, target_start)
                derived_end = max(target_end, feature_end)
                if (derived_start == feature_start):
                    derived_end += args.feature_expansion_pad
                else:
                    derived_start -= args.feature_expansion_pad
                
            elif (args.feature_expansion_format == 'justify_target'):
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
                if (args.feature_expansion_format in ['center_feature', 'center_target', 'center_both']):
                    new_pad = (size_difference+1)//2
                    derived_start -= new_pad
                    derived_end += new_pad
                    self.logger.info("Added {} nt of padding bases to either side of DERIVED FEATURE".format(new_pad))
                
                elif (args.feature_expansion_format == 'justify_feature'):
                    if (derived_start == feature_start):
                        derived_end += size_difference
                    else:
                        derived_start -= size_difference
                    self.logger.info("Added {} nt of padding bases to one side of DERIVED FEATURE".format(size_difference))
                
                elif (args.feature_expansion_format == 'justify_target'):
                    if (derived_start == target_start):
                        derived_end += size_difference
                    else:
                        derived_start += size_difference
                    self.logger.info("Added {} nt of padding bases to one side of DERIVED FEATURE".format(size_difference))
            
            self.logger.info("derived_start = {}".format(derived_start))
            self.logger.info("  derived_end = {}".format(derived_end))
            self.logger.info("derived_end - derived_start = {}".format(derived_end - derived_start))
            derived_sets.add((derived_start, derived_end))
        
        # Create a counter for the derived feature
        count = 0
        for derived_start, derived_end in sorted(derived_sets, key=lambda x: x[1]-x[0]): # in order from smallest to largest
            # If derived feature length is within the desired size range, then create the derived feature
            if (args.feature_expansion_lengths[0] <= derived_end - derived_start <= args.feature_expansion_lengths[1]):
                if (max_upstream_coord <= derived_start <= derived_end <= max_downstream_coord):
                    new_name = self.name + '_derived-' + str(count)
                    new_attributes = self.attributes.copy()
                    new_attributes[args.tag] = new_name
                    new_feature = Feature(self.contig, derived_start, derived_end, self.strand, name=new_name, source=self.source, feature_type=self.feature_type, score=self.score, frame=self.frame, attributes=new_attributes, origin=Feature.DERIVED, parent=self, gene=self.gene)
                    Feature.features[new_name] = new_feature
                    count += 1
                    self.logger.info("DERIVED FEATURE '{}' created.".format(new_name))
                else:
                    self.logger.info("DERIVED FEATURE would overlap with excluded feature.")
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
                    self.logger.info("WARNING: Selected feature '{}' overlaps with excluded feature '{}'".format(self.name, exf_obj.name))
        
        self.logger.info('feature: {}'.format(self.name))
        self.logger.info( 'bounds: {}..{}'.format(self.start, self.end))
        self.logger.info(' limits: {}..{}'.format(max_upstream_coord, max_downstream_coord))
        
        # Gradually expand feature size until the minimum number of targets is found
        feature_sequence = contig_sequence[feature_start:feature_end]
        targetsl = targets.Target.get_targets(args, feature_sequence) # Does both orientations (+/-)
        self.logger.info('targets: {}'.format(len(targetsl)))
        count = 0
        while ((len(targetsl) < minimum_targets_per_feature) and (((feature_start, feature_end) != (0, contig_length)) or ((feature_start, feature_end) != (max_upstream_coord, max_downstream_coord)))):
            feature_start = max(0, feature_start-expansion_size, max_upstream_coord)
            feature_end = min(contig_length, feature_end + expansion_size, max_downstream_coord)
            
            feature_sequence = contig_sequence[feature_start:feature_end]
            targetsl = targets.Target.get_targets(args, feature_sequence) # Does both orientations (+/-)
            self.logger.info('targets: {}'.format(len(targetsl)))
            count += 1
        
        # Save the new feature as a Feature object
        if (count > 0):
            new_name = self.name + '_derived'
            new_attributes = self.attributes.copy()
            new_attributes[args.tag] = new_name
            new_feature = Feature(self.contig, feature_start, feature_end, self.strand, name=new_name, source=self.source, feature_type=self.feature_type, score=self.score, frame=self.frame, attributes=new_attributes, origin=Feature.DERIVED, parent=self)
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
    
    @classmethod
    def get_homologous_region_with_acceptable_polymorphism(cls, args, contigs, mismatches=0, gaps=0, errors=0):
        '''
        Get the closest conserved region with minimal polymorphisms
        The 'Feature.homologs' attribute must have been set the features.
        Also, 'Feature.features' contains only selected features.
        
        :param args:
        :param contigs: Dict with key=name, value=sequence
        :param mismatches: Number of permitted mismatches
        :param gaps: Number of permitted gaps
        :param errors: Number of total permitted errors (both mismatches and gaps)
        '''
        ### ### ### See 'expand_feature_for_homology()' for a prior (incomplete) implementation of this ### ### ### 
        
        # The first thing is to get all the sequences that should be aligned using MSA
        from . import aligners
        
        # Get the MSA aligner object
        aligner = None
        for a in aligners.ms_aligners:
            if (a.name == 'mafft'):
                aligner = a
                break
        
        # Get all homologous features
        homologs = set()
        for fname, f in Feature.features.items():
            if (f.origin == Feature.INPUT):
                #contig_sequence = contigs[f.contig]
                
                # We assume 'f.homologs' have all been populated completely and correctly
                for s in f.homologs:
                    #homologs.add(tuple(sorted(s)))
                    homologs.add(tuple(s))
        
        # Alternative way to do this ('feature2gene' is an input argument) (untested)
        #homologs = {}
        #for fname, f in Feature.features.items():
        #    fp = f.get_parent() # short for 'feature_parent'
        #    g = feature2gene[fp.name]
        #    homologs.setdefault(g, set()).add(fp)
        
        # Convert homologs into a list
        homologs = list(homologs)
        
        # Calculate the length of the homology arm to search for
        length = (max(args.excise_donor_lengths)+1)//2
        
        # 'args.feature_expansion_lengths[1]' is the max distance up/down-stream to expand the feature
        # This is the length of the homology region to search within
        
        # Iterate through each homolog group (now tuples)
        for features in homologs:
            gene = features[0].get_gene()
            logging.info('Working on gene: {}'.format(gene))
            logging.info('homolog set: {}'.format(features))
            
            
            for side in ['us', 'ds']:
                logging.info('Working on MSA for side: {}'.format(side))
                us_sequences = []
                us_names = []
                us_contigs = []
                us_starts = []
                us_ends = []
            
                if side in ['us', 'US', 'left', "5'", 'upstream']:
                    # Find the US region of all homologous features
            
                    for f in features:
                        us_start = max(0, f.start-args.feature_expansion_lengths[1]-length)
                        us_end = min(f.start, len(contigs[f.contig]))
                        
                        us_sequences.append(contigs[f.contig][us_start:us_end])
                        us_names.append('{}:{}:{}..{}'.format(f.name, f.contig, us_start, us_end))
                        us_contigs.append(f.contig)
                        us_starts.append(us_start)
                        us_ends.append(us_end)
                        
                elif side in ['ds', 'DS', 'right', "3'", 'downstream']:
                    # Find the DS region of all homologous features

                    for f in features:
                        us_start = max(0, min(f.end, len(contigs[f.contig])))
                        us_end = min(f.end+args.feature_expansion_lengths[1]+length, len(contigs[f.contig]))
                        
                        us_sequences.append(contigs[f.contig][us_start:us_end])
                        us_names.append('{}:{}:{}..{}'.format(f.name, f.contig, us_start, us_end))
                        us_contigs.append(f.contig)
                        us_starts.append(us_start)
                        us_ends.append(us_end)
                
                # Write sequences to FASTA file
                msa_input_path = utils.write_merged_fasta((us_names, us_sequences), os.path.join(args.folder, 'msa-{}-{}-input.fasta'.format(utils.slugify(gene), side)))
                
                # Perform MSA
                msa_output_path = aligner.align(msa_input_path, 'msa-{}-{}-output.fasta'.format(utils.slugify(gene), side), args.folder, args.processors)
                
                # Parse MSA output to get records
                # We assume the ordering of sequences in 'records' mirrors the input sequence ordering
                records = list(aligner.load(msa_output_path))
                
                # Print records to log
                cls.logger.info('records = [')
                for r in records:
                    cls.logger.info('  {}'.format(r))
                cls.logger.info(']')
                
                # Calculate conservation of MSA
                conservation_list, consensus_list, counts_list, freqs_list = cls.calculate_msa_conservation(records, ambiguities=False)
                
                # Log the MSA
                for line in cls.format_msa(records, conservation_list, consensus_list, counts_list):
                    cls.logger.info(line)
                
                if side in ['us', 'US', 'left', "5'", 'upstream']:
                    # We find right-most acceptable homology region
                    #window_iterator = range(len(records[0].sequence)-length, -1, -1)
                    window_iterator = range(len(records[0].sequence)-length+1)[::-1]
                elif side in ['ds', 'DS', 'right', "3'", 'downstream']:
                    # We find left-most acceptable homology region
                    window_iterator = range(len(records[0].sequence)-length+1)
                
                # Identify all regions of acceptable 'length' that fulfill stringency requirements
                # An alternative method:
                #   Extract the longest region with acceptable homology from the MSA, 
                #   then check to see if the region is long enough (>length)
                good_windows = []
                for w in window_iterator: # w is the position in the alignment of the candidate window
                    mismatch_count = 0
                    gap_count = 0
                    error_count = 0
                    window_check_passed = False
                    
                    for i in range(w, w+length):
                        nt_list = [r.sequence[i] for r in records]
                        
                        # Handle mismatches (treats any level of conservation as equivalent)
                        # TODO: Make it so conservation level within mismatch positions actually matters
                        if (len(nt_list)-counts_list[i]-nt_list.count('-') > 0):
                            error_count += 1
                            mismatch_count += 1
                        
                        # Handle gaps
                        if ('-' in nt_list):
                            error_count += 1
                            gap_count += 1
                        
                        if (mismatch_count <= mismatches):
                            mismatch_check_passed = True
                        else:
                            mismatch_check_passed = False
                        
                        if (gap_count <= gaps):
                            gap_check_passed = True
                        else:
                            gap_check_passed = False
                        
                        if (error_count <= errors):
                            error_check_passed = True
                        else:
                            error_check_passed = False
                        
                        if all([mismatch_check_passed, gap_check_passed, error_check_passed]):
                            window_check_passed = True
                        else:
                            window_check_passed = False
                            break
                    
                    # This is for multi-allelic ('all') dDNA
                    # For allele-specific ('exclusive') dDNA, each Feature within the homology group(gene) should have more than the number of permitted mismatches/gaps/errors
                    # For allele-agnostic ('any') dDNA, then no MSA is performed
                    # TODO: Write an if/else statement for selecting between 'all', 'exclusive', and 'any' dDNA homologies based on input parameters
                    #       Also, write implementations for 'exclusive' and 'any'
                    if (args.donor_specificity == 'all'): # Multi-allelic
                        pass
                    elif (args.donor_specificity == 'exclusive'): # Uni-alleleic
                        pass
                    elif (args.donor_specificity == 'any'): # Allele-agnostic
                        pass
                    
                    if window_check_passed:
                        # TODO: The window at (w, w+length) will be the MAXIMUM homology length. Thus, if there is
                        #       a position with the alignment ['A', '-', '-', '-'], then it will be 'A' in the window,
                        #       where the majority/consensus would list it as '-' (to be removed)
                        #       Thus, the window must be filtered to remove majority '-' characters
                        #       BEFORE the window of size 'length' is searched
                        #records, conservation_list, consensus_list, counts_list, freqs_list = self.filter_msa_minority_gaps(records, conservation_list, consensus_list, counts_list, freqs_list)
                        #fixed_window = []
                        #for i in range(w, w+length):
                        #    nt_list = [r.sequence[i] for r in records]
                        #    if (nt_list.count('-') < len(nt_list)/2):
                        #        fixed_window.append(consensus_list[i])
                        
                        WindowRecord = namedtuple('WindowRecord', ['msa_start', 'msa_end', 'msa_sequence', 'contig', 'start', 'end', 'sequence', 'feature', 'mismatches', 'gaps', 'errors'])
                        
                        window_list = []
                        
                        # We need to convert the window's MSA coordinates into genomic coordinates
                        for ri, r in enumerate(records):
                            window_start = us_starts[ri] + w - r.sequence[:w].count('-')
                            window_end = us_starts[ri] + w + length - r.sequence[:w+length].count('-')
                            wr = WindowRecord(
                                msa_start=w,
                                msa_end=w+length,
                                msa_sequence=r.sequence[w:w+length],
                                #msa_consensus=''.join(consensus_list[w:w+length]),
                                contig=us_contigs[ri],
                                start=window_start,
                                end=window_end,
                                sequence=contigs[us_contigs[ri]][window_start:window_end],
                                feature=features[ri].name,
                                mismatches=mismatch_count,
                                gaps=gap_count,
                                errors=error_count,
                            )
                            window_list.append(wr)
                        
                        # Print to log
                        cls.logger.info('window_list = [')
                        for wr in window_list:
                            cls.logger.info('  {}'.format(wr))
                        cls.logger.info(']')
                        
                        cls.logger.info('')
                        
                        # We add all this info to our master list of good windows
                        good_windows.append((''.join(consensus_list[w:w+length]), window_list))
                    
            

        
        
        
        
        
        
        
    @classmethod
    def calculate_msa_conservation(cls, records, ambiguities=True):
        """
        Input 'records' should be a list of tuples: [('header', 'info', 'ACGGAA...'), ...]
        Returns a string with length equal to multiple sequence alignment.
        The characters used to represent the degree of conservation are
          '*'  (100%)        all residues or nucleotides in that column are identical.
          '|'  (>75%)        conserved substitutions have been observed.
          ':'  (>50%)        semi-conserved substitutions have been observed.
          '.'  (>25%)        poor conservation
          ' '  (>0% or gap)  no conservation
        Positions with gaps '-' will not receive a conservation code
        This function handles ambiguous characters in both input records and output consensus.
        """
        amb_to_nts = {
            'a': ['a'],
            'c': ['c'],
            'g': ['g'],
            't': ['t'],
            'r': ['a', 'g'],
            'y': ['c', 't'],
            'm': ['a', 'c'],
            'k': ['g', 't'],
            'w': ['a', 't'],
            's': ['c', 'g'],
            'b': ['c', 'g', 't'],
            'd': ['a', 'g', 't'],
            'h': ['a', 'c', 't'],
            'v': ['a', 'c', 'g'],
            'n': ['a', 'c', 'g', 't'],
            
            'A': ['A'],
            'C': ['C'],
            'G': ['G'],
            'T': ['T'],
            'R': ['A', 'G'],
            'Y': ['C', 'T'],
            'M': ['A', 'C'],
            'K': ['G', 'T'],
            'W': ['A', 'T'],
            'S': ['C', 'G'],
            'B': ['C', 'G', 'T'],
            'D': ['A', 'G', 'T'],
            'H': ['A', 'C', 'T'],
            'V': ['A', 'C', 'G'],
            'N': ['A', 'C', 'G', 'T'],
        }
        
        # Make IUPAC ambiguity where each key is a tuple of valid DNA nt characters,
        # and values are the ambiguous characters
        nts_to_amb = {}
        for amb, nts in amb_to_nts.items():
            for p in permutations(nts):
                nts_to_amb[p] = amb
        
        # Conservation stats
        cstring = []
        cmost_count = []
        cmost_char = []
        cfreq = []
        
        seqs = [r.sequence for r in records]
        for nts in zip(*seqs):
            # We don't count gaps for determining consensus residue
            nogap_nts = [x for x in nts if x != '-']
            
            # If a record has an ambiguous character, then it counts as a fraction of a
            # canonical character. For instance, it at a position, the alignment is ['Y', 'C', 'C'],
            # then the counts should be [('C', 2.5), ('T', 0.5)]
            #ordered_counts = Counter(nogap_nts).most_common()
            ordered_counts = sorted(nucleotides.get_nt_count(nogap_nts, case_sensitive=False).items(), key=lambda x: x[1], reverse=True)
            nt, count = ordered_counts[0]
            
            # 'cmost_char' will return ambiguity codes if there are ties.
            # For instance, if the 'ordered_counts' at a position is:
            #   [('A', 5), ('G', 5), ('C', 2), ('T', 1)]
            # Then 'cmost_char' should be 'R'
            if ambiguities:
                most_common_chars = [x[0] for x in ordered_counts if x[1] == count]
                cmost_char.append(nts_to_amb[tuple(most_common_chars)])
                cmost_count.append(len(most_common_chars)*count)
                cfreq.append(len(most_common_chars)*count/len(nts))
            else:
                cmost_char.append(nt)
                cmost_count.append(count)
                cfreq.append(count/len(nts))
            
            if ('-' in nts):
                cstring.append(' ')
            else:
                if (count == len(nts)):
                    cstring.append('*')
                elif (count/len(nts) >= 0.75):
                    cstring.append('|')
                elif (count/len(nts) >= 0.5):
                    cstring.append(':')
                elif (count/len(nts) >= 0.25):
                    cstring.append('.')
                else:
                    cstring.append(' ')
        
        return cstring, cmost_char, cmost_count, cfreq
    
    @classmethod
    def filter_msa_minority_gaps(cls, records_list, conservation_list, consensus_list, counts_list, freqs_list):
        fixed_records_list = deepcopy(records_list)
        fixed_conservation_list = deepcopy(conservation_list)
        fixed_consensus_list = deepcopy(consensus_list)
        fixed_counts_list = deepcopy(counts_list)
        fixed_freqs_list = deepcopy(freqs_list)
        
        # Loop through each position in the MSA
        for i in range(len(conservation_list)-1, -1, -1): # Count down (instead of up)
            nt_list = [r.sequence[i] for r in records_list]
            
            # If there is a minority-gap character, then do nothing
            # If there is a majority-gap character, then remove that position
            if (nt_list.count('-') > len(nt_list)/2):
                # Remove the majority gap character from each record
                for r in fixed_records_list:
                    r.sequence = r.sequence[:i] + r.sequence[i+1:]
                # Remove from the other lists:
                fixed_conservation_list.pop(i)
                fixed_consensus_list.pop(i)
                fixed_counts_list.pop(i)
                fixed_freqs_list.pop(i)
        
        # Return all filtered lists
        return fixed_records_list, fixed_conservation_list, fixed_consensus_list, fixed_counts_list, fixed_freqs_list
    
    @classmethod
    def format_msa(cls, records, cstring, cmost_char, cmost_count):
        '''
        Prints multiple sequence alignment with conservation, consensus, and counts lines
        For the 'counts' line, 'X' means there are >10, and
        '/' means the consensus residue has a fractional count (e.g. 2.5)
        '''
        lines = []
        header_len = max(len(y) for y in [r.header for r in records]+['conservation', 'consensus', 'counts'])
        for r in records:
            lines.append(('{:>'+str(header_len)+'} {}').format(r.header, r.sequence))
        lines.append(('{:>'+str(header_len)+'} {}').format('conservation', ''.join(cstring)))
        lines.append(('{:>'+str(header_len)+'} {}').format('consensus', ''.join(cmost_char)))
        
        counts_list = []
        for c in cmost_count:
            if (c < 10):
                if c.is_integer():
                    counts_list.append(str(int(c)))
                else:
                    counts_list.append('/')
            else:
                counts_list.append('X')
        lines.append(('{:>'+str(header_len)+'} {}').format('counts', ''.join(counts_list)))
        
        return lines
    
    def __repr__(self):
        labs = ['name', 'gene', 'location']
        vals = [
            self.name,
            self.get_gene(),
            '{}:{}:{}..{}'.format(self.contig, self.strand, self.start, self.end)
        ]
        return self.__class__.__name__ + '(' + ', '.join('='.join(map(str, x)) for x in zip(labs, vals)) + ')'
