#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/feature.py

# Import standard packages
import logging
import itertools
import os

# Import non-standard packages
import regex

# Import included AddTag-specific modules
#from .__init__ import Target
#from . import targets
from . import nucleotides

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
                h_set = set()
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
    
    @classmethod
    def expand_all_features(cls, args, contigs):
        # TODO: Modify this so that features aren't expanded in isolation. Instead, all features of a homologous
        #       group are expanded together.
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
                # All derived features of from the same parent will be in a list together
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
                    verdict = all(((x[0] <= args.max_homology_errors) and (x[1] <= args.max_homology_errors)) for x in v)
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
                for a in aligners.aligners:
                    if (a.name == 'mafft'):
                        aligner = a
                        break
                
                # Align sequences using MSA
                aln_filename = aligner.align(us_query_filename, None, side+'.aln', folder, args.threads)
                
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
                    self.logger.info("Added {} nt of padding bases to either side of DERIVED FEATURE".format(new_pad))
                
                elif (args.feature_expansion_method == 'justify_feature'):
                    if (derived_start == feature_start):
                        derived_end += size_difference
                    else:
                        derived_start -= size_difference
                    self.logger.info("Added {} nt of padding bases to one side of DERIVED FEATURE".format(size_difference))
                
                elif (args.feature_expansion_method == 'justify_target'):
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

    def __repr__(self):
        labs = ['name', 'gene', 'location']
        vals = [
            self.name,
            self.get_gene(),
            '{}:{}:{}..{}'.format(self.contig, self.strand, self.start, self.end)
        ]
        return self.__class__.__name__ + '(' + ', '.join('='.join(map(str, x)) for x in zip(labs, vals)) + ')'
