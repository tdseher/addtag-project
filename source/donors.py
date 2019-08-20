#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/donors.py

# Import standard packages
import os
import time
import logging
from collections import defaultdict

# Import non-standard packages
import regex

# Import included AddTag-specific modules
from . import utils
from . import nucleotides
from . import motifs
from .feature import Feature
from .targets import Target
from .thermodynamics.oligo import Primer, PrimerPair

logger = logging.getLogger(__name__)

class Donor(object):
    prefix = 'Donor'
    sequences = {}
    indices = {}
    
    logger = logger.getChild('Donor')
    
    @classmethod
    def get_contig_dict(cls):
        contigs = {}
        for sequence, obj in cls.sequences.items():
            contigs[obj.name] = sequence
        return contigs
    
    @classmethod
    def generate_fasta(cls, filename, sep=':'):
        """
        Creates a FASTA file
        >id feature:contig:orientation:start1..end1:mAT:start2..end2 ...
        """
        with open(filename, 'w') as flo:
            for sequence, obj in sorted(cls.sequences.items(), key=lambda x: int(x[1].name.split('-')[1])):
            #for sequence, obj in cls.sequences.items():
                #don_entry = tuple(["dDNA-"+str(i), feature, contig, '+', start1, start2, end1, end2, dDNA])
                print(' '.join(['>'+obj.name, 'spacers='+str(len(obj.spacers))] + sorted([obj.format_location(x, sep) for x in obj.locations])), file=flo)
                print(sequence, file=flo)
        cls.logger.info(cls.__name__ + ' dDNA FASTA generated: {!r}'.format(filename))
        if (len(cls.sequences) == 0):
            cls.logger.info(cls.__name__ + ' no sequences written: Possible error!')
        return filename
    
    def search_other_locations(self, contigs):
        """
        Function to see if the [(start1, end1),insert,(start2,end2)] region
        exists (exactly) on other chromosomes or on the same chromosome in
        another position.
        """
        for loc in self.locations:
            feature, contig, orientation, *segments = loc
            if (len(segments) == 1):
                # this has the format [(start,end)]
                start, end = segments[0]
                query = contigs[contig][start:end]
                finds = []
                for cname, cseq in contigs.items():
                    f = cseq.find(query) # Use 'regex' if you want to allow for some errors
                    if (f >= 0):
                        finds.append((cname, f, f+len(query)))
                    # Also need to search the reverse complement
                ############ NEED TO FINISH THIS PART #############
                pass
            elif (len(segments) == 3):
                # this has the format [(start1, end1),insert,(start2,end2)]
                
                query1 = contigs[contig][segments[0][0]:segments[0][1]]
                query2 = contigs[contig][segments[2][0]:segments[2][1]]
                
                
                
                ############ NEED TO FINISH THIS PART #############
                pass
        
        ############ NEED TO FINISH THIS PART #############
        return None
    
    def format_location(self, location, sep=':'):
        feature, contig, orientation, *segments = location
        output_list = [feature, contig, orientation]
        for x in segments:
            if isinstance(x, str):
                x = (x,)
            output_list.append('..'.join(map(str, x)))
        return sep.join(output_list)
    
    def is_allele_specific(self):
        """
        Returns 'True' if there is only 1 'Feature' in 'self.locations'
        Otherwise, returns 'False'.
        """
        if (len(self.get_features()) == 1):
            return True
        else:
            return False
    
    def is_allele_agnostic(self, homologs):
        """
        Returns 'True' if all homolog features are within 'self.locations'.
        Otherwise, returns 'False'.
        """
        f_set = set(self.get_features())
        h_set = homologs[self.locations[0][0]]
        
        if (len(h_set.difference(f_set)) == 0):
            return True
        else:
            return False
    
    def get_features(self):
        """Return (sorted) list of all feature names this Donor maps to"""
        return sorted(set(x[0] for x in self.locations))
    
    def get_location_features(self):
        """Return (sorted) list of all feature names this Donor maps to"""
        return sorted(set(x[0] for x in self.locations))
    
    def __init__(self, feature, contig, orientation, sequence, *segments, spacer=None):
    #def __init__(self, feature, contig, orientation, segment1, insert, segment2, sequence, spacer=None):
        # List of all genomic locations and features this dDNA corresponds to
        location = (feature, contig, orientation, *segments)
        self.locations = set()
        
        #self.contig = contig, # The name of the contig (not the sequence)
        #self.orientation = orientation
        #self.segment1 = segment1 # (start1, end1)
        #self.insert = insert
        #self.segment2 = segment2 # (start2, end2)
        self.sequence = sequence
        
        # List to hold all spacer objects that map to this dDNA
        self.spacers = set()
        
        # Get the index number
        self.index = len(self.sequences)
        self.name = self.prefix + '-' + str(self.index)
        
        # Add this to the non-redundant list of excision dDNAs if it doesn't already exist
        # otherwise, add this locus to the pre-existing one
        the_dDNA = self.sequences.setdefault(self.sequence, self)
        if spacer:
            the_dDNA.spacers.add(spacer)
        the_dDNA.locations.add(location)
        self.indices.setdefault(the_dDNA.name, the_dDNA)
        
    # >dDNA-0 C2_10210C_B:Ca22chr2B_C_albicans_SC5314:+:2096471..2096521:CCA:2097397:2097444
#    def build_sequence(self, contig_sequence, orientation, segment1, segment2, insert):
#        """Takes substrings of contig to generate the sequence"""
#        sequence = contig_sequence[segment1[0]:segment1[1]] + insert + contig_sequence[segment2[0]:segment2[1]]
#        if (orientation == '-'):
#            print(self.name, "has '-' orientation", file=sys.stderr)
#        else:
#            pass
#        return sequence
    
    #def __repr__(self):
    #    return self.__class__.__name__ + '(' + ' '.join([
    #        self.name,
    #        self.contig + ':' + self.orientation + ':' +
    #        str(self.segment1[0]) + '..' + str(self.segment1[1]) + ':' +
    #        self.insert + ':' + str(self.segment2[0]) + '..' + str(self.segment2[1]),
    #        'features=' + ','.join(self.features) or 'None',
    #        'spacers=' + str(len(self.spacers))
    #        ]) + ')'

    def __repr__(self):
        """Return a string containing a printable representation of the Target object."""
        return self.__class__.__name__ + '(' + ' '.join(
            [self.name, 'spacers='+str(len(self.spacers))] +
            [self.format_location(x) for x in sorted(self.locations)]
            ) + ')'

class ExcisionDonor(Donor):
    prefix = 'exDonor'
    sequences = {}
    indices = {}
    
    logger = logger.getChild('ExcisionDonor')
    
    # @classmethod
    # def get_targets(cls, args, sequence):
    #     targets = set()
    #     #for seq_i, sequence in enumerate(dDNAs):
    #     for orientation in ['+', '-']:
    #         if (orientation == '-'):
    #             sequence = nucleotides.rc(sequence)
    #         
    #         #for i in range(len(args.parsed_motifs)):
    #         for mymotif in OnTargetMotif.motifs:
    #             #spacers, pams, side = args.parsed_motifs[i]
    #             spacers, pams, side = mymotif.parsed_list
    #             #compiled_regex = args.compiled_motifs[i]
    #             #matches = nucleotides.motif_search(sequence, spacers, pams, side)
    #             matches = nucleotides.motif_search2(sequence, side, mymotif.compiled_regex)
    #             for seq, start, end, spacer, pam in matches:
    #                 if (orientation == '-'):
    #                     start, end = len(sequence) - end, len(sequence) - start
    #                 t_upstream = sequence[start-10:start]
    #                 t_downstream = sequence[end:end+10]
    #                 filtered_targets = target_filter(seq, spacer, pam, t_upstream, t_downstream, args)
    #                 for filt_seq, filt_spacer, filt_pam in filtered_targets:
    #                     targets.add((orientation, start, end, t_upstream, t_downstream, filt_seq, side, filt_spacer, filt_pam, mymotif.motif_string, tuple([tuple(x) if isinstance(x, list) else x for x in mymotif.parsed_list])))
    #     return sorted(targets) # becomes a list
    
    @classmethod
    def mintag_exhaustive_site_search(cls, args, f, contig_sequence, orientation):
        """
        Method for generating all potential mAT given the us/ds trims and
        insert size. Also called "Brute force" method
        
        This code is run once for each feature
        """
        for us_trim in range(args.excise_upstream_feature_trim[0], args.excise_upstream_feature_trim[1]+1):
            for ds_trim in range(args.excise_downstream_feature_trim[0], args.excise_downstream_feature_trim[1]+1):
                for us_hom in range(args.excise_upstream_homology[0], args.excise_upstream_homology[1]+1):
                    for ds_hom in range(args.excise_downstream_homology[0], args.excise_downstream_homology[1]+1):
                        for insert_length in range(args.excise_insert_lengths[0], args.excise_insert_lengths[1]+1):
                            if (args.excise_donor_lengths[0] <= us_hom+insert_length+ds_hom <= args.excise_donor_lengths[1]):
                                start1, end1 = f.start-us_hom-us_trim, f.start-us_trim
                                start2, end2 = f.end+ds_trim, f.end+ds_hom+ds_trim
                                upstream = contig_sequence[start1:end1]
                                downstream = contig_sequence[start2:end2]
                                #upstream = contigs[contig][start - args.excise_donor_homology[1]:start]
                                #downstream = contigs[contig][end:end + args.excise_donor_homology[1]]
                                
                                # when insert_length = 0, then the kmers are [''] (single element, empty string)
                                for mAT in nucleotides.kmers(insert_length):
                                    # Add this candidate dDNA to the list of all candidate dDNAs
                                    dDNA = upstream + mAT + downstream
                                    
                                    #if (args.excise_donor_lengths[0] <= len(dDNA) <= args.excise_donor_lengths[1]):
                                    #dDNAs.append(dDNA)
                                    new_targets = Target.get_targets(args, dDNA) # [(orientation, start, end, filt_seq, side, filt_spacer, filt_pam), ...]
                                    
                                    # The ReversionTarget must overlap the junction. These are valid cases:
                                    # dDNA    uuuuuuuuuuuuiiiidddddddddddd
                                    # valid      ----------
                                    # valid             ---------
                                    # valid                  ------
                                    # invalid   ----------
                                    # invalid                 -------
                                    # dDNA    uuuuuuuuuuuudddddddddddd
                                    # valid        --------
                                    # valid              --------
                                    
                                    for t in new_targets:
                                        # Target is a tuple: (orientation, start, end, upstream, downstream, sequence, side, spacer, pam, motif)
                                        # If it is not one of the two failure states, then add it (that is, it MUST overlap the junction)
                                        if not ((t[2] < len(upstream)) or (t[1] >= len(upstream) + len(mAT))):
                                            cls(f.name, f.contig, orientation, dDNA, (start1, end1), mAT, (start2, end2), spacer=t)
    
    def heuristic_site_search(self, args, features, contigs):
        """
        Method for intelligently searching a subset of potential mAT in order
        to maximize chance of finding something useful with the minimum
        search space.
        
        This tries to minimize the number of calculations to find a useful dDNA.
        """
        # Automatically expand the search space by increasing the following:
        #  - upstream trim
        #  - mAT insert size
        #  - downstream trim
        # until the time it spent looking exceeds the time threshold
        weight_threshold = 0.95
        best_weight = 0
        auto_expand = True
        kmer_size = 2
        
        f_contig, f_start, f_end, f_strand = features[feature]
        c_contig = contigs[contig]
        c_strand = '+'
        
        if auto_expand:
            # If the best reTarget calculated is below the threshold, then
            # increase the mAT length, and calculate again.
            #
            # ALSO: If the upstream homology region is allele-specific
            #       or the downstream homology region is allele-specific
            #       Then expand the feature to overcome the specificity, thereby
            #       removing polymorphisms
            while (best_weight < weight_threshold):
                # ===upstream=homology===/us-trim/insert/ds-trim===downstream=homology
                #        CCAAACC ACGGAACGAC GAA AACGACGAGG GGCTAGAGACTAG
                #        -ignore -------spacer---------PAM ignore-------
                
                # Generate proposed spacer containing the junction
                # Generate entire dDNA sequence, placing the spacer in the center
                start1, end1 = f_start-us_hom-us_trim, f_start-us_trim
                start2, end2 = f_end+ds_trim, f_end+ds_hom+ds_trim
                upstream = my_contig[start1:end1]
                downstream = my_contig[start2:end2]
                dDNA = upstream + mAT + downstream
                
                # Reference list of all possible kmers
                kmers = list(nucleotides.permute_genotypes(["ACGT"]*kmer_size)) # ['AA', 'AC', ... 'TT']
                kmers = dict(zip(kmers, [0]*len(kmers))) # {'AA': 0, 'AC': 0, ..., 'TT': 0}
                
                seq = 'AATATGGCTCGATGAGATCTCGACTAGTGC'
                # Count all kmers for proposed spacer
                for i in range(len(seq)-(kmer_size-1)):
                    kmers[seq[i:i+kmer_size]] += 1
                
                sorted(a, key=lambda x:a[x]) # Evaluate in alphabetical order for ties
                sorted(a, key=lambda x:(a[x], random.random())) # Evaluate in random order for ties
    
    @classmethod
    def bartag(cls, args, f, contig_sequence, orientation, bartags):
        """
        Method for generating dDNA with bartags.
        Does not check for ki-gRNA target sequences.
        
        This code is run once for each feature
        Does not respect:
          args.excise_upstream_homology
          args.excise_downstream_homology
          args.excise_upstream_feature_trim
          args.excise_downstream_feature_trim
          args.excise_insert_lengths
          args.excise_donor_lengths
        """
        
        # Use this instead of 'args.excise_donor_lengths'
        #slen = 100
        slen = min(args.excise_donor_lengths)
        
        # when insert_length = 0, then the kmers are [''] (single element, empty string)
        for bartag in bartags:
            # Add this bartag dDNA substring to the center of the up/down-stream DNA
            # dDNA = upstream[-slen//2-(-len(bartag)//2):] + bartag + downstream[:slen//2-len(bartag)//2]
            # start1=15, start2=5
            # dDNA = sequence[15+(-slen//2-(-len(bartag)//2)):15] + bartag + sequence[5:5+(slen//2-len(bartag)//2)]
            us_start, us_end = f.start+(-slen//2-(-len(bartag)//2)), f.start
            ds_start, ds_end = f.end, f.end+(slen//2-len(bartag)//2)
            upstream = contig_sequence[us_start:us_end]
            downstream = contig_sequence[ds_start:ds_end]
            dDNA = upstream + bartag + downstream
            
            # We don't calculate targets
            #new_targets = Target.get_targets(args, dDNA) # [(orientation, start, end, filt_seq, side, filt_spacer, filt_pam), ...]
            
            # We create the new ExcisionDonor object (with no spacers)
            cls(f.name, f.contig, orientation, dDNA, (us_start, us_end), bartag, (ds_start, ds_end), spacer=None)
    
    @classmethod
    def addtag(cls, args, f, contig_sequence, orientation):
        """
        Method for generating dDNAs with complete SPACER+PAM targeting sites
        with flanking homology arms
        """
        # Use this instead of 'args.excise_donor_lengths'
        #slen = 100
        slen = min(args.excise_donor_lengths)
        
        genome_composition = nucleotides.get_seq_dist(contig_sequence, 1, 4) # step_size=1, kmer_size=4
        
        for m in motifs.OnTargetMotif.motifs: # <------ Do I need to add OffTargetMotif.motifs here as well?
            addtag_seqs = [] # holds at most, 10x1000= 10000 sequences
            
            # Keep generating and testing addtags until the score raises above the arbitrary threshold of 0.9
            # Or until the maximum number of batches have been run
            max_batches = 10
            current_batch = 0
            best_score = 0
            while ((best_score < 0.9) and (current_batch < max_batches)):
                current_batch+= 1
                # Generate batches of 1000 addtags:
                addtag_batch = []
                for i in range(1000):
                    # Generate a pseudorandom motif
                    addtag_batch.append(m.generate_sequence(compositions=genome_composition, complement=True, samples=100))
                
                # Write the addtag sequences to a file, and align them to the genome
                #### This happens in 'main()' ####
                
                # Calculate the scores of the addtag sequences
                #### This happens in 'main()' ####
                
                # Add the batched oligos to the list of all generated ones
                addtag_seqs += addtag_batch
                
                best_score = 0
            
            # Only keep the top N best (i.e. 100)
            #for addtag in sorted(addtag_seqs, key=lambda x: x.score, reverse=True)[:100]:
            # Since we haven't scored these yet, then let's just keep all of them
            for addtag in addtag_seqs:
                # Get the contig locations of the homology regions
                us_start, us_end = f.start+(-slen//2-(-len(addtag)//2)), f.start
                ds_start, ds_end = f.end, f.end+(slen//2-len(addtag)//2)
                
                # Get the sequences for the homology regions
                upstream = contig_sequence[us_start:us_end]
                downstream = contig_sequence[ds_start:ds_end]
                
                # Stitch together the homology regions with the addtag (exogenous gRNA target)
                dDNA = upstream + addtag + downstream
                
                ##### Extra condition: ADDTAG must NOT be present in any of the ki-dDNAs #####
                
                # The target must match the 'addtag' sequence completely
                targets = Target.get_targets(args, dDNA, start=len(upstream), end=len(upstream)+len(addtag)) # [(orientation, start, end, filt_seq, side, filt_spacer, filt_pam), ...]
                
                for t in targets:
                    #                     0            1      2    3         4           5         6     7       8    9
                    # Target is a tuple: (orientation, start, end, upstream, downstream, sequence, side, spacer, pam, motif)
                    # The target must match the addtag start and end positions within the dDNA exactly
                    if ((t[1] == len(upstream)) and (t[2] == len(upstream)+len(addtag))):
                        # Add this candidate dDNA to the list of all candidate dDNAs
                        cls(f.name, f.contig, orientation, dDNA, (us_start, us_end), addtag, (ds_start, ds_end), spacer=t)
        #### End 'addtag()' ####
    
    @classmethod
    def unitag(cls, args, f, contig_sequence, orientation, unitag):
        """
        Method for generating dDNAs with complete SPACER+PAM targeting sites
        with flanking homology arms
        
        This function assumes the input 'unitag' has already been calculated
        """
        # Use this instead of 'args.excise_donor_lengths'
        #slen = 100
        slen = min(args.excise_donor_lengths)
        
        # Get the contig locations of the homology regions
        us_start, us_end = f.start+(-slen//2-(-len(unitag)//2)), f.start
        ds_start, ds_end = f.end, f.end+(slen//2-len(unitag)//2)
        
        # Get the sequences for the homology regions
        upstream = contig_sequence[us_start:us_end]
        downstream = contig_sequence[ds_start:ds_end]
        
        # Stitch together the homology regions with the addtag (exogenous gRNA target)
        dDNA = upstream + unitag + downstream
        
        ##### Extra condition: UNITAG must NOT be present in any of the ki-dDNAs #####
        
        # The target must match the 'addtag' sequence completely
        # This will search all OnTargetMotif motifs
        targets = Target.get_targets(args, dDNA, start=len(upstream), end=len(upstream)+len(unitag)) # [(orientation, start, end, filt_seq, side, filt_spacer, filt_pam), ...]
        
        for t in targets:
            #                     0            1      2    3         4           5         6     7       8    9
            # Target is a tuple: (orientation, start, end, upstream, downstream, sequence, side, spacer, pam, motif)
            # The target must match the addtag start and end positions within the dDNA exactly
            if ((t[1] == len(upstream)) and (t[2] == len(upstream)+len(unitag))):
                # Add this candidate dDNA to the list of all candidate dDNAs
                cls(f.name, f.contig, orientation, dDNA, (us_start, us_end), unitag, (ds_start, ds_end), spacer=t)
        #### End 'unitag()' ####
    
    @classmethod
    def generate_donors(cls, args, contigs, feature2gene):
        """
        Creates the DNA oligo with the structure:
        [upstream homology][unique gRNA][downstream homology]
        that excises the target feature
        """
        
        # Behave differently depending on user input (mintag/addtag/unitag/bartag)
        if (args.ko_dDNA == 'mintag'):
            # Generate the full set of potential dDNAs
            for feature_name, f in Feature.features.items():
            #for feature in features:
                #contig, start, end, strand = features[feature]
                # start & end are 0-based indices, inclusive/exclusive
                
                # assumes start < end
                # DNA 5' of feature is upstream
                # DNA 3' of feature is downstream
                
                contig_sequence = contigs[f.contig]
                orientation = '+' # ????? what is this for? I don't remember
                
                # For each potential dDNA, evaluate how good it is
                #cls.mintag_exhaustive_site_search(args, feature, contig, start, end, strand, contig_sequence, orientation)
                cls.mintag_exhaustive_site_search(args, f, contig_sequence, orientation)
            
        elif (args.ko_dDNA == 'addtag'): ######################## NEED TO ADD FLANKTAG PROCESSING ########################
            # Generate one unique tag for each feature (same take for homologs)
            # using the complement of the DNA sequence composition for the organism
            
            for feature_name, f in Feature.features.items():
                contig_sequence = contigs[f.contig]
                orientation = '+' # I forgot what this is for
                cls.addtag(args, f, contig_sequence, orientation)
            
        elif (args.ko_dDNA == 'unitag'): ######################## NEED TO ADD FLANKTAG PROCESSING ########################
            # Generate a single tag that is the same for all features
            
            # Generate a random tag (that is the complement sequence composition)
            # Evaluate all the tags
            # Pick the best tag
            # Add that best tag to all dDNAs
            
            # Pick an arbitrary motif (this is just the first-defined one)
            m = motifs.OnTargetMotif.motifs[0] # <------ Do I need to add OffTargetMotif.motifs here as well?
            
            # Create a single, random unitag according to the chosen motif (for testing)
            unitag = m.generate_sequence(compositions=genome_composition, complement=True, samples=100)
            
            # Create the dDNAs using the chosen unitag
            for feature_name, f in Feature.features.items():
                contig_sequence = contigs[f.contig]
                orientation = '+' # I forgot what this is for
                cls.unitag(args, f, contig_sequence, orientation, unitag)
        
        elif (args.ko_dDNA == 'bartag'): ######################## NEED TO ADD FLANKTAG PROCESSING ########################
            # if number of features x number of bartags > 200, then exit with warning
            if (args.bartag_number * len(Feature.features) > 200):
                raise Exception("Too many features or bartags. Must be <= 200.\nYou have: features × bartags = " + str(len(Feature.features)) + " × " + str(args.bartag_number) + " = " + str(args.bartag_number * len(Feature.features)))
            
            # Generate up to 200 barcodes
            cls.logger.info('Calculating barcodes...')
            barcodes, barcode_fail_n, barcode_success_n = bartag.generate_min_distance(args.bartag_motif, args.bartag_distance, max_successes=200)
            cls.logger.info('Found {} barcodes: {} fails, {} successes'.format(len(barcodes), barcode_fail_n, barcode_success_n))
            cls.logger.info('List of barcodes:')
            for bc in barcodes:
                cls.logger.info('  {}'.format(bc))
            
            # Assign a barcode to each feature,
            # Generate the dDNAs
            bartag_assignments = {}
            for feature_name, f in Feature.features.items():
                # Add a barcodes equal to the number specified via the command line
                bartag_assignment[feature_name] = [barcodes.pop() for i in range(args.bartag_number)]
                
                contig_sequence = contigs[f.contig]
                orientation = '+' # I forgot what this is for
                cls.bartag(args, f, contig_sequence, orientation, bartag_assignment[feature_name])
        
        else: ######################## NEED TO ADD FLANKTAG PROCESSING ########################
            # 'args.ko_dDNA' is the path to a FASTA file
            cls.logger.info('Processing ko-dDNA with FASTA as input.')
            
            # Load the FASTA file into memory
            ko_contigs = utils.old_load_fasta_file(args.ko_dDNA)
            
            if (len(ko_contigs) == 1):
                # If only a single sequence appears in FASTA file, then treat it as a unitag
                tag = list(ko_contigs.values())[0]
            
                for feature_name, f in Feature.features.items():
                    contig_sequence = contigs[f.contig]
                    orientation = '+' # I forgot what this is for
                    cls.unitag(args, f, contig_sequence, '+', tag)
                    cls.unitag(args, f, contig_sequence, '+', nucleotides.rc(tag))
            else:
                # If there are multiple sequences in FASTA file, then
                # cross-reference their primary sequence headers with
                # the (parent) feature names and gene names (from homologs file)
                for feature_name, f in Feature.features.items():
                    fp = f.get_parent() # short for 'feature_parent'
                    try:
                        if fp.name in ko_contigs:
                            tag = ko_contigs[fp.name]
                        else:
                            tag = ko_contigs[feature2gene[fp.name]]
                    except KeyError:
                        raise Exception("The ko-dDNA file '{}' has no sequence with a primary header that matches '{}' or '{}'".format(args.ko_dDNA, fp.name, feature2gene[fp.name]))
                    
                    contig_sequence = contigs[f.contig]
                    orientation = '+' # I forgot what this is for
                    cls.unitag(args, f, contig_sequence, '+', tag) # Insert TAG in '+' orientation
                    cls.unitag(args, f, contig_sequence, '+', nucleotides.rc(tag)) # Insert TAG in '-' orientation
        
        #### End 'generate_donors()' ####
    
    @classmethod
    def generate_alignments(cls):
        # Eventually, this function should return something like the following:
        #                                                                  ┌insert
        #                 ┌──────────────────ds homology─────────────────┐┌┴┐┌──────────────────us homology──────────────────┐
        # exDonor-42 dDNA ACCATAACGTTTACTTGTTTAATATGCTATTGATATCTATATTTTTTTCCCTATGTGTAGTGCTTGTATATGCGTGTGTGATGAGAATAAGATGAATAGA
        # pam<spacer gRNA                                                 CCCTATGTGTAGTGCTTGTATAT
        for ind, obj in cls.indices.items():
            cls.logger.info(obj.name)
            segment_string = ''
            loc = next(iter(obj.locations)) # Pull an arbitrary location record
            for segment in loc[3:]:
                if isinstance(segment, str):
                    length = len(segment)
                else:
                    length = segment[1] - segment[0]
                segment_string += cls.make_label(length, '')
            cls.logger.info(segment_string)
            cls.logger.info(obj.sequence)
            for s in obj.spacers:
                orientation = s[0]
                start = s[1]
                end = s[2]
                spacer = s[7]
                pam = s[8]
                if (orientation == '+'):
                    cls.logger.info(' '*start + spacer + pam)
                else:
                    cls.logger.info(' '*start + nucleotides.rc(spacer+pam))
    
    @staticmethod
    def make_label(length, label):
        """Helper function for 'generate_alignments()'"""
        if (length == 0):
            return ''
        if (length == 1):
            return '╥'
        if (length > 1):
            out = ['─'] * length
            out[0] = '┌'
            out[-1] = '┐'
            if (length >= len(label) + 4):
                start = int(length/2 - len(label)/2)
                for i, c in enumerate(label):
                    out[start+i] = c
            return ''.join(out)
    
    def get_inserts(self):
        """
        Returns a list constructed of the mAT inserts for locations
        """
        return [x[4] for x in self.locations]
    
    def get_trims(self):
        """
        Get the length of the upstream and downstream trims for each location
        """
        trims = []
        for loc in self.locations:
            #contig, start, end, strand = features[loc[0]]
            f = Feature.features[loc[0]]
            #mAT = loc[4]
            left_trim = f.start - loc[3][1]
            right_trim = loc[5][0] - f.end
            #trims.append(len(loc[4]) + loc[5][0]-l[3][1])
            trims.append((left_trim, right_trim))
        return trims
    
    def get_inserts_and_trims(self):
        """
        Returns list of insert size, and us/ds trims with the following format:
          [(insert, us-trim, ds-trim), ...]
        """
        return_list = []
        for loc in self.locations:
            f = Feature.features[loc[0]]
            #contig, start, end, strand = features[loc[0]]
            mAT = loc[4]
            left_trim = f.start - loc[3][1]
            right_trim = loc[5][0] - f.end
            return_list.append((mAT, left_trim, right_trim))
        return return_list

class ReversionDonor(Donor):
    prefix = 'reDonor'
    sequences = {}
    indices = {}
    
    logger = logger.getChild('ReversionDonor')
    
    @classmethod
    def generate_donors(cls, args, contigs, feature2gene):
        """
        Creates the DNA oligo with the structure:
        [amp-F primer]                                                               [amp-R primer]
        [upstream homology][us expanded feature][feature][ds expanded feature][downstream homology]
        """
        
        cls.logger.info("Generating 'ReversionDonor' objects")
        
        if (args.revert_amplification_primers):
            #amplicon_size = (2*min(args.revert_homology_length), 2*max(args.revert_homology_length))
            
            #subset_size = 1000 # 500 # Temporarily limit number of Primer objects used as input to the pair() function
            pp_subset_size = 1000 # Temporarily limit the number of PrimerPairs that are used to create ReversionDonor objects
            #temp_folder = os.path.join(args.temp_folder, 'addtag', os.path.basename(args.folder))
            
            # Make a folder to store the *.gb Genbank flat files
            os.makedirs(os.path.join(args.folder, 'reversion-gb'), exist_ok=True)
            
            # Make list of genes, and the expected input feature numbers
            gene_dict = defaultdict(int)
            
        #    us_gene_p_list_dict = {}
        #    ds_gene_p_list_dict = {}
            for feature_name, f in Feature.features.items():
                fp = f.get_parent() # short for 'feature_parent'
                g = feature2gene[fp.name]
                gene_dict[g] += 1 # Count how many times a feature maps to this gene in only the desired data
                
            #    # Make defaultdicts to count homologous Primer objects
            #    us_gene_p_list_dict[g] = defaultdict(list)
            #    ds_gene_p_list_dict[g] = defaultdict(list)
            
            # There should be at least one ReversionDonor for each feature
            for feature_name, f in Feature.features.items():
                cls.logger.info("Scanning feature '{}:{}:{}:{}..{}' for amplification primers".format(feature_name, f.contig, f.strand, f.start, f.end))
                
                my_contig = contigs[f.contig]
                
                #feature_length = f.end - f.start
                #feature_sequence = my_contig[f.start:f.end]
                
                fp = f.get_parent() # short for 'feature_parent'
                g = feature2gene[fp.name]
                
                us_extend_seq = my_contig[f.start:fp.start]
                ds_extend_seq = my_contig[fp.end:f.end]
                mid_feature_seq = my_contig[fp.start:fp.end]
                if isinstance(args.ki_dDNA, str):
                    # If a file is specified, then it has the knock-in DNA that should be stitched to flanking homology arms
                    #o_feature_name = None
                    #for f2g_f, f2g_g in feature2gene.items():
                    #    if ((fp.name == f2g_f) or (fp.name == f2g_g)):
                    #        o_feature_name = f2g_g
                    #        break
                    
                    ki_contigs = utils.old_load_fasta_file(args.ki_dDNA)
                    
                    try:
                        if fp.name in ki_contigs:
                            mid_feature_seq = ki_contigs[fp.name]
                        else:
                            mid_feature_seq = ki_contigs[feature2gene[fp.name]]
                    except KeyError:
                        raise Exception("The ki-dDNA file '{}' has no sequence with a primary header that matches '{}' or '{}'".format(args.ki_dDNA, fp.name, feature2gene[fp.name]))
                else:
                    # Otherwise, generate KI dDNA that are wild type
                    pass
                
                orientation = '+'
                
                # Won't work if a region spans the start/end of a circular plasmid:
                #   RRRFFFFFF------------RRR
                # Also won't work if the Feature is too close to the beginning of the contig
                region_F_start, region_F_stop = f.start-max(args.revert_homology_length), f.start-min(args.revert_homology_length)
                region_F_start, region_F_stop = max(0, region_F_start), max(0, region_F_stop)
                region_R_start, region_R_stop = f.end+min(args.revert_homology_length), f.end+max(args.revert_homology_length)
                
                region_F = my_contig[region_F_start:region_F_stop]
                region_R = my_contig[region_R_start:region_R_stop]
                
                AMP_F = 10
                AMP_R = 11
                
                cls.logger.info('Adding upstream_F primers to queue...')
                Primer.scan(my_contig, gene=g, locus=fp.name, genome=0, region=AMP_F, contig=f.contig, orientation='+', start=region_F_start, end=region_F_stop, name='AmpF', primer_size=(19,36))
                cls.logger.info('Adding downstream_R primers to queue...')
                Primer.scan(my_contig, gene=g, locus=fp.name, genome=0, region=AMP_R, contig=f.contig, orientation='-', start=region_R_start, end=region_R_stop, name='AmpR', primer_size=(19,36))
            
            cls.logger.info("Total 'Primer' objects: {}".format(len(Primer.sequences)))
            
            # Make simple dict to lookup the feature size given the gene/locus/contig
            gg2feature_size = {}
            for feature_name, f in Feature.features.items():
                f_gene = f.get_gene()
                f_locus = f.get_parent().name
                f_genome = 0
                f_contig = f.contig
                
                f_start, f_end = f.start, f.end
                
                gg2feature_size.setdefault((f_gene, f_locus, f_genome, f_contig), list()).append(f_end-f_start)
            
            cls.logger.info("gene-to-feature:")
            for k, v in gg2feature_size.items():
                cls.logger.info("  {} {}".format(k, v))
            
            # Do required stuff
            for g, g_count in gene_dict.items():
                ampF_list = []
                ampR_list = []
                for pi, (seq, p) in enumerate(Primer.sequences.items()):
                    if any(loc[3] == AMP_F for loc in p.locations):
                        ampF_list.append(p)
                    if any(loc[3] == AMP_R for loc in p.locations):
                        ampR_list.append(p)
                
                feature_sizes = []
                for fname, f in Feature.features.items():
                    if (f.get_gene() == g):
                        feature_sizes.append(f.end-f.start)
                
                
                cls.logger.info('Using fixed primer calculation cutoffs')
                cutoffs = {
                    'length': (19, 28),
                    'last5gc_count': (1, 3),
                    'gc_clamp_length': (1, 2),
                    'gc': (0.4, 0.6),
                    'max_run_length': 4,
                    'max_3prime_complementation_length': 3,
                    'min_delta_g': -4.0,
                    'tm': (52, 65),
                    'max_tm_difference': 2.5,
                    'amplicon_size_range': (2*min(args.revert_homology_length)+min(feature_sizes), 2*max(args.revert_homology_length)+max(feature_sizes)),
                }
                # Add the invariant cutoff parameters
                cutoffs['o_oligo'] = args.selected_oligo
                cutoffs['folder'] = os.path.join(args.temp_folder, 'addtag', os.path.basename(args.folder))
                
                # Perform primer calculations
                cls.logger.info('Performing calculations on primers...')
                for alist in [ampF_list, ampR_list]:
                    start_time = time.time()
                    time_expired = False
                    for pi, p in enumerate(alist):
                        if (pi % 1000 == 0):
                            if ((time.time()-start_time) > args.primer_scan_limit):
                                time_expired = True
                    
                        if not time_expired:
                            cpass = p.summarize(p.checks)
                            if ((p.checks[0] == None) or (not cpass)):
                                p.progressive_check(cutoffs)
                
                ##### Some debug code #####
                nnn = 0
                ttt = 0
                for seq, p in Primer.sequences.items():
                    ttt += 1
                    if ((p.checks[0] != None) and p.summarize(p.checks)):
                        nnn += 1
                cls.logger.info("  Summary of all 'Primer' objects:")
                cls.logger.info("    Number primers with all checks passed = {}".format(nnn))
                cls.logger.info("    Number primers in 'Primer.sequences' = {}".format(ttt))
                ##### Some debug code #####    
                
                #desired_p1_loc_set = set()
                #desired_p2_loc_set = set()
                #desired_p1_loc_set.add((g, fp.name, 0, AMP_F, f.contig, '+'))
                #desired_p2_loc_set.add((g, fp.name, 0, AMP_R, f.contig, '-'))
                ampF_list = []
                ampR_list = []
                for pi, (seq, p) in enumerate(Primer.sequences.items()):
                    if ((p.checks[0] != None) and p.summarize(p.checks)):
                        
                        if (args.donor_specificity == 'any'):
                            
                            if any(((loc[0] == g) and (loc[3] == AMP_F)) for loc in p.locations):
                                ampF_list.append(p)
                            if any(((loc[0] == g) and (loc[3] == AMP_R)) for loc in p.locations):
                                ampR_list.append(p)
                            
                            #p_loc_set = set(loc[:-2] for loc in p.locations)
                            #if (desired_p1_loc_set.intersection(p_loc_set) == desired_p1_loc_set):
                            #    ampF_list.append(p)
                            #if (desired_p2_loc_set.intersection(p_loc_set) == desired_p2_loc_set):
                            #    ampR_list.append(p)
                            ##if any((loc[3] == AMP_F) for loc in p.locations): # allow for any feature/locus
                            ##    ampF_list.append(p)
                            ##if any((loc[3] == AMP_R) for loc in p.locations): # allow for any feature/locus
                            ##    ampR_list.append(p)
                            
                        elif (args.donor_specificity == 'all'):
                            p1_features = set()
                            p2_features = set()
                            for loc in p.locations:
                                if (loc[0] == g):
                                    if (loc[3] == AMP_F):
                                        p1_features.add(loc[1])
                                    if (loc[3] == AMP_R):
                                        p2_features.add(loc[1])
                            if (len(p1_features) == g_count):
                                ampF_list.append(p)
                            if (len(p2_features) == g_count):
                                ampR_list.append(p)
                        
                        elif (args.donor_specificity == 'exclusive'):
                            p1_features = set()
                            p2_features = set()
                            for loc in p.locations:
                                if (loc[0] == g):
                                    if (loc[3] == AMP_F):
                                        p1_features.add(loc[1])
                                    if (loc[3] == AMP_R):
                                        p2_features.add(loc[1])
                            if ((len(p1_features) > 0) and (len(p2_features) == 0)):
                                ampF_list.append(p)
                            if ((len(p2_features) > 0) and (len(p1_features) == 0)):
                                ampR_list.append(p)
                
                
                
                ampF_list = sorted(ampF_list, reverse=True)
                ampR_list = sorted(ampR_list, reverse=True)
                cls.logger.info('  len(ampF_list) = {}'.format(len(ampF_list)))
                cls.logger.info('  len(ampR_list) = {}'.format(len(ampR_list)))
                
            #    for pF in ampF_list:
            #        us_gene_p_list_dict[g][pF.sequence].append(pF)
            #    
            #    for pR in ampR_list:
            #        ds_gene_p_list_dict[g][pR.sequence].append(pR)
                
                
                # Allele-specific 'PrimerPair' calculations
                #if (args.primer_specificity == 'any'):
                #if args.allele_specific_primers:
                # >>> INDENT >>>
                cls.logger.info('Calculating: primer pairs...')
                check_count = 0
                # For each pair, we reset the timer
                start_time = time.time()
                time_expired = False
                
                uf_dr_paired_primers = []
                for i1, p1 in enumerate(ampF_list[:pp_subset_size]):
                    if ((time.time()-start_time) > args.primer_pair_limit):
                        time_expired = True
                    
                    if not time_expired:
                        for i2, p2 in enumerate(ampR_list[:pp_subset_size]):
                            # If it doesn't exist, add PrimerPair to database.
                            # Otherwise, do nothing.
                            PrimerPair(p1, p2)
                            
                            # Get the 'PrimerPair' object from the dict
                            pair = (p1.sequence, p2.sequence)
                            pp = PrimerPair.pairs.get(pair)
                            
                            if pp:
                                # Run checks
                                pp.progressive_check(cutoffs)
                                check_count += 1
                                
                                # If the 'PrimerPair' passes the checks, then add it
                                if ((pp.checks[0] != None) and Primer.summarize(pp.checks)):
                                    pp.weight = pp.get_weight(minimize=gg2feature_size)
                                    
                                    uf_dr_paired_primers.append(pp)
                cls.logger.info('  checked {} primer pairs'.format(check_count))
                
                uf_dr_paired_primers = sorted(uf_dr_paired_primers, reverse=True)
                
                cls.logger.info('  len(uf_dr_paired_primers) = {}'.format(len(uf_dr_paired_primers)))
                
                
                
                
                
                
                
                cls.logger.info('\t'.join(['ReversionDonor', 'feature', 'weight', 'PrimerPair']))
            
                for feature_name, f in Feature.features.items():
                    
                    fp = f.get_parent() # short for 'feature_parent'
                    
                    if (feature2gene[fp.name] == g):
                        
                        ###### alignment ######
                        pp_labels_list = []
                        ###### alignment ######/
                        
                        ### Copied from earlier in function ###
                        my_contig = contigs[f.contig]
                        
                        us_extend_seq = my_contig[f.start:fp.start]
                        ds_extend_seq = my_contig[fp.end:f.end]
                        mid_feature_seq = my_contig[fp.start:fp.end]
                        if isinstance(args.ki_dDNA, str):
                            # If a file is specified, then it has the knock-in DNA that should be stitched to flanking homology arms
                            #o_feature_name = None
                            #for f2g_f, f2g_g in feature2gene.items():
                            #    if ((fp.name == f2g_f) or (fp.name == f2g_g)):
                            #        o_feature_name = f2g_g
                            #        break
                            
                            ki_contigs = utils.old_load_fasta_file(args.ki_dDNA)
                            
                            try:
                                if fp.name in ki_contigs:
                                    mid_feature_seq = ki_contigs[fp.name]
                                else:
                                    mid_feature_seq = ki_contigs[feature2gene[fp.name]]
                            except KeyError:
                                raise Exception("The ki-dDNA file '{}' has no sequence with a primary header that matches '{}' or '{}'".format(args.ki_dDNA, fp.name, feature2gene[fp.name]))
                        else:
                            # Otherwise, generate KI dDNA that are wild type
                            pass
                        
                        orientation = '+'
                        
                        # Won't work if a region spans the start/end of a circular plasmid:
                        #   RRRFFFFFF------------RRR
                        # Also won't work if the Feature is too close to the beginning of the contig
                        region_F_start, region_F_stop = f.start-max(args.revert_homology_length), f.start-min(args.revert_homology_length)
                        region_F_start, region_F_stop = max(0, region_F_start), max(0, region_F_stop)
                        region_R_start, region_R_stop = f.end+min(args.revert_homology_length), f.end+max(args.revert_homology_length)
                        
                        region_F = my_contig[region_F_start:region_F_stop]
                        region_R = my_contig[region_R_start:region_R_stop]
                        ### Copied from earlier in function ###
                        
                        
                        #for pp_count, pp in enumerate(uf_dr_paired_primers[:pp_subset_size]):
                        for pp in uf_dr_paired_primers:
                            
                            # Pick the right-most p1, and the left-most p2 (corresponding to the smallest amplicon)
                            p1_se_list = []
                            for loc in pp.forward_primer.locations:
                                if ((loc[0] == g) and (loc[3] == AMP_F) and (loc[4] == f.contig)):
                                    p1_se_list.append([loc[6], loc[7], loc])
                            p1_location = max(p1_se_list)[2]
                            
                            p2_se_list = []
                            for loc in pp.reverse_primer.locations:
                                if ((loc[0] == g) and (loc[3] == AMP_R) and (loc[4] == f.contig)):
                                    p2_se_list.append([loc[6], loc[7], loc])
                            p2_location = min(p2_se_list)[2]
                            
                            start1, end1 = p1_location[6], f.start
                            start2, end2 = f.end, p2_location[7]
                            
                            #start1, end1 = region_F_start + pp.forward_primer.position, f.start
                            #start2, end2 = f.end, region_R_start+pp.reverse_primer.position+len(pp.reverse_primer.sequence)
                            upstream_seq = my_contig[start1:end1]
                            downstream_seq = my_contig[start2:end2]
                            
                            #dDNA = upstream_seq + feature_sequence + downstream_seq # This works, but doesn't take foreign knock-in DNA into account
                            dDNA = upstream_seq + us_extend_seq + mid_feature_seq + ds_extend_seq + downstream_seq
                            
                            new_obj = cls(feature_name, f.contig, orientation, dDNA, (start1, end2), spacer=None)
                            
                            rd = cls.sequences[dDNA]
                            
                            # Also need to list their compatible exDonors, and by extension, their exTargets
                            cls.logger.info('\t'.join([rd.name, feature_name, str(pp.get_joint_weight()), str(pp)])) # <<<<< In pp, the 'amplicon_size' attribute is calculated INCORRECTLY for some derived features. Also, the x-shift seems off in the alignment.
                            #cls.logger.info('  (new_obj == rd): {}'.format(new_obj == rd))
                            #cls.logger.info('  contig: {}, start1: {}, end2: {}'.format(f.contig, start1, end2))
                            #cls.logger.info('  F: {}'.format(sorted(pp.forward_primer.locations)))
                            #cls.logger.info('  R: {}'.format(sorted(pp.reverse_primer.locations)))
                            
                            ###### alignment ######
                            pp_labels_list.append('{: 18.15f}'.format(pp.get_joint_weight()) + ' ' + rd.name)  # <<<<< In pp, the 'amplicon_size' attribute is calculated INCORRECTLY for some derived features. Also, the x-shift seems off in the alignment
                            ###### alignment ######
                            
                            # ############# Write to Genbank flat file #############
                            # colors = {
                            #     'grey': '#DDDDDD',
                            #     'gray': '#DDDDDD',
                            #     'pink': '#FF9CCD',
                            #     'bpink': '#FF9C9A', # salmon
                            #     'dpink': '#DDB4DD', # grey-pink
                            # }
                            # segment_start, segment_end = max(0, f.start-2000), min(len(my_contig), f.end+2000)
                            # genome_segment = my_contig[segment_start:segment_end]
                            # gb = utils.GenBankFile(f.contig + ':' + str(segment_start+1) + '-' + str(segment_end), genome_segment) # chr:start-end
                            # #gb.add_annotation('source', 0, len(genome_segment), '+', organism=args.fasta[0], mol_type='genomic DNA') # only uses the first listed args.fasta
                            # gb.add_annotation('feature', f.start-segment_start, f.end-segment_start, f.strand, label='feature', ApEinfo_revcolor=colors['grey'], ApEinfo_fwdcolor=colors['grey'])
                            # gb.add_annotation('region', region_F_start-segment_start, region_F_stop-segment_start, '+', label='upstream_primer_region', ApEinfo_revcolor=colors['dpink'], ApEinfo_fwdcolor=colors['dpink'])
                            # gb.add_annotation('region', region_R_start-segment_start, region_R_stop-segment_start, '-', label='downstream_primer_region', ApEinfo_revcolor=colors['dpink'], ApEinfo_fwdcolor=colors['dpink'])
                            # gb.add_annotation('amplicon', start1-segment_start, end2-segment_start, '+', label='ampf_ampr_amplicon', ApEinfo_revcolor=colors['bpink'], ApEinfo_fwdcolor=colors['bpink']) # salmon
                            # gb.add_annotation('primer', start1 - segment_start, start1 - segment_start + len(pp.forward_primer.sequence), '+', label='ampf_primer', ApEinfo_revcolor=colors['pink'], ApEinfo_fwdcolor=colors['pink']) # pink
                            # gb.add_annotation('primer', end2 - segment_start-len(pp.reverse_primer.sequence), end2 - segment_start, '-', label='ampr_primer', ApEinfo_revcolor=colors['pink'], ApEinfo_fwdcolor=colors['pink']) # pink
                            # gb.write(os.path.join(args.folder, 'reversion-gb', 'pp-'+str(pp_count)+'.gb')) # pp_labels_list[-1]
                            # ############# Write to Genbank flat file #############
                    
                
                        ###### alignment ######
                        label_list = ['us_region', 'us_skipped', 'us_extend_feature', 'feature', 'ds_extend_feature', 'ds_skipped', 'ds_region']
                        sequence_list = [region_F, my_contig[region_F_stop:f.start], us_extend_seq, mid_feature_seq, ds_extend_seq, my_contig[f.end:region_R_start], region_R]
                        aln_out = nucleotides.make_labeled_primer_alignments(label_list, sequence_list, '{:>18}'.format('weight')+' '+f.contig, pp_labels_list, uf_dr_paired_primers[:pp_subset_size], left_pos=region_F_start)
                        
                        cls.logger.info(feature_name)
                        for oline in aln_out:
                            cls.logger.info(oline)
                        ###### alignment ######
        
        else: # This is when (args.revert_amplification_primers == False)
            # There should be at least one ReversionDonor for each feature
            for feature_name, f in Feature.features.items():
                cls.logger.info("Scanning feature '{}:{}:{}:{}..{}' for amplification primers".format(feature_name, f.contig, f.strand, f.start, f.end))
                
                my_contig = contigs[f.contig]
                
                #feature_length = f.end - f.start
                #feature_sequence = my_contig[f.start:f.end]
                
                fp = f.get_parent() # short for 'feature_parent'
                us_extend_seq = my_contig[f.start:fp.start]
                ds_extend_seq = my_contig[fp.end:f.end]
                mid_feature_seq = my_contig[fp.start:fp.end]
                if isinstance(args.ki_dDNA, str):
                    # If a file is specified, then it has the knock-in DNA that should be stitched to flanking homology arms
                    #o_feature_name = None
                    #for f2g_f, f2g_g in feature2gene.items():
                    #    if ((fp.name == f2g_f) or (fp.name == f2g_g)):
                    #        o_feature_name = f2g_g
                    #        break
                    
                    ki_contigs = utils.old_load_fasta_file(args.ki_dDNA)
                    
                    try:
                        if fp.name in ki_contigs:
                            mid_feature_seq = ki_contigs[fp.name]
                        else:
                            mid_feature_seq = ki_contigs[feature2gene[fp.name]]
                    except KeyError:
                        raise Exception("The ki-dDNA file '{}' has no sequence with a primary header that matches '{}' or '{}'".format(args.ki_dDNA, fp.name, feature2gene[fp.name]))
                else:
                    # Otherwise, generate KI dDNA that are wild type
                    pass
                
                region_F_start, region_F_stop = f.start-max(args.revert_homology_length), f.start-min(args.revert_homology_length)
                region_F_start, region_F_stop = max(0, region_F_start), max(0, region_F_stop)
                region_R_start, region_R_stop = f.end+min(args.revert_homology_length), f.end+max(args.revert_homology_length)
                
                start1, end1 = region_F_start, f.start
                start2, end2 = f.end, region_R_stop
                upstream_seq = my_contig[start1:end1]
                downstream_seq = my_contig[start2:end2]
                
                #dDNA = upstream_seq + feature_sequence + downstream_seq # This works, but doesn't take foreign knock-in DNA into account
                dDNA = upstream_seq + us_extend_seq + mid_feature_seq + ds_extend_seq + downstream_seq
                
                orientation = '+'
                
                cls(feature_name, f.contig, orientation, dDNA, (start1, end2), spacer=None)
            
        ### End generate_donors() method ###
    
    @classmethod
    def generate_donors_old(cls, args, contigs, feature2gene):
        """
        Creates the DNA oligo with the structure:
        [amp-F primer]                                                               [amp-R primer]
        [upstream homology][us expanded feature][feature][ds expanded feature][downstream homology]
        """
        
        cls.logger.info("Generating 'ReversionDonor' objects")
        
        if (args.revert_amplification_primers):
            tm_max_difference = 4.0
            amplicon_size = (2*min(args.revert_homology_length), 2*max(args.revert_homology_length))
            tm_range = (52, 64)
            primer_length_range = (19, 32)
            min_delta_g = -5.0
            
            #subset_size = 1000 # 500 # Temporarily limit number of Primer objects used as input to the pair() function
            pp_subset_size = 1000 # Temporarily limit the number of PrimerPairs that are used to create ReversionDonor objects
            temp_folder = os.path.join('/dev/shm/addtag', os.path.basename(args.folder))
            
            # Make a folder to store the *.gb Genbank flat files
            os.makedirs(os.path.join(args.folder, 'reversion-gb'), exist_ok=True)
            
            # Make list of genes, and the expected input feature numbers
            gene_dict = defaultdict(int)
            us_gene_p_list_dict = {}
            ds_gene_p_list_dict = {}
            for feature_name, f in Feature.features.items():
                fp = f.get_parent() # short for 'feature_parent'
                g = feature2gene[fp.name]
                gene_dict[g] += 1
                
                # Make defaultdicts to count homologous Primer objects
                us_gene_p_list_dict[g] = defaultdict(list)
                ds_gene_p_list_dict[g] = defaultdict(list)
            
            # There should be at least one ReversionDonor for each feature
            for feature_name, f in Feature.features.items():
                cls.logger.info("Scanning feature '{}:{}:{}:{}..{}' for amplification primers".format(feature_name, f.contig, f.strand, f.start, f.end))
                
                my_contig = contigs[f.contig]
                
                #feature_length = f.end - f.start
                #feature_sequence = my_contig[f.start:f.end]
                
                fp = f.get_parent() # short for 'feature_parent'
                g = feature2gene[fp.name]
                
                us_extend_seq = my_contig[f.start:fp.start]
                ds_extend_seq = my_contig[fp.end:f.end]
                mid_feature_seq = my_contig[fp.start:fp.end]
                if isinstance(args.ki_dDNA, str):
                    # If a file is specified, then it has the knock-in DNA that should be stitched to flanking homology arms
                    #o_feature_name = None
                    #for f2g_f, f2g_g in feature2gene.items():
                    #    if ((fp.name == f2g_f) or (fp.name == f2g_g)):
                    #        o_feature_name = f2g_g
                    #        break
                    
                    ki_contigs = utils.old_load_fasta_file(args.ki_dDNA)
                    
                    try:
                        if fp.name in ki_contigs:
                            mid_feature_seq = ki_contigs[fp.name]
                        else:
                            mid_feature_seq = ki_contigs[feature2gene[fp.name]]
                    except KeyError:
                        raise Exception("The ki-dDNA file '{}' has no sequence with a primary header that matches '{}' or '{}'".format(args.ki_dDNA, fp.name, feature2gene[fp.name]))
                else:
                    # Otherwise, generate KI dDNA that are wild type
                    pass
                
                orientation = '+'
                
                # Won't work if a region spans the start/end of a circular plasmid:
                #   RRRFFFFFF------------RRR
                # Also won't work if the Feature is too close to the beginning of the contig
                region_F_start, region_F_stop = f.start-max(args.revert_homology_length), f.start-min(args.revert_homology_length)
                region_F_start, region_F_stop = max(0, region_F_start), max(0, region_F_stop)
                region_R_start, region_R_stop = f.end+min(args.revert_homology_length), f.end+max(args.revert_homology_length)
                
                region_F = my_contig[region_F_start:region_F_stop]
                region_R = my_contig[region_R_start:region_R_stop]
                
                cls.logger.info('Calculating upstream_F primers...')
                upstream_F = sorted(
                    args.selected_oligo.scan(region_F, 'left',  primer_size=primer_length_range, tm_range=tm_range, min_delta_g=min_delta_g, folder=temp_folder, time_limit=args.primer_scan_limit),
                    key=lambda x: x.weight,
                    reverse=True
                )
                cls.logger.info('  len(upstream_F) = {}'.format(len(upstream_F)))
                #upstream_F = upstream_F[:subset_size]
                #if ((subset_size != None) and (len(upstream_F) > subset_size)):
                #    cls.logger.info('upstream_F: skipping {}/{} calculated primers'.format(max(0, len(upstream_F)-subset_size), len(upstream_F)))
                for pF in upstream_F:
                    us_gene_p_list_dict[g][pF.sequence].append(pF)
                
                cls.logger.info('Calculating downstream_R primers...')
                downstream_R = sorted(
                    args.selected_oligo.scan(region_R, 'right', primer_size=primer_length_range, tm_range=tm_range, min_delta_g=min_delta_g, folder=temp_folder, time_limit=args.primer_scan_limit),
                    key=lambda x: x.weight,
                    reverse=True
                )
                cls.logger.info('  len(downstream_R) = {}'.format(len(downstream_R)))
                #downstream_R = downstream_R[:subset_size]
                #if ((subset_size != None) and (len(downstream_R) > subset_size)):
                #    cls.logger.info('downstream_R: skipping {}/{} calculated primers'.format(max(0, len(downstream_R)-subset_size), len(downstream_R)))
                for pR in downstream_R:
                    ds_gene_p_list_dict[g][pR.sequence].append(pR)
                
                
                
                
                if args.allele_specific_primers:
                    
                    cls.logger.info('Calculating: primer pairs...')
                    uf_dr_paired_primers = sorted(
                        args.selected_oligo.pair(upstream_F, downstream_R, amplicon_size_range=amplicon_size, tm_max_difference=tm_max_difference, intervening=(f.start-region_F_stop)+(region_R_start-f.end), min_delta_g=min_delta_g, folder=temp_folder, time_limit=args.primer_pair_limit),
                        key=lambda x: x.get_joint_weight(),
                        reverse=True
                    )
                    cls.logger.info('  len(uf_dr_paired_primers) = {}'.format(len(uf_dr_paired_primers)))
                    
                    # Ammend the intervening sequence
                    for pp in uf_dr_paired_primers:
                        pp.intervening = region_R_start - region_F_stop
                    
                    #if ((pp_subset_size != None) and (len(uf_dr_paired_primers) > pp_subset_size)):
                    #    cls.logger.info('Skipping {}/{} primer pairs'.format(max(0, len(uf_dr_paired_primers)-pp_subset_size), len(uf_dr_paired_primers)))
                    
                    
                    ###### alignment ######
                    pp_labels_list = []
                    ###### alignment ######
                    
                    
                    cls.logger.info('\t'.join(['ReversionDonor', 'feature', 'weight', 'PrimerPair']))
                    #for pp_count, pp in enumerate(uf_dr_paired_primers[:pp_subset_size]):
                    for pp_count, pp in enumerate(uf_dr_paired_primers):
                        start1, end1 = region_F_start + pp.forward_primer.position, f.start
                        start2, end2 = f.end, region_R_start+pp.reverse_primer.position+len(pp.reverse_primer.sequence)
                        upstream_seq = my_contig[start1:end1]
                        downstream_seq = my_contig[start2:end2]
                        
                        #dDNA = upstream_seq + feature_sequence + downstream_seq # This works, but doesn't take foreign knock-in DNA into account
                        dDNA = upstream_seq + us_extend_seq + mid_feature_seq + ds_extend_seq + downstream_seq
                        
                        cls(feature_name, f.contig, orientation, dDNA, (start1, end2), spacer=None)
                        
                        rd = cls.sequences[dDNA]
                        
                        # Also need to list their compatible exDonors, and by extension, their exTargets
                        cls.logger.info('\t'.join([rd.name, feature_name, str(pp.get_joint_weight()), str(pp)])) # <<<<< In pp, the 'amplicon_size' attribute is calculated INCORRECTLY for some derived features. Also, the x-shift seems off in the alignment.
                        
                        ###### alignment ######
                        pp_labels_list.append('{: 18.15f}'.format(pp.get_joint_weight()) + ' ' + rd.name)  # <<<<< In pp, the 'amplicon_size' attribute is calculated INCORRECTLY for some derived features. Also, the x-shift seems off in the alignment
                        ###### alignment ######
                        
                        # ############# Write to Genbank flat file #############
                        # colors = {
                        #     'grey': '#DDDDDD',
                        #     'gray': '#DDDDDD',
                        #     'pink': '#FF9CCD',
                        #     'bpink': '#FF9C9A', # salmon
                        #     'dpink': '#DDB4DD', # grey-pink
                        # }
                        # segment_start, segment_end = max(0, f.start-2000), min(len(my_contig), f.end+2000)
                        # genome_segment = my_contig[segment_start:segment_end]
                        # gb = utils.GenBankFile(f.contig + ':' + str(segment_start+1) + '-' + str(segment_end), genome_segment) # chr:start-end
                        # #gb.add_annotation('source', 0, len(genome_segment), '+', organism=args.fasta[0], mol_type='genomic DNA') # only uses the first listed args.fasta
                        # gb.add_annotation('feature', f.start-segment_start, f.end-segment_start, f.strand, label='feature', ApEinfo_revcolor=colors['grey'], ApEinfo_fwdcolor=colors['grey'])
                        # gb.add_annotation('region', region_F_start-segment_start, region_F_stop-segment_start, '+', label='upstream_primer_region', ApEinfo_revcolor=colors['dpink'], ApEinfo_fwdcolor=colors['dpink'])
                        # gb.add_annotation('region', region_R_start-segment_start, region_R_stop-segment_start, '-', label='downstream_primer_region', ApEinfo_revcolor=colors['dpink'], ApEinfo_fwdcolor=colors['dpink'])
                        # gb.add_annotation('amplicon', start1-segment_start, end2-segment_start, '+', label='ampf_ampr_amplicon', ApEinfo_revcolor=colors['bpink'], ApEinfo_fwdcolor=colors['bpink']) # salmon
                        # gb.add_annotation('primer', start1 - segment_start, start1 - segment_start + len(pp.forward_primer.sequence), '+', label='ampf_primer', ApEinfo_revcolor=colors['pink'], ApEinfo_fwdcolor=colors['pink']) # pink
                        # gb.add_annotation('primer', end2 - segment_start-len(pp.reverse_primer.sequence), end2 - segment_start, '-', label='ampr_primer', ApEinfo_revcolor=colors['pink'], ApEinfo_fwdcolor=colors['pink']) # pink
                        # gb.write(os.path.join(args.folder, 'reversion-gb', 'pp-'+str(pp_count)+'.gb')) # pp_labels_list[-1]
                        # ############# Write to Genbank flat file #############
                        
                    
                    ###### alignment ######
                    label_list = ['us_region', 'us_skipped', 'us_extend_feature', 'feature', 'ds_extend_feature', 'ds_skipped', 'ds_region']
                    sequence_list = [region_F, my_contig[region_F_stop:f.start], us_extend_seq, mid_feature_seq, ds_extend_seq, my_contig[f.end:region_R_start], region_R]
                    aln_out = nucleotides.make_labeled_primer_alignments(label_list, sequence_list, '{:>18}'.format('weight')+' '+f.contig, pp_labels_list, uf_dr_paired_primers[:pp_subset_size])
                    
                    cls.logger.info(feature_name)
                    for oline in aln_out:
                        cls.logger.info(oline)
                    ###### alignment ######
            
            
            ##### BEGIN (non-allele-specific processing) #####
            
            # Go through each Feature (requires homology table)
            # And keep only primers present in all Features of the same gene
            us_gene_primers = {}
            ds_gene_primers = {}
            gene_count = {}
            for feature_name, f in Feature.features.items():
                fp = f.get_parent() # short for 'feature_parent'
                g = feature2gene[fp.name]
                if not g in gene_count:
                    gene_count[g] = 0
                
                upstream_F = []
                for pF_seq, pF_list in us_gene_p_list_dict[g].items():
                    if (len(pF_list) == gene_dict[g]):
                        upstream_F.append(pF_list[gene_count[g]])
                
                downstream_R = []
                for pF_seq, pR_list in ds_gene_p_list_dict[g].items():
                    if (len(pR_list) == gene_dict[g]):
                        downstream_R.append(pR_list[gene_count[g]])
                
                # Do primer pairing for each homolog group (gene)
                # And NOT for each feature
                
                
                
                cls.logger.info("Scanning feature '{}:{}:{}:{}..{}' for multi-allele amplification primers".format(feature_name, f.contig, f.strand, f.start, f.end))
                
                my_contig = contigs[f.contig]
                
                us_extend_seq = my_contig[f.start:fp.start]
                ds_extend_seq = my_contig[fp.end:f.end]
                mid_feature_seq = my_contig[fp.start:fp.end]
                if isinstance(args.ki_dDNA, str):
                    # If a file is specified, then it has the knock-in DNA that should be stitched to flanking homology arms
                    ki_contigs = utils.old_load_fasta_file(args.ki_dDNA)
                    
                    try:
                        if fp.name in ki_contigs:
                            mid_feature_seq = ki_contigs[fp.name]
                        else:
                            mid_feature_seq = ki_contigs[feature2gene[fp.name]]
                    except KeyError:
                        raise Exception("The ki-dDNA file '{}' has no sequence with a primary header that matches '{}' or '{}'".format(args.ki_dDNA, fp.name, feature2gene[fp.name]))
                else:
                    # Otherwise, generate KI dDNA that are wild type
                    pass
                
                orientation = '+'
                
                # Won't work if a region spans the start/end of a circular plasmid:
                #   RRRFFFFFF------------RRR
                # Also won't work if the Feature is too close to the beginning of the contig
                region_F_start, region_F_stop = f.start-max(args.revert_homology_length), f.start-min(args.revert_homology_length)
                region_F_start, region_F_stop = max(0, region_F_start), max(0, region_F_stop)
                region_R_start, region_R_stop = f.end+min(args.revert_homology_length), f.end+max(args.revert_homology_length)
                
                
                
                                    
                cls.logger.info('Calculating: primer pairs...')
                uf_dr_paired_primers = sorted(
                    args.selected_oligo.pair(upstream_F, downstream_R, amplicon_size_range=amplicon_size, tm_max_difference=tm_max_difference, intervening=(f.start-region_F_stop)+(region_R_start-f.end), min_delta_g=min_delta_g, folder=temp_folder, time_limit=args.primer_pair_limit),
                    key=lambda x: x.get_joint_weight(),
                    reverse=True
                )
                cls.logger.info('  len(uf_dr_paired_primers) = {}'.format(len(uf_dr_paired_primers)))
                
                # Ammend the intervening sequence
                for pp in uf_dr_paired_primers:
                    pp.intervening = region_R_start - region_F_stop
                
                ###### alignment ######
                pp_labels_list = []
                ###### alignment ######
                
                cls.logger.info('\t'.join(['ReversionDonor', 'feature', 'weight', 'PrimerPair']))
                for pp_count, pp in enumerate(uf_dr_paired_primers):
                    start1, end1 = region_F_start + pp.forward_primer.position, f.start
                    start2, end2 = f.end, region_R_start+pp.reverse_primer.position+len(pp.reverse_primer.sequence)
                    upstream_seq = my_contig[start1:end1]
                    downstream_seq = my_contig[start2:end2]
                    
                    #dDNA = upstream_seq + feature_sequence + downstream_seq # This works, but doesn't take foreign knock-in DNA into account
                    dDNA = upstream_seq + us_extend_seq + mid_feature_seq + ds_extend_seq + downstream_seq
                    
                    cls(feature_name, f.contig, orientation, dDNA, (start1, end2), spacer=None)
                    
                    rd = cls.sequences[dDNA]
                    
                    # Also need to list their compatible exDonors, and by extension, their exTargets
                    cls.logger.info('\t'.join([rd.name, feature_name, str(pp.get_joint_weight()), str(pp)]))  # <<<<< In pp, the 'amplicon_size' attribute is calculated INCORRECTLY for some derived features. Also, the x-shift seems off in the alignment
                    
                    ###### alignment ######
                    pp_labels_list.append('{: 18.15f}'.format(pp.get_joint_weight()) + ' ' + rd.name)  # <<<<< In pp, the 'amplicon_size' attribute is calculated INCORRECTLY for some derived features. Also, the x-shift seems off in the alignment
                    ###### alignment ######
                
                ###### alignment ######
                label_list = ['us_region', 'us_skipped', 'us_extend_feature', 'feature', 'ds_extend_feature', 'ds_skipped', 'ds_region']
                sequence_list = [region_F, my_contig[region_F_stop:f.start], us_extend_seq, mid_feature_seq, ds_extend_seq, my_contig[f.end:region_R_start], region_R]
                aln_out = nucleotides.make_labeled_primer_alignments(label_list, sequence_list, '{:>18}'.format('weight')+' '+f.contig, pp_labels_list, uf_dr_paired_primers[:pp_subset_size])
                
                cls.logger.info(feature_name)
                for oline in aln_out:
                    cls.logger.info(oline)
                ###### alignment ######
                
                
                
                gene_count[g] += 1
                
            
            ##### END (non-allele-specific processing) #####
            
        else: # This is when (args.revert_amplification_primers == False)
            # There should be at least one ReversionDonor for each feature
            for feature_name, f in Feature.features.items():
                cls.logger.info("Scanning feature '{}:{}:{}:{}..{}' for amplification primers".format(feature_name, f.contig, f.strand, f.start, f.end))
                
                my_contig = contigs[f.contig]
                
                #feature_length = f.end - f.start
                #feature_sequence = my_contig[f.start:f.end]
                
                fp = f.get_parent() # short for 'feature_parent'
                us_extend_seq = my_contig[f.start:fp.start]
                ds_extend_seq = my_contig[fp.end:f.end]
                mid_feature_seq = my_contig[fp.start:fp.end]
                if isinstance(args.ki_dDNA, str):
                    # If a file is specified, then it has the knock-in DNA that should be stitched to flanking homology arms
                    #o_feature_name = None
                    #for f2g_f, f2g_g in feature2gene.items():
                    #    if ((fp.name == f2g_f) or (fp.name == f2g_g)):
                    #        o_feature_name = f2g_g
                    #        break
                    
                    ki_contigs = utils.old_load_fasta_file(args.ki_dDNA)
                    
                    try:
                        if fp.name in ki_contigs:
                            mid_feature_seq = ki_contigs[fp.name]
                        else:
                            mid_feature_seq = ki_contigs[feature2gene[fp.name]]
                    except KeyError:
                        raise Exception("The ki-dDNA file '{}' has no sequence with a primary header that matches '{}' or '{}'".format(args.ki_dDNA, fp.name, feature2gene[fp.name]))
                else:
                    # Otherwise, generate KI dDNA that are wild type
                    pass
                
                region_F_start, region_F_stop = f.start-max(args.revert_homology_length), f.start-min(args.revert_homology_length)
                region_F_start, region_F_stop = max(0, region_F_start), max(0, region_F_stop)
                region_R_start, region_R_stop = f.end+min(args.revert_homology_length), f.end+max(args.revert_homology_length)
                
                start1, end1 = region_F_start, f.start
                start2, end2 = f.end, region_R_stop
                upstream_seq = my_contig[start1:end1]
                downstream_seq = my_contig[start2:end2]
                
                #dDNA = upstream_seq + feature_sequence + downstream_seq # This works, but doesn't take foreign knock-in DNA into account
                dDNA = upstream_seq + us_extend_seq + mid_feature_seq + ds_extend_seq + downstream_seq
                
                orientation = '+'
                
                cls(feature_name, f.contig, orientation, dDNA, (start1, end2), spacer=None)
            
        ### End generate_donors() method ###
    
    @classmethod
    def generate_donors_old_old(cls, args, contigs):
        """
        Creates the DNA oligo with the structure:
        [upstream homology][original feature][downstream homology]
        """
        # There should be one ReversionDonor for each feature
        for feature_name, f in Feature.features.items():
        #for feature in features:
            # First, get the homology blocks up- and down-stream of the feature
            #contig, start, end, strand = features[feature]
            
            feature_length = f.end - f.start
            
            my_contig = contigs[f.contig]
            orientation = '+'
            
            for us_hom in range(args.revert_upstream_homology[0], args.revert_upstream_homology[1]+1):
                for ds_hom in range(args.revert_downstream_homology[0], args.revert_downstream_homology[1]+1):
                    if (args.revert_donor_lengths[0] <= us_hom+feature_length+ds_hom <= args.revert_donor_lengths[1]):
                        # to implement:
                        #   make revert_upstream_homology and revert_downstream_homology exclude the trim sequences
                        if (f.strand == '+'):
                            start1, end1 = f.start-us_hom, f.start
                            start2, end2 = f.end, f.end+ds_hom
                        else:
                            start1, end1 = f.start-ds_hom, f.start
                            start2, end2 = f.end, f.end+us_hom
                        upstream = my_contig[start1:end1]
                        downstream = my_contig[start2:end2]
                        feature_sequence = my_contig[f.start:f.end]
                        
                        dDNA = upstream + feature_sequence + downstream
                        #for re_seq, re_target in ReversionTarget.indices.items():
                        #    if feature in [x[0] for x in re_target.locations]:
                        #        cls(feature, contig, orientation, (start1, end1), '..'.join(map(str, [start, end])), (start2, end2), dDNA, re_target)
                        #cls(feature, contig, orientation, (start1, end1), '..'.join(map(str, [start, end])), (start2, end2), dDNA, None)
                        cls(feature_name, f.contig, orientation, dDNA, (start1, end2), spacer=None)
