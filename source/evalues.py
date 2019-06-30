#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/evalues.py

# List general Python imports
import sys
import math

# Import included AddTag-specific modules
import nucleotides
import utils
#from algorithms.distance_matrix import load_scores

# Calculating the E-value requires three things:
#  1) Scoring matrix (match, mismatch, gap_open, gap_extension) scores
#     In terms of the scoring matrix, this simple version of the E-value
#     calculation only considers:
#       * match scores
#       * mismatch scores
#     And it does not take into account:
#       x gap_open_score
#       x gap_extension_score
#  2) Estimating database-dependent parameters
#       * Nucleotide composition of database
#       * Database length
#  3) Estimating alignment-dependent parameters
#       * Alignment score of query against subject
#       * Query length

# References:
#  * Karlin S, et al. Methods for assessing the statistical significance of
#    molecular sequence features by using general scoring schemes. Proceedings
#    of the National Academy of Sciences. 1990;87(6):2264-8.
#  * Korf I, et al. BLAST: O'Reilly Media, Inc.; 2003.
#  * Gertz EM. BLAST Scoring Parameters. National Center for Biotechnology
#    Information. 2005.

######### Duplicate function from source/algorithms/distance_matrix.py #########
def load_scores(file_path, sep="\t"):
    """Load the alignment scores"""
    import regex
    scoring = {}
    with open(file_path, 'r') as flo:
        for line in flo:
            line = line.rstrip()
            if (len(line) > 0):
                m = regex.match(r'^\s*#', line)
                if not m:
                    sline = line.split(sep)
                    if (sline[2] == ''):
                        scoring[(sline[0], sline[1])] = float(sline[3])
                    else:
                        scoring[(sline[0], sline[1], sline[2])] = float(sline[3])
    return scoring
################################################################################

def nats_to_bits(nats):
    '''Convert nats to bits'''
    return nats/math.log(2)

def bits_to_nats(bits):
    '''Convert bits to nats'''
    return bits*math.log(2)

##### These 3 fxns don't apply to sum statistics #####
def raw_score_to_bit_score(raw_score, k, l):
    # raw_score:
    # k:
    # l: lambda in nats
    return nats2bits(l*raw_score - math.log(k))

def raw_score_to_expect(raw_score, k, l, m, n):
    # raw_score:
    # k:
    # l: lambda in nats
    # m: effective length of the query
    # n: effective length of the database
    return k*m*n*math.exp(-l*raw_score)

def bit_score_to_expect(bit_score, m, n):
    '''
    Note that k doesn't appear in this equation; it has already been accounted
    for when deriving the normalized (nat) score (raw_score_to_bit_score).
    
    The formula is similar to that used to calculate an Expect from a nat score.
    However, the base of the exponent is 2 rather than e because you're using
    bits rather than nats.
    '''
    # bit_score:
    # m: effective length of the query
    # n: effective length of the database
    return m*n* 2**(-bit_score)

def expected_raw_score_if_e_value_is_1(k, m, n, l):
    """
    Calculate the raw alignment score if
    e-value = 1
    """
    # k: ???
    # m: effective length of query 
    # n: effective length of database 
    # l: lambda in nats
    return  math.log(k*m*n)/l;
######################################################

def expected_HSP_length(k, m, n, h):
    # k: ???
    # m: actual length of query 
    # n: actual length of database 
    # h: average nats/aligned pair
    return math.log(k*m*n)/h;

def effective_query_length(m, expected_HSP_length, k):
    # m: actual sequence length of query
    # expected_HSP_length: ???
    # k: gapped k

    m_prime = m - expected_HSP_length
    
    # Cap the minimum effective query length to 1/k.
    # Notice that no effective length of the query or the database can ever be
    # less than 1/k. Setting an effective length to 1/k basically amounts to
    # ignoring a short sequence for statistical purposes; in cases when both
    # m and n are less than 1/k, BLAST searches are ill-advised.
    if (m_prime < 1/k):
        return 1/k
    else:
        return m_prime

def effective_db_length(n, expected_HSP_length, num_seqs_in_db, k):
    # n: actual, summed length of all subjects (across all sequences)
    # k: gapped k

    n_prime = n - (num_seqs_in_db*expected_HSP_length)
    
    # Cap the minimum effective query length to 1/k.
    # Notice that no effective length of the query or the database can ever be
    # less than 1/k. Setting an effective length to 1/k basically amounts to
    # ignoring a short sequence for statistical purposes; in cases when both
    # m and n are less than 1/k, BLAST searches are ill-advised.
    if (n_prime < 1/k):
        return 1/k;
    else:
        return n_prime

def effective_db_length_from_contigs(contigs, k, expected_HSP_length):
    n_primer = sum(map(len, contigs.values())) - (len(contigs)*expected_HSP_length)
    # Cap the minimum effective query length to 1/k.
    # Notice that no effective length of the query or the database can ever be
    # less than 1/k. Setting an effective length to 1/k basically amounts to
    # ignoring a short sequence for statistical purposes; in cases when both
    # m and n are less than 1/k, BLAST searches are ill-advised.
    if (n_prime < 1/k):
        return 1/k;
    else:
        return n_prime

def search_space(m, n):
    return m*n

def e_value_to_p_value(e_value):
    '''Convert e-value to p-value'''
    return 1-math.exp(-e_value);

def p_value_to_e_value(p_value):
    '''convert p-value to e-value'''
    return -math.log(1-p_value)

def sumScore(raw_scores, k, l, m, n, g):
    # raw_scores: raw scores are in a list/tuple
    # k:
    # l: lambda in nats
    # m: effective length of query sequence
    # n: effective length of sbjct sequence
    # g: gap size; for NCBI-BLAST this value is 50.
    
    r = len(raw_scores)
    if (r == 1):
        raise Exception("do not take sum for a single score!")
    
    total_raw_score = 0
    for individual_raw_score in raw_scores:
        total_raw_score += individual_raw_score
    
    n_score = l*total_raw_score;
    return n_score - math.log(k*m*n)  - (r -1)*(math.log(k)+2*math.log(g)) - math.log(math.factorial(r))

def real_subject_length(contigs):
    s = 0
    for k, v in contigs.items():
        s += len(v)
    return s

def pairwiseSumP(sumS, r):
    # sumS: the sum score
    # r: number of HSPs being combined
    return (math.exp(-sumS) * sumS**math.factorial(r-1)) / (math.factorial(r)*math.factorial(r -1))

def rTestsCorrectedP (r, sum_p, beta):
    # beta: gap decay constant
    return sum_p / (beta**(r-1) * (1-beta))

def dbSizeCorrectedExpect(sumP, effective_length_db, sbjct_seq_length):
    # sumP: ???
    # effective_length_db: different than the book (???)
    # sbjct_seq_length: ???
    return (effective_length_db/sbjct_seq_length)*sumP # Different than the book

def drange(start, stop, step=1, rounding=None):
    """
    A generator for a range that can use decimal steps
    Only works for increasing numbers
    """
    # For more info in floating point error in Python, see:
    #  https://docs.python.org/3/tutorial/floatingpoint.html
    i = 0
    r = start
    while (r < stop):
        if (rounding == None):
            yield r
        else:
            yield round(r, rounding)
        i += 1
        r = start + i * step

def alignment_probabilities(query_compositions, db_compositions):
    # Assumes both compositions have identical set of dict keys
    #query_compositions = {'A':0.05, 'C':0.05, 'G':0.45, 'T':0.45}
    #db_compositions = {'A':0.375, 'C':0.125, 'G':0.125, 'T':0.375}
    
    aprob = {}
    
    for qk, qv in query_compositions.items():
        for dbk, dbv in db_compositions.items():
            aprob[(qk, dbk)] = qv * dbv
    
    return aprob

def sign(x):
    """
    Returns:
       1 if x >= 0
      -1 if x <  0
    """
    if (x < 0):
        return -1
    else:
        return 1

def is_positive(x):
    """
    Returns:
      True   if x >= 0
      False  if x <  0
    """
    return x >= 0

def calculate_k(matrix, my_lambda, H):
    '''
    The input includes Lambda, H, and an array of probabilities for each score.
    
    Based on the input scoring matrix, this function will calculate K
    according to one of 4 different equations.
    
    Returns K (float)
    '''
    # l: lambda in nats
    # h: relative entropy?
    # matrix: dictionary containing scoring matrix
    
    def score_probabilities(scores):
        '''
        Convert scoring matrix into a dict containing the probabilities of each score.
        keys are scores, values are probabilities
        returns:
          {1: 0.25, -1: 0.25, -1.2: 0.5}
        '''
        d = {}
        for k, v in scores.items():
            d.setdefault(v, []).append(k)
        p = {}
        n = len(scores)
        for k, v in d.items():
            p[k] = len(v)/n
        return p
    
    def full_score_probabilities(sprob, divisor, rounding=2):
        '''
        Input: sprob: a sparse dict (e.g. sprob={-2: 0.7, 0: 0.1, 3: 0.2})
               divisor: an int/float
               rounding: number of decimal places to round (should be same as divisor)
        '''
        
        fsprob = {}
        low = min(sprob)
        high = max(sprob)
        
        for d in drange(low, high+divisor, divisor, rounding=rounding):
            fsprob[d] = 0.0
        for k, v in sprob.items():
            fsprob[k] = v
        
        return fsprob
    
    def matrix_average(scores):
        s = 0
        for k, v in scores.items():
            s += v
        return s/len(scores)
    
    def gcd_float(a, b, decimals=2):
        """
        Finds the greatest common divisor of 'a' and 'b', which can
        be floating points. Decimals will be rounded to the 'decimals' place
        """
        # Convert to positive numbers, and round to decimal place
        raa = round(abs(a), decimals)
        rab = round(abs(b), decimals)
        
        # Multiply to turn into an integer
        iraa = int(raa * 10**decimals)
        irab = int(rab * 10**decimals)
        
        # Get common divisor
        igcd = math.gcd(iraa, irab)
        
        # Convert to float
        gcd = igcd / 10**decimals
        
        return gcd
    
    def simple_gcd(x, y):
        while y != 0:
            (x, y) = (y, x % y)
        return x
    
    def nonzero_e_to_the_x_minus_1(x):
        '''
        Returns (e**x) - 1, but prevents 0 from being returned
        Will always return a non-zero number
        '''
        r = math.exp(x)-1
        if (r <= 0):
            return x
        else:
            return r
    
    def score_probabilities_of_length(sp, previous_permuted_sp, low, high, divisor, length, rounding=2):
        '''
        This is an iterative function that permutes through every possible
        score combination of 'length'. It will return a dict with the probability of
        any given score. Scores are keys, probabilities are values.
        
        For example:
          sp = {-3: 0.6, -1: 0.1, 0: 0.1, 2: 0.2}
          spj = sp
          spj = score_probabilities_of_length(sp, spj, min(sp), max(sp), 1, 1)
                {-3: 0.6, -2: 0.0, -1: 0.1, 0: 0.1, 1: 0.0, 2: 0.2}
          spj = score_probabilities_of_length(sp, spj, min(sp), max(sp), 1, 2)
                {-6: 0.36, -5: 0.0, -4: 0.12, -3: 0.12, -2: 0.01, -1: 0.26, 0: 0.01, 1: 0.04, 2: 0.04, 3: 0.0, 4: 0.04}
        '''
        # Populate as empty
        permuted_sp = {}
        for d in drange(low*length, high*length+divisor, divisor, rounding=2):
            permuted_sp[d] = 0.0
        
        if (length == 1):
            for s, p in sp.items():
                permuted_sp[s] = p
        
        elif (length > 1):
            # Add the probabilities
            for prev_s, prev_p in previous_permuted_sp.items():
                for s, p in sp.items():
                    new_score = round(s + prev_s, rounding)
                    new_probability = p * prev_p
                    permuted_sp[new_score] += new_probability
        
        return permuted_sp
    
    def calculate_sigma_bar(sp, my_lambda, low, high, divisor, iteration_limit=1000, sum_limit=0.001):
        '''
        Calculates sigma_bar
        '''
        spj = sp
        sigma_inner = 1
        sigma_bar = 0.0
        j = 0
        while ((j < iteration_limit) and (sigma_inner > sum_limit)):
            j += 1
            #max_cum_score += high
            #min_cum_score += low
            
            # Calculate Pj(i)
            spj = score_probabilities_of_length(sp, spj, low, high, divisor, j)
            
            sum1 = 0.0 # The sum of the probability of negative score of length j
            sum2 = 0.0 # The sum of the probability of positive scores of length j
            for s, p in spj.items():
                if (s < 0):
                    sum1 += p * math.exp(s*my_lambda)
                else:
                    sum2 += p
            
            sigma_inner = sum1 + sum2
            sigma_bar += sigma_inner/j
        
        return sigma_bar
    
    # Convert scoring matrix into a dict (called sfp/sprob) containing the probabilities of each score
    sprob = score_probabilities(matrix)
    sprob_average = matrix_average(matrix) # The average score (expected score) of the matrix (should be <0)
    
    low = min(sprob) # Identify the lowest possible score (must be negative)
    high = max(sprob) # Identify the highest possible score (must be positive?)
    my_range = high - low # Find the range in scores
    
    # Check if input arguments are valid
    if ((my_lambda <= 0) or (H <= 0)):
        # Theory dictates that H and lambda must be positive, so
        # return -1 to indicate an error
        #return -1.0
        raise Exception("H and lambda must be positive")
    
    # Karlin-Altschul theory works only if the expected score is negative
    if (sprob_average >= 0):
        #return -1.0
        raise Exception ("The expected score (average) in the scoring scheme must be negative")
    
    # Find the minimum divisor between all scores
    divisors = {}
    for i1, s1 in enumerate(sprob):
        for i2, s2 in enumerate(sprob):
            if (i1 <= i2):
                divisors[(s1, s2)] = gcd_float(s1, s2)
    divisor = min(divisors.values()) # if high=2, low=-6, then divisor=2
    
    fsprob = full_score_probabilities(sprob, divisor, rounding=2)
    
    # Pj(i) = the probability of a local alignment of length j with score i.
    # There are only finitely many i for which P1(i) > 0.
    # l = min{i | P1(i) > 0} # l = minimum possible score (that is defined)
    # u = max{i | P1(i) > 0} # u = maximum possible score (that is defined)
    
    # the following must hold true:
    #  l < 0
    #  u > 0
    
    # The function P1(i) may be represented by an array p with initial index zero and length=(u-l+1) such that:
    #  P1(i) = p[i-l]
    # In our coding notation:
    #  P1(i) = prob_array[i-low]
    # def P1(i):
    #     '''probability that the alignment score is i'''
    #     return prob_array[i-low]
    # Thus, the function P1(i) is the same as:
    #   P1(i) = fsprob[i]
    
    # If (l != -1) and (u != 1):
    #  sigma_bar = ....
    #
    #  Pj(i) = sum(from k=range(-inf, +inf):  P1(k) * Pj-1(i-k)  )
    #        = sum over all possible k ( [Probability of a single nt alignment of score k] * [Probability of alignment of length=j-1 of score i-k])
    # This assumes that the probability of appending an item to the end of a
    # sequence is independent of what has gone before (hence the probability
    # multiplication).
    #
    # There are only finitely many i for which P1(i) > 0,
    # and thus an inductive argument using the Pj(i) equation shows that for
    # any j, there are only finitely many i for which Pj(i) > 0.
    
    # Therefore for any fixed value of j:
    #   sigma_inner = ...
    # which may be computed in a finite number of operations.
    
    # On the other hand, the sigma_outer computation must approximate
    # the sum over all j = 1, ..., inf.
    # The sum over all j is truncated when terms in the sum become sufficiently small.
    # (It is limited by iteration_limit and sum_limit)
    
    # delta is the greatest common divisor of all scores that have nonzero probability.
    # Usually d = 1.
    
    # Once sigma_bar has been computed, K is obtained from the formula:
    #  K = (delta * lambda * math.exp(-2 * sigma_bar)) / (H * (1 - math.exp(-delta * lambda)))
    
    # The BLASTn code uses:
    #  K = -math.exp(-2*sigma_bar) / ((H/lambda)*nonzero_e_to_the_x_minus_1(-my_lambda))
    
    # Are these two K calculations the same??? YES, they are identical!
    
    if ((high == divisor) or (low == -divisor)):
        # Don't compute sigma_bar
        if ((high == divisor) and (low == -divisor)):
            K = (fsprob[divisor] - fsprob[-divisor])**2 / fsprob[-divisor]
        elif ((high == divisor) and (low != -divisor)):
            K = (H/(divisor*my_lambda)) * (1 - math.exp(-divisor*my_lambda))
        elif ((high != divisor) and (low == -divisor)):
            K = ((my_lambda * (1-math.exp(-divisor*my_lambda))) / (divisor*H)) * sprob_average*sprob_average
    else:
        # Compute sigma_bar
        sigma_bar = calculate_sigma_bar(sprob, my_lambda, low, high, divisor)
        
        # K from Gertz 2005:
        K = (divisor * my_lambda * math.exp(-2 * sigma_bar)) / (H * (1 - math.exp(-divisor * my_lambda)))
        
        # The BLASTn code uses (which is identical):
        #K2 = -math.exp(-2*sigma_bar) / ((H/my_lambda)*nonzero_e_to_the_x_minus_1(-my_lambda))
    
    return K

def expected_score(matrix, compositions):
    '''
    The expected score of a scoring matrix is the sum of its raw scores weighted
    by their frequencies of occurrence. The expected score should always be
    negative.
    '''
    # matrix={('A', 'A'): 1, ('A', 'C'): -1, ...}
    # compositions={'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25}
    expected_score = 0
    for k, v in matrix.items():
        expected_score += compositions[k[0]]*compositions[k[1]]*v
    return expected_score

def calculate_compositions(sequence, kmer=1):
    d = {}
    for k in nucleotides.kmers(kmer):
        d[k] = 0
    for i in range(len(sequence) - (kmer-1)):
        seq = sequence[i:i+kmer].upper() # Make upper-case
        
        for char in seq:
            if ((char != 'A') and (char != 'C') and (char != 'G') and (char != 'T')):
                seqs = nucleotides.disambiguate_iupac(seq)
                seqs_len = len(seqs)
                for s in seqs:
                    d[s] += 1/seqs_len
                break
        else:
            d[seq] += 1
    return d

def normalize_compositions(compositions):
    f = {}
    s = sum(compositions.values())
    if (s == 0):
        for k, v in compositions.items():
            f[k] = 0.0
    else:
        for k, v in compositions.items():
            f[k] = v/s
    return f

def contig_compositions(contigs, kmer=1, frequencies=True):
    d = {}
    for k in nucleotides.kmers(kmer):
        d[k] = 0
    
    if (len(contigs) == 0):
        raise Exception("Variable does not contain sequences when it should!")
    
    for k, v in contigs.items():
        comp = calculate_compositions(v, kmer=kmer)
        for k2, v2 in comp.items():
            d[k2] += v2
    if frequencies:
        return normalize_compositions(d)
    else:
        return d

#def calculate_lambda(match, mismatch, matrix, compositions=None, precision=0.001):
def calculate_lambda(match, mismatch, compositions=None, precision=0.001):
    # to calculate lambda, you need
    #  sequence composition
    #  average match_score
    #  average mismatch_score

    if (compositions == None):
        compositions = {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25}
    # Normalize compositions if necessary
    compositions = normalize_compositions(compositions)

    def score(a, b):
        if (a == b):
            return match
        else:
            return mismatch
    
    expected_score = 0
    for ki, vi in compositions.items():
        for kj, vj in compositions.items():
            expected_score += vi*vj*score(ki, kj)
            #expected_score += vi*vj*matrix[(ki, kj)]
    
    if ((match <= 0) or (expected_score >= 0)):
        # Illegal scores
        raise Exception("Illegal match/mismatch scores. Cannot calculate lambda")
    
    # Initial estimates for lambda
    my_lambda = 1
    high = 2
    low = 0
    
    # Calculate lambda
    while (high - low > precision):
        # Calculate the sum of all normalized scores
        my_sum = 0;
        for ki, vi in compositions.items():
            for kj, vj in compositions.items():
                my_sum += vi*vj*math.exp(my_lambda*score(ki, kj))
                #my_sum += vi*vj*math.exp(my_lambda*matrix[(ki, kj)])
        
        # Refine guess at lambda
        if (my_sum > 1): 
            high = my_lambda;
            my_lambda = (my_lambda + low)/2
        else:
            low = my_lambda
            my_lambda = (my_lambda + high)/2
    
    # Compute target frequency and H
    targetID = 0;
    for k, v in compositions.items():
        targetID += (v**2) * math.exp(my_lambda * match)
    
    # H is the relative entropy of the scoring matrix, and is calculated from
    # normalized scores?
    H = my_lambda * match*targetID + my_lambda*mismatch*(1 -targetID);
    
    return expected_score, my_lambda, H, targetID

def filter_scores(matrix):
    """
    Returns copy of matrix with gap alignment scores removed
    """
    new_matrix = {}
    for k, v in matrix.items():
        if (len(k) == 2):
            new_matrix[k] = v
    return new_matrix

def matrix_averages(scores):
    '''
    Returns:
     * the average score of entire matrix
     * the average match score
     * the average non-match score (mismatch)
    '''
    total_s = 0
    total_n = 0
    
    match_s = 0
    match_n = 0
    
    nonmatch_s = 0
    nonmatch_n = 0
    
    for k, v in scores.items():
        total_s += v
        total_n += 1
        
        if (k[0] == k[1]):
            match_s += v
            match_n += 1
        else:
            nonmatch_s += v
            nonmatch_n += 1
    
    # entire_table_average, match_average, non-match_average
    return total_s/total_n, match_s/match_n, nonmatch_s/nonmatch_n

def test():
    """Test"""
    if (len(sys.argv[1:]) != 3):
        print("USAGE python3 evalues.py query.fasta subject.fasta scores.txt", file=sys.stderr)
        sys.exit(1)
    query_file = sys.argv[1]
    subject_file = sys.argv[2]
    scores_file = sys.argv[3]
    
    qcontigs = utils.old_load_fasta_file(query_file)
    scontigs = utils.old_load_fasta_file(subject_file)
    SCORES = load_scores(scores_file)
    
    # SCORES = {
    #     # Matches
    #     ('A', 'A'): 1,
    #     ('C', 'C'): 1,
    #     ('G', 'G'): 1,
    #     ('T', 'T'): 1,
    #     # Transitions (Purines: A<->G, Pyrimidines: C<->T)
    #     ('A', 'G'): -1,
    #     ('G', 'A'): -1,
    #     ('C', 'T'): -1,
    #     ('T', 'C'): -1,
    #     # Transversions
    #     ('A', 'T'): -1.2, # Weak
    #     ('T', 'A'): -1.2, # Weak
    #     ('C', 'G'): -1.2, # Strong
    #     ('G', 'C'): -1.2, # Strong
    #     ('A', 'C'): -1.2, # Amino
    #     ('C', 'A'): -1.2, # Amino
    #     ('G', 'T'): -1.2, # Keto
    #     ('T', 'G'): -1.2, # Keto
    #     # Insertions
    #     ('-', 'A', 'open'): -4,
    #     ('-', 'C', 'open'): -4,
    #     ('-', 'G', 'open'): -4,
    #     ('-', 'T', 'open'): -4,
    #     ('-', 'A', 'extend'): -2,
    #     ('-', 'C', 'extend'): -2,
    #     ('-', 'G', 'extend'): -2,
    #     ('-', 'T', 'extend'): -2,
    #     # Deletions
    #     ('A', '-', 'open'): -4,
    #     ('C', '-', 'open'): -4,
    #     ('G', '-', 'open'): -4,
    #     ('T', '-', 'open'): -4,
    #     ('A', '-', 'extend'): -2,
    #     ('C', '-', 'extend'): -2,
    #     ('G', '-', 'extend'): -2,
    #     ('T', '-', 'extend'): -2,
    #     # Gaps should never be aligned to each other
    #     ('-', '-', 'open'): -4,
    #     ('-', '-', 'extend'): -2,
    # }
    
    SCORES = filter_scores(SCORES)
    
    # scontigs = {
    #     'contig_1A': 'GCCATACTTTCGGTCAAAATACTCAATACATTGTGTAAGTGCTCGCATCAATGCGGACAAAAAAACGGCCTCATAGTATACCTATGCTACCTAAGATTAG',
    #     'contig_2A': 'TTGCAGAATAGGAAGTGAAGAAAAGTTCCTGCAAACTGTGCACTCGATGTTTAAACACGAGATCGGTAGTTACGACAAGGATAAATCATACCAAATAAGA',
    #     'contig_3A': 'GCTCATCCGAAATCGGCCATTTACAAAATTATTCGGGGCGCTAGTAAATTGTGAGTCATATACATGTCGAGGACAATTTCGCTATTTTAAAATCTTAGTT',
    #     'contig_4A': 'TTCATCCGCTAAGTAATAAGTTTTGCGTTATCGAACTCGAGTCTAATAATATCAGAAAACTCTATATCCTCCTTTTAATCTATGAATTACTCTACAAACA',
    #     'contig_5A': 'TATATTCTCTGTGCTAAAGTATAGTTCCAATATGTATTAGGTAAGTTATCTGGACTAATATTTTATCAAAATAAGAGTAGTTTTGGTTAGCGCCTACTAA',
    #     'contig_6A': 'ATTGAAAAAATATCCTCCTAGTATATCATCCCACAACAATTTCTTATTCCTCCCTTCTTCGTTTATGTGCTAACATCAACGGGACTCTATAGGAGCCTAG',
    #     'contig_7A': 'CATGCTTCCCTCCCACGTTTGTTTATTATTTATCTACTCGTAATCTCATTGATACTCAGAGTTAACACATGTGTTTGCCTGACGAATACCTTGAATAAGG',
    #     'contig_RA': 'GAGGAAGTGTACAGATTGTCTACCAGATTCAGTAATTATACTTAACAAATATACATAGATTATAAGTCTAAATAGTGAGAAGAAAGCCGCGATATAATTA',
    #     'contig_1B': 'ACCTTATATGAATTTATGTATTAACCGGGTAATGACCTGCTGTGAGTTCCCTAGTATGTAGCTTGTAAAGATTCTCAGAGTCGTTAATGTACTTGTACTT',
    #     'contig_2B': 'AAAGCGCTAGAAGTAACAATTACAATAGTTTCGTTAGAATCATTTTAAGTATGAAAGATAATAACGGCAGAGCAGGTTAAAAACTTCAGTTCAACGGCAA',
    #     'contig_3B': 'GCAAATTTTAGTTTATTGATATTCGGAGAATATAAAATAGTTAAATAAGTGTTATCACTCATAAAGTACTTTTAGAGGTCTGGAAAGTTCGTTTGGACAT',
    #     'contig_4B': 'TGGTGCGGGCGGAATGTCGTTTTCGTCTGTGGCAAGAGATGCATAATTTTTCGTGCATGATTTAATACAATGCTGCTGTTTATATAACTCTGAAAGCATG',
    #     'contig_5B': 'ATAAACCAGTTCTTAGTTAGATGTCTTAACAGAAAAAATACTCTAAAAAAAACTGCCGAAGAGATTCAGTGATGTATTTGTGTTTTGTGGTTCATGTTTA',
    #     'contig_6B': 'CATAAGGGGCCTAGAACTTATGTTACATGTAGTAATAGATTATTGTTACAAAAAGCTTAACGATATTATGATCTTTTCTCCATTTATAACTCAATTAAAG',
    #     'contig_7B': 'CTCCCATCGTTACATTGGAATCTTCGTTTTAGATTTGTAGACGGTATTGTCTTTTAAAACACTAGATAATGAGGGGCCTACGAAAAAATAACTTAGTATG',
    #     'contig_RB': 'CAAAACTGTACTATTTGGCTTCAGTCTACTTACGTCTAGGATCCCCCCTCAATTCACAAATTGAAGTATTTTTGAGAAAATAACAGGTTTAGTAAAGATA',
    #     'contig_M': 'CTTGTATAAATATAACATAGTTGATGAAATAACAGGATATGACACTACTAGTGTACGGATGTGTGTTTATAACAAAAGATCGCAGCCTTAGTGCGTGAGC',
    # }
    
    dbcomp = contig_compositions(scontigs, kmer=1, frequencies=True)
    n = real_subject_length(scontigs) # 28605418
    print("Number db sequences: {}".format(len(scontigs)))
    print("DB length: {}".format(n))
    print("DB composition: {}".format(dbcomp))
    print("")
    
    for qname, qseq in qcontigs.items():
        #qseq = 'GCACGCCTAAGCGCTATAAGACGCTACGATACGGAGCTCAGCTAGCGACT' # 50 nt
        #qseq = 'TTATATATATATATAAAATAGTTTAGCCTACGACTACGATCTCACGCAGCAGCTACCTATCGCGAACGCATACTCTACGCGAGCATCACACATGGCCGATTTGTTAATTATGTGGCTATT' # 120 nt
        #qseq = 'AAATTAGCTGAAAAGGCAGACGTATTAACTGTTGAAATCGAACATGTTGATGTTGATGCTAGCTATAGGCGCGCATCGCTCTCGCGCGAGAGAGCGCTCTCGCGCGCGCATATATATAGGAGAGGTATTCCTGTTGCAACTGTTGCTATTAACAATAGTACTAATGCTGCATTGTTGGCA' # 180nt
        raw_score = len(qseq)//3
        qcomp = normalize_compositions(calculate_compositions(qseq))
        
        m = len(qseq)
        # Default NCBI blastn parameters: match=1, mismatch=-2, gapopen=0, gapextension=-2.5
        average, average_match, average_mismatch = matrix_averages(SCORES)
        #ES, LAMBDA, H, TARGETID = calculate_lambda(average_match, average_mismatch, compositions['A'], compositions['C'], compositions['G'], compositions['T'])
        ES, LAMBDA, H, TARGETID = calculate_lambda(average_match, average_mismatch, dbcomp)
        K = calculate_k(SCORES, LAMBDA, H)
        normalized_score = LAMBDA * raw_score
        print("Query name: {}".format(qname))
        print("Query length: {}".format(m))
        
        print("Raw score: {}".format(raw_score))
        print("Normalized score: {}".format(normalized_score))
        print("Expected score: {}".format(ES))
        print("Lambda: {} nats ({} bits)".format(LAMBDA, LAMBDA/math.log(2)))
        print("H: {} nats ({} bits)".format(H, H/math.log(2)))
        print("K: {}".format(K))
        print("Expected %ID: {}".format(TARGETID*100))
        
        eHSP = expected_HSP_length(K, m, n, H)
        
        em = effective_query_length(m, eHSP, K)
        en = effective_db_length(n, eHSP, len(scontigs), K)
        ess = search_space(em, en)
        
        E_VALUE = raw_score_to_expect(raw_score, K, LAMBDA, em, en)
        P_VALUE = e_value_to_p_value(E_VALUE)
        
        print("Expected HSP length: {}".format(eHSP))
        print("Effective query length: {}".format(em))
        print("Effective subject length: {}".format(en))
        print("e-value: {}".format(E_VALUE))
        print("p-value: {}".format(P_VALUE))
        print("Effective search space: {}".format(ess))
        
        r1 =  expected_raw_score_if_e_value_is_1(K, em, en, LAMBDA)
        print("Expected raw score if e-value = 1: {}".format(r1))
        print("e-value if raw score = {}: {}".format(r1, raw_score_to_expect(r1, K, LAMBDA, em, en)))
    
class EstimateVariables(object):
    def __init__(self, scoring_matrix, scontigs):
        self.dbcomp = contig_compositions(scontigs, kmer=1, frequencies=True)
        self.n = real_subject_length(scontigs) # 28605418
        
        self.SCORES = filter_scores(scoring_matrix)
        self.number_scontigs = len(scontigs)
        
        self.average, self.average_match, self.average_mismatch = matrix_averages(self.SCORES)
        self.ES, self.LAMBDA, self.H, self.TARGETID = calculate_lambda(self.average_match, self.average_mismatch, self.dbcomp)
        self.K = calculate_k(self.SCORES, self.LAMBDA, self.H)
    
    #def calculate_evalue(match, mismatch, gap_open, gap_extension, alignment_score, alignment_nonclipped_length):
    def calculate_evalue(self, alignment_score, query_sequence):
        # Traditionally,
        #  E = k*m*n* e**-(lambda*S)
        #    or
        #  E = k*m*n* 2**-(lambda*S)
        # Traditional E-value scores ungapped alignments
        
        # H: The 'hit', or alignment
        # S(H): The score (in bits or nats) of the alignment
        # Normalized score: S'(H) = lambda * S(H)
        # P-value: The probability of obtaining at least one hit in a search space
        #          of a given size with a bit-score as large or larger than S'(H), in a random model
        
        # Parameter definitions
        #   E-value: Expected number of alignments with a score at least S'(H) given the query 'm' and a search size of 'n' (in a random model)
        #         k: A small adjustment (k) takes into account the fact that optimal local
        #            alignment scores for alignments that start at different places in the two sequences may be highly correlated. For
        #            example, a high-scoring alignment starting at residues 1, 1 implies a pretty high alignment score for an alignment
        #            starting at residues 2, 2 as well. The value of k is usually around 0.1, and its impact on the statistics of alignment
        #            scores is relatively minor, so don't bother with its derivation.
        #         m: length of query
        #         n: length of subject (database)
        #  lambda*S: the normalized score (of the alignment)
        
        m = len(query_sequence)
        
        eHSP = expected_HSP_length(self.K, m, self.n, self.H)
        em = effective_query_length(m, eHSP, self.K)
        en = effective_db_length(self.n, eHSP, self.number_scontigs, self.K)
        E_VALUE = raw_score_to_expect(alignment_score, self.K, self.LAMBDA, em, en)
        
        return E_VALUE

def test2():
    if (len(sys.argv[1:]) != 2):
        print("USAGE python3 evalues.py subject.fasta scores.txt", file=sys.stderr)
        sys.exit(1)
    subject_file = sys.argv[1]
    scores_file = sys.argv[2] # algorithms\distance_scores.txt
    
    scontigs = utils.old_load_fasta_file(subject_file)
    score_matrix = load_scores(scores_file)
    
    ev = EstimateVariables(score_matrix, scontigs)
    seq = 'ACAGGTGTGTGGCTCTCAGAGCGCGCCTCTCCGATCCAGCAATAGGAATATTGTGAGGATCCGGATTTAC'
    for raw_score in [10, 20, 30, 40]:
        for i in [10, 20, 30, 40, 50, 60, 70]:
            print(raw_score, i, ev.calculate_evalue(raw_score, seq[:i]))
    
    print('===================')
    seq = 'CAATATATGAGAACTGGTGAAGG'
    raw_score = -2 # cigar_string = '1X1=1X1=1X7=1X3=1X5=1X'
    print(raw_score, len(seq), ev.calculate_evalue(raw_score, seq)) # evalue = 425224026.9218612
    raw_score = 46 # cigar_string = '23='
    print(raw_score, len(seq), ev.calculate_evalue(raw_score, seq)) # evalue = 2.671432045203441e-05

if (__name__ == '__main__'):
    #test()
    test2()
