#!/usr/bin/env python3

import regex
import random
import time

def levenshtein_distance(seq1, seq2):
    '''
    A fast and memory efficient implementation
    by Hjelmqvist, Sten
    '''
    # degenerate cases
    if seq1 == seq2:
        return 0
    if len(seq1) == 0:
        return len(seq2)
    if len(seq2) == 0:
        return len(seq1)
  
    # create two work vectors of integer distances
    #int[] v0 = new int[seq2.Length + 1];
    #int[] v1 = new int[seq2.Length + 1];
    v0 = []
    v1 = []
  
    # initialize v0 (the previous row of distances)
    # this row is A[0][i]: edit distance for an empty seq1
    # the distance is just the number of characters to delete from seq2
    # for (int i = 0; i < v0.Length; i++)
    # v0[i] = i;
    for i in range(len(seq2)+1):
        v0.append(i)
        v1.append(0)
 
    for i in range(len(seq1)): 
        # calculate v1 (current row distances) from the previous row v0
        # first element of v1 is A[i+1][0]
        # edit distance is delete (i+1) chars from seq1 to match empty seq2
        v1[0] = i + 1
  
        # use formula to fill in the rest of the row
        for j in range(len(seq2)):
            cost = 0 if seq1[i] == seq2[j] else 1;
            v1[j + 1] = min(v1[j]+1, v0[j+1]+1, v0[j]+cost)
  
        # copy v1 (current row) to v0 (previous row) for next iteration
        for j in range(len(seq2)+1):
            v0[j] = v1[j]
  
    return v1[len(seq2)]

def random_motif_sequence(motif_string):
    # Convert the IUPAC sequence to an equivalent regex
    iupac = {
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
        
        '.': ['A', 'C', 'G', 'T'],
    }
    
    seq = ''
    
    matches = regex.findall(r'(.)(?:\{(\d*)(?:(,)(\d*))?\})?', motif_string)
    # example: 'ACN{2}CCN{6,8}G{3,}A{,4}'
    # matches = [('A', '', '', ''), ('C', '', '', ''), ('N', '2', '', ''), ('C', '', '', ''), ('C', '', '', ''), ('N', '6', ',', '8'), ('G', '3', ',', ''), ('A', '', ',', '4')]
    for m in matches:
        qmin, qmax = 1, 1
        if m[1]:
            qmin = int(m[1])
        if m[3]:
            qmax = int(m[3])
        if qmin > qmax:
            qmax = qmin
        
        if (m[0] not in ['>', '<', '/', '\\', '|']):
            n = random.randint(qmin, qmax)
            for j in range(n):
                seq += random.choice(iupac[m[0]])
    
    return seq
    

def generate_min_distance(motif_string, dist, max_fails=1000, max_successes=200):
    seqs = []
    distances = {}
    
    cur_fails = 0
    cur_successes = 0
    while ((cur_fails < max_fails) and (cur_successes < max_successes)):
        # Create a random sequence according to the motif
        seq = random_motif_sequence(motif_string)
        
        # Add current random sequence if it is sufficient distance from all other previously-generated sequences
        for s in seqs:
            current_dist = distances.setdefault((seq, s), levenshtein_distance(seq, s))
            if (current_dist < dist):
                cur_fails += 1
                break
        else:
            seqs.append(seq)
            cur_successes += 1
        
            #if (len(seqs) % 10 == 0):
            #    print(len(seqs), flush=True)
        
    return seqs, cur_fails, cur_successes

if (__name__ == '__main__'):
    
    print("\t".join(["edit distance", "sequence length", "time (seconds)", "fails", "successes"]))
    
    for edit_distance in [0, 1, 2, 3]:
        for seq_length in [5, 10, 15, 20]:
            motif = 'N{'+str(seq_length)+'}'
            start = time.time()
            x, x_fails, x_successes = generate_min_distance(motif, edit_distance, max_successes=200)
            print("\t".join([str(y) for y in [edit_distance, seq_length, time.time()-start, x_fails, x_successes]]))
