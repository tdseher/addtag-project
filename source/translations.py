#!/usr/bin/env python3
import regex

class CodonTable(object):
    # table = {N: (name, codons, initiations), ...}
    tables = load_tables_file("translations.txt")
    
    def __init__(self, table):
        """Create a new instance of a the CodonTable object"""
        self.name = self.tables[table][0]
        self.table = self.tables[table][1]
        self.initiation = self.tables[table][2]
    
    def codon_translate(self, sequence):
        """
        Returns the input DNA sequence translated with start, stop, and initiation sites
        over a single reading frame only.
        Input:
          DNA sequence
        Output:
          Returns a list of lists: aa, nt, nt_start, nt_stop, aa_initiation_indices
          Each element is an ORF (with no length filter)
        """
        # Returns a new element in the return list for every stop codon '*' found?
        
        ret_list = []
        r = ['', '', 1, 0, []]
        
        for c in range(0, len(sequence), 3):
            nt = sequence[c:c+3]
            
            try:
                aa = self.table[nt]
            except KeyError:
                # X = any amino acid
                if (len(nt) == 3):
                    aa = 'X'
                else:
                    break
            
            r[0] += aa
            r[1] += nt
            r[3] = c + 3
            if (nt in self.initiation):
                r[4].append(len(r[0]) - 1)
            if (aa == '*'):
                ret_list.append(r)
                r = ['', '', c + 1 + 3, 0, []]
        
        if (len(r[0]) > 0):
            ret_list.append(r)
        
        return ret_list
    
    def get_orfs(self, nt_sequence, minimum_aa_length=80, initiation_trim=True, ignore_gaps=True, require_initiation=True):
        """
        Cannot handle nucleotide IUPAC ambiguities
        """
        
        # sequence =            NNN NNN NNN NNN NNN NNN NNN NNN NNN
        # contig_frames[frame] = X   *   X   M   X   X   X   *   X
        # init_indices_frames[frame] = [9]
        # fragments =          [X, X M X X X, X]
        
        if ignore_gaps:
            seq = nt_sequence.replace('-', '')
            seq = seq.replace('.', '')
            seq = seq.replace(' ', '')
        else:
            seq = nt_sequence
        seq = seq.upper()
        
        sequences = {
            1: seq,
            2: seq[1:],
            3: seq[2:],
            -1: self.rc(seq),
            -2: self.rc(seq)[1:],
            -3: self.rc(seq)[2:]
        }
        
        orfs = {
            1: [],
            2: [],
            3: [],
            -1: [],
            -2: [],
            -3: []
        }
        
        for frame in sorted(sequences):
            #print frame
            put_orf_list = self.codon_translate(sequences[frame])
            #for p in put_orf_list:
            #    print frame, p
            
            for aa, nt, nt_start, nt_stop, aa_indices in put_orf_list:
                if (initiation_trim and (len(aa_indices) > 0)):
                    short_aa = aa[aa_indices[0]:]
                    nt_start += 3 * aa_indices[0]
                    short_nt = nt[3 * aa_indices[0]:]
                else:
                    short_aa = aa
                    short_nt = nt
                
                # correct for frame:
                if (frame > 0):
                    nt_start += frame - 1
                    nt_stop += frame - 1
                elif (frame < 0):
                    nt_start, nt_stop = nt_stop, nt_start
                    nt_start = len(seq) - nt_start + frame + 2
                    nt_stop = len(seq) - nt_stop + frame + 2
                
                if (len(short_aa) >= minimum_aa_length):
                    if require_initiation:
                        if (short_nt[:3] in self.initiation):
                            orfs[frame].append([short_aa, short_nt, nt_start, nt_stop])
                    else:
                        orfs[frame].append([short_aa, short_nt, nt_start, nt_stop])
        
        return orfs
    
    @staticmethod
    def rc(sequence):
        """Returns the reverse complement of the input DNA sequence"""
        complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB') # 'w','s', and 'n' are left out because their rc() is the same
        rcseq = sequence.translate(complements)[::-1]
        return rcseq
    
    @classmethod
    def load_tables_file(cls, filename):
        with open(filename, 'r') as flo:
            lines = []
            for i, line in enumerate(flo):
                if (i % 7 == 0):
                    t, name, codons, initiations = cls.load_table(lines)
                    cls.tables[t] = (codons, initiations)
                    lines = []
                else:
                    lines.append(line.rstrip())
    
    @classmethod
    def load_table(cls, lines):
        # Parse the definition line
        m = regex.match(r'^(\d+)\.\s+([^(]*).*$', lines[0])
        t, name = m.groups()
        t = int(t)
        
        # Parse
        d = {}
        for line in lines[1:]:
            m = regex.match(r'^\s*([^= ]*)\s*=\s*(.*)$', line)
            label, data = m.groups()
            d[label] = data
        
        for i in range(64):
            c =  d["Base1"][i] + d["Base2"][i] + d["Base3"][i]
            codons[c] = d["AAs"][i]
            initiations[c] = d[""]
        
        return (t, name, codons, initiations)

def test():
    """Code to test the classes and functions in 'translations.py'"""
    
    print("=== translation test ===")
    seq = "GCACTACATAGTAATTTTTTGTATGGGAAAAAATGTATGGCAATYTTTTGTATGGTAATTTTTTGTG\
          TGAGAATACATAGTAATTTTTTGTATGGTAAAAAATACATAGTAAAAAATGTATGGTAATTTTTTGTA\
          TGKGTGCACTACATAGTAATTTTTTGTATGGGAAAAAATGCTAGTAATTTTTTTGTATGGCAATTTTT\
          TGTGTGCGTGYAATACATAGTAATTTTTTGTATGGTAAAAAATACTAGTAATTTTTTTGTGTGTGTAG\
          AAAAAACACTAACATTTTTTTGTATGGTAATTTTTTGTGTGTGAGAATACTAGTAATTTTTTGTATGT\
          GTGCACTACATARTAATTTTTTGTATGGGAAAAAATACATAGTAAAAAATGTATGGCAATTTTTTTGT\
          GTGTGTACTACATAGTAAATTTTTGTGTGTGTGAGAATACATAGTAATTTTTTGTATGGGAAAAAATA\
          CATAGTAATTTTTTTGTATGGTGTATTATGTAATATTTTTTTGTATGGCAATTTTTTGTATGCGTATA\
          ATACATAGTAAATTTTTTGTATGGCAATTTTTTTGTATGGTATTTTTTTGTATGGGAAAAAACCCTAA\
          CATTTTTTGTATGGTGATTTTTTTGTATGGCAATTTTTATTGTGTGGCAATTTTTTTGTATGGTGTTT\
          TTTCTTGTGTTCGATCCCAGCTTCTCCTTATTTATACACCCTCTGGTTCTCCCCCTTTTTGTTAGACA\
          CCACCCAACTTTATCACCATCACTGGATAGAGCAGGGTCGATTACCCTTGGTTCATGCAAAAAAAAAA\
          TATACCTTACTAAAAAAAATAGATTGCAGCACAATAGTTTCGCGTATGGTCTCCCACTACACTACTCG\
          GTATTGCTCTTAGCAGCTTAACTACGGTTGATCGAACGGGGAACGGTGCTTTCTGCTAGATATGGCCG\
          CAACCGAAAAAATAAAAAAATGCGAACCTTATCTGCTGTCTCTCACGTGACCAAATTCACTTTTAGGG\
          CTCTAGCCCCACAACCACCAAAGTGTATGTGCTGTCGCTGCAGGGGAGGGGTAATCGGGGTGCCCAGA\
          ATTGTGTGGAGCCATTTTTTGAGCCGGAAAGTTGGGTGGCTGTGGCACAAAACGGGAATATGTATGTC\
          CGGGTGGCCAGTTGACTGGGTTGTAGCAGTGCTGCAGCTACGAATGTTAGAGACAAAATGTAGTCCAG\
          GGCCGGCCAGACGCAGTGTGCGTGCTGGTTGTGCAATTATACTAGCACATCTAGGTGTTTTAGGAAGA\
          AACGTCCACCACCAAAAAATTAATTTCCAAATATTGGCACTTTTATTACACCAGTGGTGTTACACAGC\
          CCCCAAAACACTATCCCGGTTGTTAGATTGAAGTTTGTCTAACAAAAATTGGAAGTTCTATTTTTTAC\
          TTTTTGTACAAAATTTGGCAAGAAAATTGAAGTAAAATATTTTTATATAAATTTCAAACTAAAAATAT\
          CCAAAACATAAAATAATACACAAAATGAGATGCTCTATTTAGCCAAACCAACCATTACGGGTGTGTTG\
          TTTGGTGTTGTCTGACCATGGGTATACCATTTGTTAGTGTATAGCTGCACTGTTTTATGTTTTTCGAT".replace(" ", "")
    
    # Make the list of codon-translation tables
    tables = []
    for t in [1, 3, 12]:
        tables.append(CodonTable(t))
    
    # Iterate through the CodonTable objects
    for t in tables:
        # identify all ORFs using this CodonTable
        orfs = t.get_orfs(seq, minimum_length=60)
        
        # iterate through the identified ORFs
        for f in sorted(orfs): # f is the frame [-3, -2, -1, 1, 2, 3]
            for o in orfs[f]: # o looks like this: ['RELTTTHG*', 'CGCGAGCTAACTACCACTCATGGGTAG', 20, 46]
                print(f, o)
    

if (__name__ == '__main__'):
    test()