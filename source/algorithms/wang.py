#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/wang/__init__.py

# Metric for "efficiency"
#  Wang
#   Range: 0-100. SVM model trained on human cell culture data on guides
#   from >1000 genes. The Xu score can be considered an updated version of
#   this score, as the training data overlaps a lot. Delivery: lentivirus.
#   See Wang et al.: http://http//www.ncbi.nlm.nih.gov/pmc/articles/PMC3972032/

def seqToVec(seq, offsets={"A":0,"C":1,"G":2,"T":3}):
    """ convert a x bp sequence to a 4 * x 0/1 vector
    >>> seqToVec("AAAAATTTTTGGGGGCCCCC")
    [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0]
    """
    assert(len(seq)==20)
    row = [0]*len(seq)*4
    pseudoOffset = offsets["A"]
    for pos, nucl in enumerate(seq):
        nucl = nucl.upper()
        # treat N, Y, etc like "A". Happens very rarely.
        nuclOffset = offsets.get(nucl, pseudoOffset)
        vecPos = (pos*len(offsets))+nuclOffset
        #if vecPos not in range(len(row)):
            #ofh = open("temp.txt", "a")
            #ofh.write(str(vecPos)+" "+seq+" "+str(row)+"pos %d, nucl %s" % (pos, nucl)+"\n")
            #assert(False)
        row[vecPos] = 1
    return row

def listToSvml(vec, res):
    """ convert a list of values to a line in svml format line like "+1 1:0.5 2:1.5 ...
    """
    parts = [str(res)]
    for i, val in enumerate(vec):
        parts.append("%d:%d" % (i+1, val))
    return " ".join(parts)

def calcWangSvmScores(seqs):
    """
    Use the wang.model file to score sequences. Input is only the 20bp guide sequence.
    Uses libsvm's svm-predict program, V2.6.
    The score is inversed, so higher scores are better, like all other scores

    The results here are off mostly by 1-5% from the results returned by the Wang et al source code.
    I never found out why, there are no parameters for "svm_predict". Should not be due to a version
    difference either, I'm using the same libsvm version as the e1071 R module.
    This is necessary for a web server implementation as e1071 in R cannot read the model from a file.

    The original implementation from the paper can be called with calcWangSvmScoresUsingR()

    See compareWangScores.py:

    The Pearson correlation between both ways to calculate the score is 97%.

    Histogram of the score differences:
    0.000000 ************************************************************ 3074
    0.050000 ********************************* 1674
    0.100000 ************ 612
    0.150000 **** 191
    0.200000 * 52
    0.250000  7
    0.300000  1
    cat out/wangDiffs.tsv | cut -f4 | tr -d '-' | grep -v diff | textHistogram stdin stdout -real -binSize=0.05

    >>> calcWangSvmScores(["ATAGACCTACCTTGTTGAAG"])
    [60]
    >>> calcWangSvmScores(["NTAGACCTACCTTGTTGAAG"])
    [60]
    """
    scores = []
    vecOrder = {"A":0, "C":1, "T":2, "G":3}

    lines = []
    for seq in seqs:
        seq = seq.upper()
        assert(len(seq)==20)
        vec = seqToVec(seq, offsets=vecOrder)
        lines.append(listToSvml(vec, 0))

    dataIn = "\n".join(lines)
    binPath = getBinPath("svm-predict")
    modelFname = join(binDir, "src", "wangSabatiniSvm", "wang.model")
    cmd = [binPath, "-b", "1", "/dev/stdin", modelFname, "/dev/stdout"]
    proc = Popen(cmd,stdout=PIPE, stdin=PIPE, stderr=None, bufsize=BUFSIZE)
    dataOut = proc.communicate(input=dataIn)[0]

    lines = dataOut.splitlines()
    for line in lines:
        if line.startswith("labels"):
            continue
        if line.startswith("Accuracy"):
            break
        score = int(100*(1.0 - float(line.split()[-1])))
        scores.append(score)

    return scores

def test():
    pass

if (__name__ == "__main__"):
    test()
