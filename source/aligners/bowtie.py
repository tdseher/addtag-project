#!/usr/bin/env python3

"""AddTag Copyright (c) 2016 Thaddeus D. Seher & Aaron Hernday"""

# source/aligners/bowtie.py

# List general Python imports
import sys
import os
import subprocess

# import non-standard package
import regex

# import AddTag-specific packages
from .. import utils

def index_reference():
    pass

def align():
    pass

# Code from CRISPOR
def findOfftargetsBowtie(queue, batchId, batchBase, faFname, genome, pamPat, bedFname):
    " align guides with pam in faFname to genome and write off-targets to bedFname "
    
    # maximum number of occurences in the genome to get flagged as repeats. 
    # This is used in bwa samse, when converting the same file
    # and for warnings in the table output.
    MAXOCC = 60000
    
    # the BWA queue size is 2M by default. We derive the queue size from MAXOCC
    MFAC = 2000000/MAXOCC
    
    # minimum off-target score of standard off-targets (those that end with NGG)
    # This should probably be based on the CFD score these days
    # But for now, I'll let the user do the filtering
    MINSCORE = 0.0
    
    # minimum off-target score for alternative PAM off-targets
    # There is not a lot of data to support this cutoff, but it seems
    # reasonable to have at least some cutoff, as otherwise we would show
    # NAG and NGA like NGG and the data shows clearly that the alternative
    # PAMs are not recognized as well as the main NGG PAM.
    # so for now, I just filter out very degenerative ones. the best solution
    # would be to have a special penalty on the CFD score, but CFS does not 
    # support non-NGG PAMs (is this actually true?)
    ALTPAMMINSCORE = 1.0
    
    
    tmpDir = batchBase+".bowtie.tmp"
    os.mkdir(tmpDir)

    # make sure this directory gets removed, no matter what
    global tmpDirsDelExit
    tmpDirsDelExit.append(tmpDir)
    if not DEBUG:
        atexit.register(delTmpDirs)

    # write out the sequences for bowtie
    queue.startStep(batchId, "seqPrep", "preparing sequences")
    bwFaFname = abspath(join(tmpDir, "bowtieIn.fa"))
    guideSeqs, qSeqs, allPamSeqs = writeBowtieSequences(faFname, bwFaFname, pamPat)

    genomePath =  abspath(join(genomesDir, genome, genome))
    oldCwd = os.getcwd()

    # run bowtie
    queue.startStep(batchId, "bowtie", "aligning with bowtie")
    os.chdir(tmpDir) # bowtie writes to hardcoded output filenames with --refout
    # -v 3 = up to three mismatches
    # -y   = try hard
    # -t   = print time it took
    # -k   = output up to X alignments
    # -m   = do not output any hit if a read has more than X hits
    # --max = write all reads that exceed -m to this file
    # --refout = output in bowtie format, not SAM
    # --maxbts=2000 maximum number of backtracks
    # -p 4 = use four threads
    # --mm = use mmap
    maxOcc = MAXOCC # meaning in BWA: includes any PAM, in bowtie we have the PAM in the input sequence
    cmd = "$BIN/bowtie -e 1000 %(genomePath)s -f %(bwFaFname)s  -v 3 -y -t -k %(maxOcc)d -m %(maxOcc)d dummy --max tooManyHits.txt --mm --refout --maxbts=2000 -p 4" % locals()
    runCmd(cmd)
    os.chdir(oldCwd)

    queue.startStep(batchId, "parse", "parsing alignments")
    pamLen = len(pamPat)
    hits = parseRefout(tmpDir, guideSeqs, pamLen)

    queue.startStep(batchId, "scoreOts", "scoring off-targets")
    # make the list of alternative PAM sequences
    altPats = offtargetPams.get(pamPat, [])
    altPamSeqs = []
    for altPat in altPats:
        altPamSeqs.extend(expandIupac(altPat))

    # iterate over bowtie hits and write to a BED file with scores
    # if the hit looks OK (right PAM + score is high enough)
    tempBedPath = join(tmpDir, "bowtieHits.bed")
    tempFh = open(tempBedPath, "w")

    offTargets = {}
    for guideIdWithMod, chrom, start, end, strand, _, tSeq in hits:
        guideId = guideIdWithMod.split(".")[0]
        guideSeq = guideSeqs[guideId]
        genomePamSeq = tSeq[-pamLen:]
        logging.debug( "PAM seq: %s of %s" % (genomePamSeq, tSeq))
        if genomePamSeq in altPamSeqs:
            minScore = ALTPAMMINSCORE
        elif genomePamSeq in allPamSeqs:
            minScore = MINSCORE
        else:
            logging.debug("Skipping off-target for %s: %s:%d-%d" % (guideId, chrom, start, end))
            continue

        logging.debug("off-target minScore = %f" % minScore )

        # check if this match passes the off-target score limit
        if pamIsCpf1(pamPat):
            otScore = 0.0
        else:
            tSeqNoPam = tSeq[:-pamLen]

            otScore = calcHitScore(guideSeq, tSeqNoPam)

            if otScore < minScore:
                logging.debug("off-target not accepted")
                continue

        editDist = getEditDist(guideSeq, tSeqNoPam)
        guideHitCount = 0
        guideId = guideId.split(".")[0] # full guide ID looks like s33+.0:A>T
        name = guideId+"|"+strand+"|"+str(editDist)+"|"+tSeq+"|"+str(guideHitCount)+"|"+str(otScore)
        row = [chrom, str(start), str(end), name]
        # this way of collecting the features will remove the duplicates
        otKey = (chrom, start, end, strand, guideId)
        logging.debug("off-target key is %s" % str(otKey))
        offTargets[ otKey ] = row

    for rowKey, row in offTargets.iteritems():
        tempFh.write("\t".join(row))
        tempFh.write("\n")

    tempFh.flush()

    # create a tempfile which is moved over upon success
    # makes sure we do not leave behind a half-written file if 
    # we crash later
    tmpFd, tmpAnnotOffsPath = tempfile.mkstemp(dir=tmpDir, prefix="annotOfftargets")
    tmpFh = open(tmpAnnotOffsPath, "w")

    # get name of file with genome locus names
    genomeDir = genomesDir # make local var
    segFname = "%(genomeDir)s/%(genome)s/%(genome)s.segments.bed" % locals()

    # annotate with genome locus names
    cmd = "$BIN/overlapSelect %(segFname)s %(tempBedPath)s stdout -mergeOutput -selectFmt=bed -inFmt=bed | cut -f1,2,3,4,8 > %(tmpAnnotOffsPath)s" % locals()
    runCmd(cmd)

    shutil.move(tmpAnnotOffsPath, bedFname)
    queue.startStep(batchId, "done", "Job completed")

    if DEBUG:
        logging.info("debug mode: Not deleting %s" % tmpDir)
    else:
        shutil.rmtree(tmpDir)

def cleanup():
    pass

def test():
    """Code to test the classes and functions in 'source/bowtie.py'"""
    
    print("=== index_reference ===")
    
    print("=== align ===")
    

if (__name__ == '__main__'):
    test()
