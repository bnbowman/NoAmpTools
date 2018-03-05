#! /usr/bin/env python

## Copyright (c) 2018, Pacific Biosciences of California, Inc.
##
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted (subject to the limitations in the
## disclaimer below) provided that the following conditions are met:
##
##  * Redistributions of source code must retain the above copyright
##    notice, this list of conditions and the following disclaimer.
##
##  * Redistributions in binary form must reproduce the above
##    copyright notice, this list of conditions and the following
##    disclaimer in the documentation and/or other materials provided
##    with the distribution.
##
##  * Neither the name of Pacific Biosciences nor the names of its
##    contributors may be used to endorse or promote products derived
##    from this software without specific prior written permission.
##
## NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
## GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
## BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
## WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
## OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
## CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
## USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
## ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
## OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
## OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
## SUCH DAMAGE.

## Author: Brett Bowman

import sys
import re
from collections import defaultdict

from pbcore.io import IndexedBamReader, PacBioBamIndex, IndexedFastaReader, FastaRecord

import ConsensusCore2 as cc

MIN_ACC   = 0.8
MIN_T     = 0.35

if len(sys.argv) != 4:
    print "ERROR:\tExpected 3 arguments but got {0}".format(len(sys.argv)-1)
    print "Usage:\tloadingDiagnostic HG19.FASTA ALIGN_BAM_PBI SCRAPS_BAM"
    raise SystemExit

indexedFasta = sys.argv[1]
inputFile    = sys.argv[2]
scrapsBam    = sys.argv[3]

cfg = cc.AlignConfig(cc.AlignParams.Default(), 1);

## Locus,ChrName,ChrIdx,GeneStart,RegionStart,RegionEnd,GeneEnd
TARGETS = [["HTT", "chr4", 3, 3075691, 3076603, 3076661, 3076815],
           ["FMR1", "chrX", 23, 146993123, 146993568, 146993629, 146994131],
           ["ALS", "chr9", 8, 27572985, 27573522, 27573541, 27574014],
           ["SCA10", "chr22", 21, 46190744, 46191234, 46191305, 46191756]]

GUIDES = {"C9orf72a": "GCAATTCCACCAGTCGCTAG",
          "C9orf72b": "GCATGATCTCCTCGCCGGCA",
          "FMR1"    : "AGAGGCCGAACTGGGATAAC",
          "HTT"     : "AGCGGGCCCAAACTCACGGT",
          "ATXN10"  : "ATACAAAGGATCAGAATCCC"}

def ScoreCas9Site( seq ):
    maxKey = None
    maxAcc = None
    for key, rna in GUIDES.iteritems():
        query = rna + "NGG"
        aln = cc.Align(seq, query, cfg)
        Ns = sum(1 for b in aln.Query() if b == 'N')
        acc = (aln.Matches() + Ns) / float(len(query))
        if maxAcc is None or acc > maxAcc:
            maxKey = key
            maxAcc = acc
    return (maxKey, maxAcc)

def ScoreCas9SiteSides( outSeq, inSeq ):
    k1, a1 = ScoreCas9Site( outSeq )
    k2, a2 = ScoreCas9Site( inSeq )
    if max([a1, a2]) < MIN_ACC:
        return ("N/A", "N/A", "N/A")
    elif a1 >= a2:
        return (k1, "OUT", a1)
    else:
        return (k2, "IN",  a2)

def HasEcoR1( seq ):
    return "T" if ("GAATTC" in seq) else "F"

def HasBamH1( seq ):
    return "T" if ("GGATCC" in seq) else "F"

def LargestAs( seq ):
    grps = [len(m.group(0)) for m in re.finditer(r"(\w)\1*", seq) if m.group(0)[0] == "A"]
    grps = [g for g in grps if g >= 10]
    return sum(grps)

def LargestTs( seq ):
    grps = [len(m.group(0)) for m in re.finditer(r"(\w)\1*", seq) if m.group(0)[0] == "T"]
    grps = [g for g in grps if g >= 10]
    return sum(grps)

def LargestAsAndTs( seq ):
    grps = [len(m.group(0)) for m in re.finditer(r"(\w)\1*", seq) if m.group(0)[0] in ["T", "A"]]
    grps = [g for g in grps if g >= 10]
    return grps

def ReadGenomeWindowsFromPBI( fns, tList ):
    # Conver the target-list to a dictionary for faster searching
    cov = defaultdict(int)
    acc = defaultdict(float)
    windows = {}
    for fn in fns:
        pbi = PacBioBamIndex( fn )
        hnIdx     = pbi.columnNames.index("holeNumber")
        tIdIdx    = pbi.columnNames.index("tId")
        tStartIdx = pbi.columnNames.index("tStart")
        tEndIdx   = pbi.columnNames.index("tEnd")
        matIdx    = pbi.columnNames.index("nM")
        missIdx   = pbi.columnNames.index("nMM")
        delIdx    = pbi.columnNames.index("nDel")
        insIdx    = pbi.columnNames.index("nIns")
        mapQvIdx  = pbi.columnNames.index("mapQV")
        ctxIdx    = pbi.columnNames.index("contextFlag")

        for row in pbi:
            # Skip secondary alignments
            if row[mapQvIdx] == 0:
                continue

            flag    = row[ctxIdx]
            if not ((flag & 1) and (flag & 2)):
                continue

            nM      = row[matIdx]
            nMM     = row[missIdx]
            nIns    = row[insIdx]
            nDel    = row[delIdx]
            tAcc    = nM / float(nM + nMM + nIns + nDel)

            # Track which ZMWs we've seen
            hn      = int(row[hnIdx])
            tId     = row[tIdIdx]
            tStart  = row[tStartIdx]
            tEnd    = row[tEndIdx]
            tCov    = tEnd - tStart

            target = "OFF"
            for tName, _, tTid, _, tRS, tRE, _ in tList:
                if tTid != tId:
                    continue
                elif tStart < tRS and tEnd > tRE:
                    target = tName
                    break

            if tCov > cov[hn]:
                cov[hn] = tCov
                windows[hn] = (hn, tId, tStart, tEnd, target)

    return sorted(v for k,v in windows.iteritems())

def ReadAdaptersFromScraps( bam ):
    adps  = defaultdict(int)    
    polyA = defaultdict(int)
    with IndexedBamReader( bam ) as handle:
        for record in handle:
            if record.scrapType != "A":
                continue
            hn  = record.holeNumber
            seq = record.peer.seq
            adps[hn] += 1
            tFrac = sum(1 for b in seq if b == "T") / float(len(seq))
            if tFrac > MIN_T:
                polyA[hn] += 1

    # Convert our counts into a T/F depending on whether there are polyAs
    res = {}
    for hn, v in adps.iteritems():
        if v >= 2:
            res[hn] = "T" if polyA[hn] >= 1 else "F"
    return res

def SortWindowsByChromosome( windows ):
    byChrom = defaultdict(list)
    for hn, tId, tS, tE, tTarg in windows:
        byChrom[tId].append( (hn, tS, tE, tTarg) )
    return byChrom

def FindOverlaps( sortedWin ):
    retval = {}
    for ch, data in sortedWin.iteritems():
        ovls = []
        if ch != 3:
            continue
        currE = data[0][1]
        currOvl = [data[0]]
        for hn, s, e, t in data[1:]:
            if s < currE:
                currOvl.append( (s, e, t) )
            else:
                if len(currOvl) > 1000:
                    ovls.append( currOvl )
                currOvl = [ (hn, s, e, t) ]
                currE   = e
        retval[ch] = ovls
    return retval

# Second, tabulate the number of usable reads/ZMWs
windows = ReadGenomeWindowsFromPBI( [inputFile], TARGETS )
adps    = ReadAdaptersFromScraps( scrapsBam )
#byChrom = SortWindowsByChromosome( windows )
#ovls = FindOverlaps( byChrom )

print "HoleNumber,Chromosome,Start,End,Target,PolyAAdp,PolyARegion,MaxPolyARegion,TotalPolyARegion,LeftEcoR1,LeftBamH1,RightEcoRI,RightBamH1,LeftRna,LeftRnaSide,LeftRnaAcc,RightRna,RightRnaSide,RightRna,HasPolyA,HasLeft,HasRight"
fa = IndexedFastaReader( indexedFasta )
for hn, tid, s, e, target in windows:
    # First skip ZMWs with no adp results, i.e. with <= 1 adp
    try:
        polyA = adps[hn]
    except:
        continue

    chrm = fa[tid]

    # Search for restriction sites near
    fiveP  = chrm.sequence[s-5:s+6]
    threeP = chrm.sequence[e-5:e+6]
    fiveEco, fiveBam   = HasEcoR1(fiveP),  HasBamH1(fiveP)
    threeEco, threeBam = HasEcoR1(threeP), HasBamH1(threeP)

    # Count and summarize any PolyA/T regions
    region  = chrm.sequence[s:e]
    AT = LargestAsAndTs( region )
    maxAT = 0 if len(AT) == 0 else max(AT)

    # Check for Guide RNA matches
    OutFiveP  = chrm.sequence[s-33:s+10]
    InFiveP   = FastaRecord("tmp", chrm.sequence[s-10:s+33]).reverseComplement().sequence
    InThreeP  = chrm.sequence[e-33:e+10]
    OutThreeP = FastaRecord("tmp", chrm.sequence[e-10:e+33]).reverseComplement().sequence
    k1, s1, a1 = ScoreCas9SiteSides( OutFiveP,  InFiveP )
    k2, s2, a2 = ScoreCas9SiteSides( OutThreeP, InThreeP )

    # Summary columns
    hasPolyA = "T" if (polyA == "T" or maxAT > 0) else "F"
    hasLeft  = "T" if (fiveEco == "T" or fiveBam == "T" or k1 != "N/A") else "F"
    hasRight = "T" if (threeEco == "T" or threeBam == "T" or k2 != "N/A") else "F"

    print "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21}".format(hn, tid, s, e, target, polyA, len(AT), maxAT, sum(AT), fiveEco, fiveBam, threeEco, threeBam, k1, s1, a1, k2, s2, a2, hasPolyA, hasLeft, hasRight)
