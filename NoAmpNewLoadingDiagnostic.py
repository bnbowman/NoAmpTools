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

import json
import re
from collections import defaultdict

import numpy as np
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt

import seaborn as sns

from pbcore.io import IndexedBamReader, PacBioBamIndex, IndexedFastaReader, FastaRecord, openDataSet

import ConsensusCore2 as cc

MIN_ACC   = 0.8

if len(sys.argv) < 4:
    print "ERROR:\tExpected at least 3 arguments but got {0}".format(len(sys.argv)-1)
    print "Usage:\tloadingDiagnostic OUTPUT_PREFIX HG19.FASTA ALIGN_BAM [ALIGN_BAM ..]"
    raise SystemExit

outputPrefix = sys.argv[1]
indexedFasta = sys.argv[2]
inputFiles   = sys.argv[3:]

cfg = cc.AlignConfig(cc.AlignParams.Default(), 1);

## Locus,ChrName,ChrIdx,GeneStart,RegionStart,RegionEnd,GeneEnd
TARGETS = [["HTT", "chr4", 4, 3074344, 3076603, 3076661, 3077078],
           ["FMR1", "chrX", 23, 146993123, 146993568, 146993629, 146994131],
           ["ALS", "chr9", 9, 27572985, 27573522, 27573541, 27574014],
           ["FUCHS", "chr18", 18, 53251995, 53253386, 53253458, 53253577],
           ["SCA10", "chr22", 22, 46190744, 46191234, 46191305, 46191756],
           ["EWINGS_Chr20", "chr20", 20, 21553989, 21556922, 21557001, 21557036],
           ["EWINGS_ChrX", "chrX", 23, 30325813, 30328875, 30328976, 30329062]]

GUIDES = {"FMR1"     : "AGAGGCCGAACTGGGATAAC",
          "FMR1_201" : "CGCGCGTCTGTCTTTCGACC",
          "HTT"      : "AGCGGGCCCAAACTCACGGT",
          "HTT_SQ1"  : "CTTATTAACAGCAGAGAACT"}

def TargetsToTargetDict( targets ):
    tDict = defaultdict(list)
    for t in targets:
        if t[2] <= 22:
            target_tId = t[2]-1
        else:
            target_tId = t[2]
        tDict[target_tId].append( t )
    return tDict

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

def ParseSingleAdapterCounts( string ):
    # All return values are tuples of [TC6-count, AltAdp-count]
    if string == ".":
        return ["F", "F"]
    
    adpId = int(string.split(',')[0])
    if adpId == 0:
        return ["T", "F"]
    elif adpId == 1:
        return ["F", "T"]
    
    # We shouldn't get here, but return 0,0 for unknown if we do
    return ["F", "F"]

def ParseAdapterTypes( record ):
    # All return values are tuples of [Left-TC6-count, Left-AltAdp-count, Right-TC6-count, Right-AltAdp-count]

    # Try to read the AD tag in this record, if any
    try:
        opt = record.peer.opt("ad")
    except:
        return ["F", "F", "F", "F"]

    # If we found an AD tag, parse the left and right sides
    left, right = opt.split(";")
    leftCounts  = ParseSingleAdapterCounts( left )
    rightCounts = ParseSingleAdapterCounts( right )

    # Concatenate and return left-side and right-side results, ordered by genomic position
    if record.isReverseStrand:
        return rightCounts + leftCounts
    else:
        return leftCounts + rightCounts

def ReadAlignedBamFile( fns, tDict ):
    # Dictionaries for tracking ZMW-level results
    cov = defaultdict(int)
    adps = {}
    windows = {}
    for fn in fns:
        for record in IndexedBamReader(fn):
            # Skip secondary alignments
            if record.MapQV == 0:
                continue

            hn     = record.holeNumber
            tId    = record.tId
            tStart = record.tStart
            tEnd   = record.tEnd
            tCov   = tEnd - tStart
            
            adpTypes = ParseAdapterTypes(record)

            # Search our target list for targets that overlap our current subread
            target = "OFF"
            for tName, _, tTid, _, tRS, tRE, _ in tDict[tId]:
                if tTid != tId:
                    continue
                elif tStart < tRS and tEnd > tRE:
                    target = tName
                    break

            # If our coverage for this subread is better than anything we've already seen
            #  for this ZMW, keep it instead
            if tCov > cov[hn]:
                cov[hn] = tCov
                windows[hn] = (hn, tId, tStart, tEnd, target)
                adps[hn] = adpTypes

    # Return a tuple containing both our cleaned up Window and Adapter results
    return (windows, adps)

def SummarizeData( indexedFasta, windows, adps ):
    summaries = []

    fa = IndexedFastaReader( indexedFasta )
    for hn, (_, tid, s, e, target) in windows.iteritems():
        # First skip ZMWs with no adp results, i.e. with <= 1 adp
        try:
            leftTc6, leftAlt, rightTc6, rightAlt = adps[hn]
        except:
            continue

        chrm = fa[tid]

        # Search for restriction sites near the ends
        fiveP    = chrm.sequence[max(s-5, 0):s+6]
        threeP   = chrm.sequence[e-5:e+6]
        fiveEco  = HasEcoR1(fiveP)
        threeEco = HasEcoR1(threeP)

        # Search for restriction sites contained within
        inside    = chrm.sequence[s+6:e-5]
        insideEco = HasEcoR1(inside)

        # Count and summarize any PolyA/T regions
        region  = chrm.sequence[s:e]
        AT = LargestAsAndTs( region )
        maxAT = 0 if len(AT) == 0 else max(AT)

        # Check for Guide RNA matches
        OutFiveP  = chrm.sequence[max(s-33, 0):s+10]
        InFiveP   = FastaRecord("tmp", chrm.sequence[max(s-10,0):s+33]).reverseComplement().sequence
        InThreeP  = chrm.sequence[e-33:e+10]
        OutThreeP = FastaRecord("tmp", chrm.sequence[e-10:e+33]).reverseComplement().sequence
        k1, s1, a1 = ScoreCas9SiteSides( OutFiveP,  InFiveP )
        k2, s2, a2 = ScoreCas9SiteSides( OutThreeP, InThreeP )

        # Summary columns
        hasPolyA = "T" if maxAT > 0 else "F"
        hasLeft  = "T" if (fiveEco == "T" or k1 != "N/A") else "F"
        hasRight = "T" if (threeEco == "T" or k2 != "N/A") else "F"

        summaries.append( (hn, tid, s, e, e-s, target, len(AT), maxAT, sum(AT), leftTc6, rightTc6, leftAlt, rightAlt, fiveEco, insideEco, threeEco, k1, s1, a1, k2, s2, a2, hasPolyA, hasLeft, hasRight) )

    return sorted(summaries)

def WriteSummaryCsv( outputPrefix, summaries ):
    with open(outputPrefix.lower() + ".loading.csv", 'w') as handle:
        handle.write("HoleNumber,Chromosome,Start,End,InsertSize,Target,PolyARegion,MaxPolyARegion,TotalPolyARegion,LeftAdpTc6,RightAdpTc6,LeftAdpAlt,RightAdpAlt,LeftEcoR1,InsideEcoR1,RightEcoR1,LeftRna,LeftRnaSide,LeftRnaAcc,RightRna,RightRnaSide,RightRna,HasPolyA,HasLeft,HasRight\n")
        for row in summaries:
            handle.write(",".join(str(s) for s in row) + "\n")

def PlotAdapterEcoR1Table( outputPrefix, summaries ):
    counts = defaultdict(int)
    total = 0.000001
    for row in summaries:
        leftTc6, leftAlt, rightTc6, rightAlt = row[9:13]
        leftEcoR1, insideEcoR1, rightEcoR1 = row[13:16]
        left  = "ALT" if leftAlt  == "T" else "TC6"
        right = "ALT" if rightAlt == "T" else "TC6"
        ecoR1 = [leftEcoR1, rightEcoR1].count("T")
        t = sorted([left, right], reverse=True) + [ecoR1]
        summaryStr = "{0}:{1} ({2}x EcoR1)".format(*t)
        counts[summaryStr] += 1
        total += 1

    cumsum = 0
    rows = []
    for k1 in ["TC6", "ALT"]:
        for k2 in ["TC6", "ALT"]:
            if k1 == "ALT" and k2 == "TC6":
                continue
            for k3 in [0, 1, 2]:
                summaryStr = "{0}:{1} ({2}x EcoR1)".format(k1, k2, k3)
                count = counts[summaryStr]
                cumsum += count
                rows.append( [summaryStr, count, "{}%".format(round(100 * count / total, 2))] )
    rows.append( ["Sum", cumsum, "{}%".format(round(100 * cumsum / total, 2))] )

    # Plot the results as a table
    fig = plt.figure(frameon=False, figsize=(6, 5.5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    t = ax.table(cellText=rows,
                colLabels=["Molecule Ends", "Count", "Fraction"],
                loc='center', cellLoc='center')
    t.set_fontsize(24)
    t.scale(1, 3)
    pltFilename = "{0}_adapter_pairs.png".format(outputPrefix.lower())
    plt.savefig(pltFilename, bbox='tight')
    plt.close()

    p = {"caption": "Table of Adapter & EcoR1 Counts",
           "image": pltFilename,
            "tags": [],
              "id": "{0} - Adapter & EcoR1 Counts".format(outputPrefix),
           "title": "{0} - AdapterEcoR1Counts".format(outputPrefix)}
    return p

def PlotAdapterOnTargetTable( outputPrefix, summaries ):
    counts = defaultdict(int)
    total = 0.000001
    for row in summaries:
        target = "False" if row[5] == "OFF" else "True"
        leftTc6, rightTc6, leftPolyA, rightPolyA = row[9:13]
        leftEcoR1, insideEcoR1, rightEcoR1 = row[13:16]
        left  = "ALT" if leftPolyA  == "T" else "TC6"
        right = "ALT" if rightPolyA == "T" else "TC6"
        ecoR1 = [leftEcoR1, rightEcoR1].count("T")
        if left == "ALT" and right == "ALT":
            continue
        elif ecoR1 == 1 and ((left == "ALT") ^ (right == "ALT")):
            counts[("TC6:ALT (1x EcoR1)", target)] += 1
        elif ecoR1 == 2 and left == "TC6" and right == "TC6":
            counts[("TC6:TC6 (2x EcoR1)", target)] += 1
        else:
            counts["OTHER"] += 1
        total += 1

    cumsum = 0
    rows = []
    for (k1, k2) in [("ALT", 1), ("TC6", 2)]:
        for k3 in ["False", "True"]:
            summaryStr = "TC6:{0} ({1}x EcoR1)".format(k1, k2)
            count = counts[(summaryStr, k3)]
            cumsum += count
            rows.append( [summaryStr, k3, count, "{}%".format(round(100 * count / total, 2))] )
    # Append a final line combining every other category
    otherCt = counts["OTHER"]
    cumsum += otherCt
    rows.append( ["Other", "", otherCt, "{}%".format(round(100 * otherCt / total, 2))] )
    # Add the final summation row
    rows.append( ["Sum", "", cumsum, "{}%".format(round(100 * cumsum / total, 2))] )

    # Plot the results as a table
    fig = plt.figure(frameon=False, figsize=(6, 3.5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    t = ax.table(cellText=rows,
                colLabels=["Molecule Ends", "Target", "Count", "Fraction"],
                loc='center', cellLoc='center')
    t.set_fontsize(24)
    t.scale(1, 3)
    pltFilename = "{0}_adapter_ontarget.png".format(outputPrefix.lower())
    plt.savefig(pltFilename, bbox='tight')
    plt.close()

    p = {"caption": "Table of Adapter & OnTarget Counts",
           "image": pltFilename,
            "tags": [],
              "id": "{0} - Adapter & OnTarget Counts".format(outputPrefix),
           "title": "{0} - AdapterOnTargetCounts".format(outputPrefix)}
    return p

def PlotInternalEcoR1Count( outputPrefix, summaries ):
    counts = defaultdict(int)
    total = 0.000001
    for row in summaries:
        insideEcoR1 = row[14]
        counts[insideEcoR1] += 1
        total += 1

    rows = []
    rows.append( ["True", counts["T"], "{}%".format(round(100 * counts["T"] / total, 2))] )
    rows.append( ["False", counts["F"], "{}%".format(round(100 * counts["F"] / total, 2))] )
    rows.append( ["Sum", counts["T"] + counts["F"], "{}%".format(round(100 * (counts["T"] + counts["F"]) / total, 2))] )

    # Plot the results as a table
    fig = plt.figure(frameon=False, figsize=(6, 2.04))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    t = ax.table(cellText=rows,
                colLabels=["Internal EcoR1 Site", "Count", "Fraction"],
                loc='center', cellLoc='center')
    t.set_fontsize(24)
    t.scale(1, 3)
    pltFilename = "{0}_internal_ecoR1.png".format(outputPrefix.lower())
    plt.savefig(pltFilename, bbox='tight')
    plt.close()

    p = {"caption": "Table of Internal EcoR1 Counts",
           "image": pltFilename,
            "tags": [],
              "id": "{0} - Internal EcoR1 Counts".format(outputPrefix),
           "title": "{0} - InternalEcoR1Counts".format(outputPrefix)}
    return p

def PlotInsertSizeHistogram( outputPrefix, summaries ):
    sizes = defaultdict(list)
    for row in summaries:
        size = row[4]
        target = row[5]
        sizes[target].append( size )

    # Convert our sizes into a jagged array, starting with off-targets
    labels = ["OFF"]
    sizeList = [np.array(sizes["OFF"])]
    for k in sorted(sizes.keys()):
        print k, len(sizes[k])
        if k != "OFF":
            labels.append( k )
            sizeList.append( np.array(sizes[k]) )

    # Plot the results as a table
    sns.kdeplot(sizes["OFF"], shade=True, label="OFF")
    for k in sorted(sizes.keys()):
        if k != "OFF" and len(sizes[k]) >= 10:
            sns.kdeplot(sizes[k], shade=True, label=k)
    plt.xlim(0, 8000)
    plt.ylim(0, 0.001)
    pltFilename = "{0}_insert_sizes.png".format(outputPrefix.lower())
    plt.savefig(pltFilename)

    p = {"caption": "Distribution of Insert Sizes",
           "image": pltFilename,
            "tags": [],
              "id": "{0} - Insert Size Distribution".format(outputPrefix),
           "title": "{0} - InsertSizeDistribution".format(outputPrefix)}
    return p

def WriteReportJson( plotList=[], tableList=[] ):
    reportDict = {"plots":plotList, "tables":tableList}
    reportStr = json.dumps(reportDict, indent=1)
    with open("report.json", 'w') as handle:
        handle.write(reportStr)

# Second, tabulate the number of usable reads/ZMWs
tDict = TargetsToTargetDict( TARGETS )
windows, adps = ReadAlignedBamFile( inputFiles, tDict )
summaries = SummarizeData( indexedFasta, windows, adps )
WriteSummaryCsv( outputPrefix, summaries )
p1 = PlotAdapterEcoR1Table( outputPrefix, summaries )
p2 = PlotAdapterOnTargetTable( outputPrefix, summaries )
p3 = PlotInternalEcoR1Count( outputPrefix, summaries )
p4 = PlotInsertSizeHistogram( outputPrefix, summaries )
WriteReportJson( [p1, p2, p3, p4] )
