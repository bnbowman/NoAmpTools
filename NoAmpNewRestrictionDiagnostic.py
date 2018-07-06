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
import json
from collections import defaultdict
import string

import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt

from pbcore.io import IndexedBamReader, PacBioBamIndex, IndexedFastaReader, FastaRecord, openDataSet

COMPLEMENT = string.maketrans("ACGTacgt-", "TGCAtgca-")

if len(sys.argv) < 4:
    print "ERROR:\tExpected at least 3 arguments but got {0}".format(len(sys.argv)-1)
    print "Usage:\tNoAmpNewRestrictionDiagnostic.py OUTPUT_PREFIX HG19.FASTA ALIGN_BAM [ALIGN_BAM ..]"
    raise SystemExit

outputPrefix = sys.argv[1]
indexedFasta = sys.argv[2]
inputFiles    = sys.argv[3:]

## Locus,ChrName,ChrIdx,GeneStart,RegionStart,RegionEnd,GeneEnd
TARGETS = [["HTT", "chr4", 3, 3075691, 3076603, 3076661, 3076815],
           ["FMR1", "chrX", 23, 146993123, 146993568, 146993629, 146994131],
           ["ALS", "chr9", 8, 27572985, 27573522, 27573541, 27574014],
           ["SCA10", "chr22", 21, 46190744, 46191234, 46191305, 46191756]]

def HasCutSite( re_site, seq ):
    return "T" if (re_site in seq) else "F"

def HasEcoR1( seq ):
    return HasCutSite("GAATTC", seq)

def HasBamH1( seq ):
    return HasCutSite("GGATCC", seq)

def HasEcoR5( seq ):
    return HasCutSite("GATATC", seq)

def HasSpe1( seq ):
    return HasCutSite("ACTAGT", seq)

def HasAvr2( seq ):
    return HasCutSite("CCTAGG", seq)

def HasBssSa1( seq ):
    return HasCutSite("CACGAG", seq)

def HasAcc1( seq ):
    # GTMKAC = [GTAGAC GTATAC GTCGAC GTCTAC]
    sites = ["GTAGAC", "GTATAC", "GTCGAC", "GTCTAC"]
    res = [HasCutSite(site, seq) for site in sites]
    return "T" if ("T" in res) else "F"

def HasKpn1( seq ):
    return HasCutSite("GGTACC", seq)

def HasBgl2( seq ):
    return HasCutSite("AGATCT", seq)

def HasSexA1( seq ):
    # ACCWGGT = [ACCAGGT ACCTGGT]
    sites = ["ACCAGGT", "ACCTGGT"]
    res = [HasCutSite(site, seq) for site in sites]
    return "T" if ("T" in res) else "F"

def ParseAdapterTypes( record ):
    # Try to read the AD tag in this record, if any
    try:
        opt = record.peer.opt("ad")
    except:
        return (None, None)

    # If we found an AD tag, parse it
    left, right = opt.split(";")
    if left == ".":
        leftType = None
    else:
        leftType = int(left.split(',')[0])
    if right == ".":
        rightType = None
    else:
        rightType = int(right.split(',')[0])

    return (leftType, rightType)

def ReadAlignedBamFile( fns, tList ):
    # Dictionaries for tracking adapter results
    adps    = defaultdict(int)    
    altAdps = defaultdict(int)

    # Dictionaries for tracking ZMW-level results
    cov = defaultdict(int)
    windows = {}
    for fn in fns:
        for record in IndexedBamReader(fn):
            # Skip secondary alignments
            if record.MapQV == 0:
                continue

            # We have nothing to learn from subreads with no adapters
            leftT, rightT = ParseAdapterTypes(record)
            if leftT is None and rightT is None:
                continue
            
            hn = record.holeNumber

            # Adapter stats are duplicated between the subreads that sandwich them, so
            #  we arbitrarily pick the right side here to avoid double-counts when 
            #  recording information about our adapters
            if rightT is not None:
                adps[hn] += 1
                if rightT == 1:
                    altAdps[hn] += 1

            # If we made it past all of our filters, parse the rest of the data we want:
            tId    = record.tId
            tStart = record.tStart
            tEnd   = record.tEnd
            tCov   = tEnd - tStart

            # Search our target list for targets that overlap our current subread
            target = "OFF"
            for tName, _, tTid, _, tRS, tRE, _ in tList:
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

    # Convert our dictionary of windows to a flat, sorted list
    windowResults = sorted(v for k,v in windows.iteritems())

    # Convert our adapter counts into a T/F depending on whether there are polyAs
    adpResults = {}
    for hn, v in adps.iteritems():
        if v >= 2:
            adpResults[hn] = "T" if altAdps[hn] >= 1 else "F"

    # Return a tuple containing both our cleaned up Window and Adapter results
    return (windowResults, adpResults)

def SearchSequence( seq ):
    res = {"EcoR1" : HasEcoR1( seq ),
           "EcoR5" : HasEcoR5( seq ),
           "BamH1" : HasBamH1( seq ),
           "Spe1"  : HasSpe1( seq ),
           "Acc1"  : HasAcc1( seq ),
           "Kpn1"  : HasKpn1( seq ),
           "Bgl2"  : HasBgl2( seq ),
           "Avr2"  : HasAvr2( seq ),
           "BssSa1": HasBssSa1( seq ),
           "SexA1": HasSexA1( seq )}
    return res

def ResultsToName( reDict ):
    # Return "Mixed" if >1 to avoid double-counting
    if reDict.values().count("T") > 1:
        return "Mixed"
    # Return the name of the RE we hit, if any
    for k, v in reDict.iteritems():
        if v == "T":
            return k
    # If we got this far we have no RE hits
    return "None"

def SummarizeRestrictionData( indexedFasta, windows, adps ):
    """Summarize the data for each ZMW, and their left and right sides"""
    results = []

    fa = IndexedFastaReader( indexedFasta )
    for hn, tid, s, e, target in windows:
        # First skip ZMWs with no adp results, i.e. with <= 1 adp
        try:
            polyA = adps[hn]
        except:
            continue
        chrm = fa[tid]

        # Search for restriction sites near the alignment edges
        fiveP  = chrm.sequence[s-5:s+6]
        threeP = chrm.sequence[e-5:e+6]
        threeP_rc = threeP.translate(COMPLEMENT)[::-1]
        left = SearchSequence( fiveP )
        right = SearchSequence( threeP_rc )

        results.append( (hn, tid, s, e, target, left, right) )

    return sorted(results)

def WriteSummaryCsv( outputPrefix, summaries ):
    """Write a summary of our results to CSV for downstream QC"""
    with open( outputPrefix.lower() + ".cut_sites.csv", 'w') as handle:
        handle.write("HoleNumber,Chromosome,Start,End,Target,LeftEcoR1,LeftEcoRV,LeftBamH1,LeftSpe1,LeftAcc1,LeftKpn1,LeftBgl2,LeftAvr1,LeftBssSa1,LeftSexA1,")
        handle.write("RightEcoRI,RightEcoRV,RightBamH1,RightSpe1,RightAcc1,RightKpn1,RightBgl2,RightAvr1,RightBssSa1,RightSexA1\n")

        for hn, tid, s, e, target, left, right in summaries:
            handle.write("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24}\n".format(hn, tid, s, e, target,
                                                                                            left['EcoR1'],  left['EcoR5'],  left['BamH1'],  left['Spe1'],  left['Acc1'],  left['Kpn1'],  left['Bgl2'],  left["Avr2"],  left["BssSa1"],  left["SexA1"],
                                                                                            right['EcoR1'], right['EcoR5'], right['BamH1'], right['Spe1'], right['Acc1'], right['Kpn1'], right['Bgl2'], right["Avr2"], right["BssSa1"], right["SexA1"]))

def TabulateRestrictionTable( summaries ):
    """Tabulate our cut-sites combinations and plot them as a table"""
    counts = defaultdict(lambda : defaultdict(int))
    for hn, tid, s, e, target, left, right in summaries:
        # Skip unclear ZMWs before tabulation to avoid double-counting
        leftName = ResultsToName( left )
        rightName = ResultsToName( right )
        if leftName == "Mixed" or rightName == "Mixed":
            continue
        counts[leftName][rightName] += 1
    return counts

def PlotRestrictionCountsTable( outputPrefix, counts ):
    """Tabulate our cut-sites combinations and plot them as a table"""
    # Convert those counts to an array-of-arrays with sums
    ids = ["None"] + sorted(SearchSequence("NNNNNNNNNNN").keys())
    rows = []
    for leftRE in ids:
        row = []
        for rightRE in ids:
            row.append( counts[leftRE][rightRE] )
        rows.append( [leftRE] + row + [sum(row)] )

    finalRow = ["Sum"] + [0] * (len(rows[0]) - 1)
    for row in rows:
        for i, val in enumerate(row):
            if i > 0:
                finalRow[i] += val
    rows.append( finalRow )

    # Plot the results as a table
    fig = plt.figure(frameon=False, figsize=(8, 6.5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    t = ax.table(cellText=rows,
                colLabels=["Left \ Right"] + ids + ["Sum"],
                loc='center', cellLoc='center')
    t.set_fontsize(24)
    t.scale(1, 3.0)
    pltFilename = "{0}_re_counts.png".format(outputPrefix.lower())
    plt.savefig(pltFilename, bbox='tight')

    p = {"caption": "Table of Restriction Enzyme Cut-Site Counts",
           "image": pltFilename,
            "tags": [],
              "id": "{0} - Restriction-Enzyme Counts".format(outputPrefix),
           "title": "{0} - RestrictionEnzymeCounts".format(outputPrefix)}
    return p

def PlotRestrictionFracsTable( outputPrefix, counts, precision=2 ):
    """Tabulate our cut-sites combinations and plot them as a table"""
    ids = ["None"] + sorted(SearchSequence("NNNNNNNNNNN").keys())
    total = 0.000001
    for leftRE in ids:
        if leftRE == "Mixed": continue
        for rightRE in ids:
            if rightRE == "Mixed": continue
            total += counts[leftRE][rightRE]

    # Convert those counts to an array-of-arrays with sums
    rows = []
    for leftRE in ids:
        row = []
        for rightRE in ids:
            row.append( round(100 * counts[leftRE][rightRE] / total, precision) )
        rows.append( [leftRE] + row + [round(sum(row), precision)] )

    finalRow = ["Sum"] + [0] * (len(rows[0]) - 1)
    for row in rows:
        for i, val in enumerate(row):
            if i > 0:
                finalRow[i] += val
    for i, val in enumerate(finalRow):
        if i > 0:
            finalRow[i] = round(finalRow[i], precision)
    rows.append( finalRow )

    # Plot the results as a table
    fig = plt.figure(frameon=False, figsize=(8, 6.5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    t = ax.table(cellText=rows,
                colLabels=["Left \ Right"] + ids + ["Sum"],
                loc='center', cellLoc='center')
    t.set_fontsize(24)
    t.scale(1, 3)
    pltFilename = "{0}_re_fractions.png".format(outputPrefix.lower())
    plt.savefig(pltFilename, bbox='tight')

    p = {"caption": "Table of Restriction Enzyme Cut-Site Fractions",
           "image": pltFilename,
            "tags": [],
              "id": "{0} - Restriction-Enzyme Fractions".format(outputPrefix),
           "title": "{0} - RestrictionEnzymeFractions".format(outputPrefix)}
    return p

def WriteReportJson( plotList=[], tableList=[] ):
    reportDict = {"plots":plotList, "tables":tableList}
    reportStr = json.dumps(reportDict, indent=1)
    with open("report.json", 'w') as handle:
        handle.write(reportStr)


# Second, tabulate the number of usable reads/ZMWs
windows, adps = ReadAlignedBamFile( inputFiles, TARGETS )
summaries = SummarizeRestrictionData( indexedFasta, windows, adps )
WriteSummaryCsv( outputPrefix, summaries )
counts = TabulateRestrictionTable( summaries )
p1 = PlotRestrictionCountsTable( outputPrefix, counts )
p2 = PlotRestrictionFracsTable( outputPrefix, counts )
WriteReportJson( [p1, p2] )
