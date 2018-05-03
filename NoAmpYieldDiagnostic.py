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
import os.path
import json
from collections import defaultdict

import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from pbcore.io import IndexedBamReader, PacBioBamIndex

BIN_SIZE  = 2000000
MIN_ALIGN = 400
MIN_SNR   = 3.75
MIN_BQ    = 40
N_COLOR   = 9

outputPrefix  = sys.argv[1]
inputFiles    = sys.argv[2:]

## Chromosome Labels from HG19
LABELS = ["chr{0}".format(l) for l in range(1,23) + ['M', 'X', 'Y']]

## Chromosome Sizes for HG19
SIZES = {"chr1": 249250621, "chr2":  243199373, "chr3":  198022430, "chr4":  191154276,
         "chr5": 180915260, "chr6":  171115067, "chr7":  159138663, "chr8":  146364022,
         "chr9": 141213431, "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
        "chr13": 115169878, "chr14": 107349540, "chr15": 102531392, "chr16":  90354753,
        "chr17":  81195210, "chr18":  78077248, "chr19":  59128983, "chr20":  63025520,
        "chr21":  48129895, "chr22":  51304566,
         "chrM":     16571,  "chrX": 155270560, "chrY":   59373566}
## Duplicate each entry, so we can access it by either name or index
for i, c in enumerate(LABELS):
    SIZES[i] = SIZES[c]

## Locus,ChrName,ChrIdx,GeneStart,RegionStart,RegionEnd,GeneEnd
TARGETS = [["HTT", "chr4", 4, 3075691, 3076603, 3076661, 3076815],
           ["FMR1", "chrX", 23, 146993123, 146993568, 146993629, 146994131],
           ["ALS", "chr9", 9, 27572985, 27573522, 27573541, 27574014],
           ["FUCHS", "chr18", 18, 53251995, 53253386, 53253458, 53253577],
           ["SCA10", "chr22", 22, 46190744, 46191234, 46191305, 46191756],
           ["EWINGS_Chr20", "chr20", 20, 21553989, 21556922, 21557001, 21557036],
           ["EWINGS_ChrX", "chrX", 23, 30325813, 30328875, 30328976, 30329062]]

def GetIndexFiles( fns ):
    indexFiles = []
    for fn in fns:
        if fn.endswith(".pbi"):
            if os.path.isfile(fn):
                indexFiles.append( fn )
            else:
                raise SystemExit("Invalid input file! Cannot access '{0}'".format(fn))
        elif fn.endswith(".bam"):
            if os.path.isfile(fn + ".pbi"):
                indexFiles.append( fn + ".pbi" )
            else:
                raise SystemExit("Input file not indexed! Cannot access '{0}.pbi'".format(fn))
        else:
            raise SystemExit("Unrecognized file format! Cannot read '{0}'".format(fn))
    return indexFiles

def GetColorPalatte( sizes ):
    colors = {}
    for i, n in enumerate( sorted(sizes.keys()) ):
        color = cm.gist_rainbow(i%(N_COLOR+1) / float(N_COLOR))
        colors[i] = color
        colors[n] = color
    return colors

def ReadCoverageDataFromPBI( fns ):
    raw = {}
    zmws = {}
    for fn in fns:
        pbi = PacBioBamIndex( fn )
        hnIdx     = pbi.columnNames.index("holeNumber")
        tIdIdx    = pbi.columnNames.index("tId")
        tStartIdx = pbi.columnNames.index("tStart")
        tEndIdx   = pbi.columnNames.index("tEnd")
        mapQvIdx  = pbi.columnNames.index("mapQV")

        for row in pbi:
            # Skip secondary alignments
            if row[mapQvIdx] == 0:
                continue

            # Don't trust and skip small alignments
            tStart = row[tStartIdx]
            tEnd   = row[tEndIdx]
            tCov   = tEnd - tStart
            if tCov < MIN_ALIGN:
                continue

            # Tabulate the number of Subreads/Zmw by bin
            hn     = row[hnIdx]
            tId    = row[tIdIdx]
            tMid   = (tStart + tEnd) / 2
            tBin   = tMid / BIN_SIZE
            try:
                raw[tId][tBin] += 1
                zmws[tId][tBin].add( hn )
            except:
                raw[tId] = defaultdict(int)
                raw[tId][tBin] += 1
                zmws[tId] = defaultdict(set)
                zmws[tId][tBin].add( hn )

    # Convert the ZMW data from dicts-of-sets into a dicts-of-ints
    zmwCounts = {}
    for tId, binDict in zmws.iteritems():
        zmwCounts[tId] = defaultdict(int)
        for tBin, hns in binDict.iteritems():
            zmwCounts[tId][tBin] = len(hns)

    return raw, zmwCounts

def ConvertAlignDictToList( dataDict, sizes ):
    data = []
    for tId in sorted(dataDict.keys()):
        # Skip the mitochondrial reference / ChrM
        if tId == 22:
            continue

        # Add each binned region as a new data row / graph bar
        size = sizes[tId]
        covDict = dataDict[tId]
        for i in range(size / BIN_SIZE):
            data.append( (tId, covDict[i]) )

        # Append a zero at the end as a spacer for the dividing line
        data.append( (tId, 0) )

    return data

def PlotCoverageData( data, colors, targets, name, outputPrefix ):
    # Separate the Chromosome from the coverage information
    chrs = [v[0] for v in data]
    cov  = [v[1] for v in data]

    # Find the edges and middles of each bin
    cutoffs = [p-1 for p in range(1, len(chrs)) if chrs[p] != chrs[p-1]] + [len(chrs)]
    medians = [(cutoffs[i]/2.0 if i == 0 else cutoffs[i-1] + (cutoffs[i]-cutoffs[i-1])/2.0) for i in range(len(cutoffs))]
    labels = list(range(1,23)) + ['X', 'Y']

    plt.figure(figsize=(24,6))
    bars = plt.bar(range(len(cov)), cov, width=1.0)
    for i in range(len(bars)):
        bars[i].set_color(colors[chrs[i]])
    for x in cutoffs:
        plt.axvline(x=x, ls='--', color='black', alpha=0.2)
    plt.xlim(0, len(cov))
    plt.xticks(medians, labels)
    plt.xlabel("Genomic Position")
    plt.ylabel("Alignment Count")
    pltTitle = "{0} - {1} Distribution by Genomic Position".format(outputPrefix, name)
    plt.title(pltTitle)

    # Add target labels to the graph and set the y-axis
    maxCov = max(cov)
    for t in targets:
        tName   = t[0]
        chrIdx = t[2]
        mid    = (t[3] + t[6]) / 2
        chrBin = mid / BIN_SIZE
        offset = len(name) * 4 - 4
        tBin    = cutoffs[chrIdx-2] + chrBin - offset + 2  # +2 to align text-middle rather than text-bottom
        plt.text(tBin, maxCov, tName, fontsize=14, rotation='vertical')

    plt.ylim(0, maxCov * 1.1)
    pltFilename = "{0}_{1}_coverage.png".format(outputPrefix.lower(), name.lower())
    plt.savefig(pltFilename)

    p = {"caption": pltTitle,
           "image": pltFilename,
            "tags": [],
              "id": "{0} - GenomeCoverageBy{1}".format(outputPrefix, name),
           "title": "{0} - Genome Coverage By {1}".format(outputPrefix, name)}
    return p

def ReadOnTargetCountsFromPBI( fns, targets ):
    # Conver the target-list to a dictionary for faster searching
    tDict = defaultdict(list)
    for t in targets:
        if t[2] <= 22:
            target_tId = t[2]-1
        else:
            target_tId = t[2]
        tDict[target_tId].append( t )

    hits = {t[0]:defaultdict(int) for t in targets}
    hns = defaultdict(int)
    bcs = defaultdict(int)

    for fn in fns:
        pbi = PacBioBamIndex( fn )
        hnIdx     = pbi.columnNames.index("holeNumber")
        tIdIdx    = pbi.columnNames.index("tId")
        tStartIdx = pbi.columnNames.index("tStart")
        tEndIdx   = pbi.columnNames.index("tEnd")
        mapQvIdx  = pbi.columnNames.index("mapQV")

        # Get barcode columns if present, otherwise None
        if "bcForward" in pbi.columnNames:
            bcIdx = pbi.columnNames.index("bcForward")
            bqIdx = pbi.columnNames.index("bcQual")
        else:
            bcIdx, bqIdx = None, None

        for row in pbi:
            # Skip secondary alignments
            if row[mapQvIdx] == 0:
                continue

            # Track which ZMWs we've seen
            hn      = int(row[hnIdx])
            tStart  = row[tStartIdx]
            tEnd    = row[tEndIdx]
            tCov    = tEnd - tStart
            hns[hn] = max(tCov, hns[hn])

            # Record the barcode
            if bqIdx and row[bqIdx] >= MIN_BQ:
                bcs[hn] = row[bcIdx]

            # Skip alignments to chromosomes w/ no hits
            tId    = row[tIdIdx]
            if tId not in tDict.keys():
                continue

            # Skip records with poor alignments
            if tCov < MIN_ALIGN:
                continue

            hnMax = "{0}_max".format(hn)
            for tName, _, _, _, rS, rE, _ in tDict[tId]:
                if tStart < rS and tEnd > rE:
                    hits[tName][hn] += 1
                    hits[tName][hnMax] = max(tCov, hits[tName][hnMax])

    # Return our counts
    if len(bcs) == 0:
        return hits, hns, None
    else:
        return hits, hns, bcs

def InvertBarcodeDict( hnCov, bcCalls ):
    res = defaultdict(set)
    for hn in hnCov.keys():
        res["ALL"].add( hn )
        if bcCalls and hn in bcCalls:
            bc = bcCalls[hn]
            res[bc].add( hn )
            res[bc].add( "{0}_max".format(hn) )
    return res

def PlotOnTargetTable( onTargetD, setName, hnSet, nZmw, totalCov, targets, sizes, outputPrefix ):
    expCov = {t[0]:t[6]-t[3] for t in targets}
    gSize  = float(sum(l for chrm, l in sizes.iteritems() if isinstance(chrm, str)))

    tZmw, tEnrich, tCcs, tSubread = 0, 0, 0, 0
    rLabs, rows = [], []

    for tName, hDict in onTargetD.iteritems():
        rLabs.append( tName )
        zmws     = len([hn for hn in hDict.keys() if isinstance(hn, int) and hn in hnSet])
        zmwFrac  = round(100 * zmws / float(nZmw), 3)
        cov      = sum(c for hn, c in hDict.iteritems() if isinstance(hn, str) and hn in hnSet)
        exp      = expCov[tName]
        enrich   = int( (cov / float(totalCov)) / (exp / gSize) )
        subreads = sum(c for hn, c in hDict.iteritems() if isinstance(hn, int) and hn in hnSet)
        ccsEst   = sum(1 for hn, c in hDict.iteritems() if (isinstance(hn, int) and hn in hnSet and c >= 3))
        rows.append( [tName, str(zmws), str(zmwFrac) + "%", str(enrich) + "-fold", str(ccsEst), str(subreads)] )

        # Accumulate sample-wide totals
        tZmw     += zmws
        tEnrich  += enrich
        tCcs     += ccsEst
        tSubread += subreads

    # Add a row of totals at the bottom
    rLabs.append( "Total" )
    tFrac = round(100 * tZmw / float(nZmw), 3)
    tEnrich = tEnrich / len(rows)
    rows.append( ["Total", str(tZmw), str(tFrac) + "%", str(tEnrich) + "-fold", str(tCcs), str(tSubread)] )

    # Plot
    fig = plt.figure(frameon=False, figsize=(8, 4.115))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    t = ax.table(cellText=rows,
                rowLabels=rLabs,
                colLabels=["Target", "nZMW", "FracZMW", "Enrichment", "estCCS", "nSubread"],
                loc='center', cellLoc='center')
    t.set_fontsize(24)
    t.scale(1, 2.72)
    if setName == "ALL":
        pltFilename = "{0}_target_table.png".format(outputPrefix.lower())
    else:
        pltFilename = "{0}_target_table.{1}.png".format(outputPrefix.lower(), setName)
    plt.savefig(pltFilename, bbox='tight')

    p = {"caption": "Table of On-Target Subreads/ZMWs By Locus",
           "image": pltFilename,
            "tags": [],
              "id": "{0} (BC #{1}) - On-Target Table".format(outputPrefix, setName),
           "title": "{0} (BC #{1}) - OnTargetTable".format(outputPrefix, setName)}
    return p

def PlotOnTargetTables( onTargetD, hnCov, bcCalls, targets, sizes, outputPrefix ):
    # Sort our hns by barcode instead of vise-versa
    bcSets = InvertBarcodeDict( hnCov, bcCalls )

    # Iterate over each
    plots = []
    for bc, hnSet in bcSets.iteritems():
        print bc, len(list(hnSet)), sorted(list(hnSet))
        if bc == "ALL":
            nZmw = len(hnCov.keys())
            tCov = sum(cov for hn, cov in hnCov.iteritems())
        else:
            nZmw = len([hn for hn in hnCov.keys() if hn in hnSet])
            tCov = sum(cov for hn, cov in hnCov.iteritems() if hn in hnSet)
        plots.append( PlotOnTargetTable( onTargetD, bc, hnSet, nZmw, tCov, targets, sizes, outputPrefix ) )

    return plots

# First plot the over-all coverage information
indexFiles = GetIndexFiles( inputFiles )
rawD, zmwD = ReadCoverageDataFromPBI( indexFiles )
raw = ConvertAlignDictToList( rawD, SIZES )
zmw = ConvertAlignDictToList( zmwD, SIZES )
colors = GetColorPalatte( SIZES )
p1 = PlotCoverageData( raw, colors, TARGETS, "Subread", outputPrefix)
p2 = PlotCoverageData( zmw, colors, TARGETS, "ZMW", outputPrefix)

# Second, tabulate the number of usable reads/ZMWs
onTarget, hnCov, bcCalls = ReadOnTargetCountsFromPBI( indexFiles, TARGETS )
p3 = PlotOnTargetTables( onTarget, hnCov, bcCalls, TARGETS, SIZES, outputPrefix )

# Finally, combine our plots into a JSON report to output for ZIA
reportDict = {"plots":[], "tables":[]}
reportDict["plots"].append( p1 )
reportDict["plots"].append( p2 )
for p in p3:
    reportDict["plots"].append( p3 )
reportStr = json.dumps(reportDict, indent=1)
with open("report.json", 'w') as handle:
    handle.write(reportStr)
