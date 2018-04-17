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
from collections import defaultdict

import matplotlib; matplotlib.use('agg')
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

plt.style.use('ggplot')

COLORS = sns.color_palette()

outputPrefix = sys.argv[1]
fns          = sys.argv[2:]

def RepeatCount( tpl ):
    totalCount = 0
    for elem in tpl.split('_'):
        seq, count = elem.split('x')[:2]
        totalCount += int(count)
    return totalCount

def RepeatCounts(df):
    return [RepeatCount(v) for v in df.columns.values[1:]]

def MaxRepeatCount( tpls ):
    maxV = 0
    for k, tplList in tpls.iteritems():
        maxV = max(maxV, max(RepeatCount(tpl) for tpl in tplList))
    maxV = ((maxV / 10) + 2) * 10
    return maxV  

def RepeatSize( tpl ):
    size = 0
    for elem in tpl.split('_'):
        seq, count = elem.split('x')[:2]
        size += len(seq) * int(count)
    return size

def RepeatSizes( columns ):
    return [RepeatSize(tpl) for tpl in columns]

def MaxRepeatSize( tpls ):
    maxV = 0
    for k, tplList in tpls.iteritems():
        maxV = max(maxV, max(RepeatSize(tpl) for tpl in tplList))
    maxV = ((maxV / 10) + 6) * 10
    return maxV 

def LogAddExp( row ):
    return reduce(np.logaddexp, [l for l in row[1:] if np.isfinite(l)])

def LogToProb( row ):
    sumLL = row[-1]
    return pd.Series([np.exp(x - sumLL) for x in row[1:-1]])

def CsvsToDataFrame( csvs ):
    data = None
    for fn in csvs:
        barcode = fn.split('.')[-3].split('-')[0]
        barcode = ''.join(c for c in barcode if c.isalnum())
        df = pd.read_csv(fn)
        counts = RepeatCounts(df)
        df["SumLL"] = df.apply(LogAddExp, axis=1)
        df = pd.DataFrame(df.apply(LogToProb, axis=1))
        colSums = np.array(df.sum(0))
        vals = colSums / sum(colSums)

        raw = {"variable": pd.Series(counts),
               "values": pd.Series(vals),
               "Barcode": barcode}
        if data is None:
            data = pd.DataFrame(raw)
        else:
            data = data.append(pd.DataFrame(raw))
    return data

def CsvsToHistogramDataFrame( csvs ):
    data = None
    for fn in csvs:
        barcode = fn.split('.')[-3].split('-')[0]
        barcode = ''.join(c for c in barcode if c.isalnum())
        df = pd.read_csv(fn)
        bestTpls = []
        for rowIdx, row in df.iterrows():
            maxCol = None
            maxVal = None
            for col, val in row.iteritems():
                if col == "ZmwId":
                    continue
                if maxVal is None or val > maxVal:
                    maxCol = col
                    maxVal = val
            if np.isfinite(maxVal):
                bestTpls.append( maxCol )
        bestSizes = RepeatSizes(bestTpls)

        raw = {"Sizes": pd.Series(bestSizes),
               "Barcode": barcode}
        if data is None:
            data = pd.DataFrame(raw)
        else:
            data = data.append(pd.DataFrame(raw))
    return data

def CsvsToCounts( csvs ):
    data = {}
    for fn in csvs:
        barcode = fn.split('.')[-3].split('-')[0]
        barcode = ''.join(c for c in barcode if c.isalnum())
        df = pd.read_csv(fn)
        counts = RepeatCounts(df)
        data[barcode] = []
        for index, row in df.iterrows():
            row    = list(row)[1:]
            maxLL  = max(row)
            maxIdx = row.index(maxLL)
            count  = counts[maxIdx]
            data[barcode].append( count )
        data[barcode] = sorted(data[barcode], reverse=True)
    return data

def PlotRepeatPMF( outputPrefix, name, tpls, csvs ):
    plotName = "{0}_{1}_zmws.png".format(outputPrefix.lower(), name.lower())
    maxRpt = MaxRepeatCount( tpls )
    data = CsvsToDataFrame( csvs )
    g = sns.FacetGrid(data, row="Barcode", size=2.0, aspect=6, xlim=(0,maxRpt))
    g = g.map(plt.plot, "variable", "values", color="darkblue")
    for i, (bc, tplList) in enumerate(sorted(tpls.iteritems())):
        ax = g.facet_axis(i, 0)
        ylim = ax.get_ylim()
        for tpl in tplList:
            ct = RepeatCount(tpl)
            ax.axvline(x=ct, color="red", ls='--')
            ax.text(ct+0.4, ylim[1]-0.07, ct, fontsize=15)
    g.set_axis_labels("Repeat Count", "Density")
    plt.savefig(plotName)
    plt.close()

    p = {"caption": "Repeat Analysis Probability Mass Plot",
           "image": plotName,
            "tags": [],
              "id": "{0} - Probability Mass Distribution for {1}".format(outputPrefix, name),
           "title": "{0} - Probability Mass Distribution for {1}".format(outputPrefix, name)}
    return p

def PlotRepeatHistogram( outputPrefix, name, tpls, csvs ):
    plotName = "{0}_{1}_histogram.png".format(outputPrefix.lower(), name.lower())
    maxSize = MaxRepeatSize( tpls )
    data = CsvsToHistogramDataFrame( csvs )
    g = sns.FacetGrid(data, row="Barcode", size=2.0, aspect=6, xlim=(0,maxSize))
    g = g.map(plt.hist, "Sizes", density=True, bins=list(range(0,maxSize+1,3)), color="darkblue")
    g.set_axis_labels("Repeat Region Sizes", "Density")
    plt.savefig(plotName)
    plt.close()
    p = {"caption": "Repeat Analysis Region Size Histogram",
           "image": plotName,
            "tags": [],
              "id": "{0} - Repeat Region Size Histogram for {1}".format(outputPrefix, name),
           "title": "{0} - Repeat Region Size Histogram for {1}".format(outputPrefix, name)}
    return p

def PlotRepeatDotPlot( outputPrefix, name, tpls, csvs ):
    plots = []

    plotNameRoot = "{0}_{1}_repeats".format(outputPrefix.lower(), name.lower())
    data = CsvsToCounts( csvs )
    for barcode, bestTpls in data.iteritems():

        x, y = [], []
        for i, tplSize in enumerate(bestTpls):
            for j in range(tplSize):
                x.append( j )
                y.append( i + 1 )

        plt.scatter(x=x, y=y, color=COLORS[1])
        plt.title("Per-ZMW Repeat Estimates for {0}\nBarcode {1}--{1}".format(name, barcode))
        plt.xlabel("Number of Repeat Units")
        plt.ylabel("ZMWs sorted by Repeat Size")

        # Add vertical lines representing each consensus
        called = sorted([int(tpl.split('x')[-1]) for tpl in tpls[barcode]])
        for call in called:
            plt.axvline(x=call, color="red", ls="--")
            plt.text(call+0.4, max(y)*0.93, call, fontsize=15, color="black")

        plotName = "{0}.{1}--{1}.png".format(plotNameRoot, barcode)
        plt.savefig(plotName)
        plt.close()

        plots.append( {"caption": "Repeat Analysis Per-ZMW Estimated Region Sizes",
                       "image": plotName,
                       "tags": [],
                       "id": "{0} - Repeat Region Size Histogram for {1} ({2}--{2})".format(outputPrefix, name, barcode),
                       "title": "{0} - Repeat Region Size Histogram for {1} ({2}--{2})".format(outputPrefix, name, barcode)}
                    )
    return plots

def SortFiles( fns ):
    httSeqs = defaultdict(list)
    httCsvs = []
    fmrSeqs = defaultdict(list)
    fmrCsvs = []

    for fn in fns:
        if fn.endswith('.fastq'):
            for line in open(fn):
                if line.startswith('@Barcode'):
                    bc = line.strip().split('_')[0][8:].split('--')[0]
                    bc = ''.join(c for c in bc if c.isalnum())
                    tpl = line.strip().split()[-1]
                    if tpl.startswith('CAGx') or tpl.startswith('CTGx'):
                        httSeqs[bc].append( tpl )
                    if tpl.startswith('CGGx') or tpl.startswith('GCCx'):
                        fmrSeqs[bc].append( tpl )
        elif "_zmws" in fn and fn.endswith('HTT.csv'):
            httCsvs.append(fn)
        elif "_zmws" in fn and fn.endswith('FMR1.csv'):
            fmrCsvs.append(fn)

    return httSeqs, httCsvs, fmrSeqs, fmrCsvs

def WriteReportJson( plotList=[], tableList=[] ):
    if plotList or tableList:
        reportDict = {"plots":plotList, "tables":tableList}
        reportStr = json.dumps(reportDict, indent=1)
        with open("report.json", 'w') as handle:
            handle.write(reportStr)
    else:
        with open("report.json", 'w') as handle:
            handle.write("{}")

httSeqs, httCsvs, fmrSeqs, fmrCsvs = SortFiles( fns )

plotList = []
if len(httSeqs) > 0 and len(httCsvs) > 0:
    p1 = PlotRepeatPMF(outputPrefix, "HTT", httSeqs, httCsvs)
    p2 = PlotRepeatHistogram(outputPrefix, "HTT", httSeqs, httCsvs)
    p3 = PlotRepeatDotPlot(outputPrefix, "HTT", httSeqs, httCsvs)
    plotList.append( p1 )
    plotList.append( p2 )
    plotList += p3
elif len(httSeqs) > 0 or len(httCsvs) > 0:
    raise Warning("Input Error! Recieved HTT sequences but there are no ZMW scores, skipping...")

if len(fmrSeqs) > 0 and len(fmrCsvs) > 0:
    p1 = PlotRepeatPMF(outputPrefix, "FMR1", fmrSeqs, fmrCsvs)
    p2 = PlotRepeatHistogram(outputPrefix, "FMR1", fmrSeqs, fmrCsvs)
    p3 = PlotRepeatDotPlot(outputPrefix, "FMR1", fmrSeqs, fmrCsvs)
    plotList.append( p1 )
    plotList.append( p2 )
    plotList += p3
elif len(fmrSeqs) > 0 or len(fmrCsvs) > 0:
    raise Warning("Input Error! Recieved FMR1 sequences but there are ZMW scores, skipping...")

WriteReportJson( plotList )
