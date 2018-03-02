#! /usr/bin/env python

## Copyright (c) 2017, Pacific Biosciences of California, Inc.
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

outputPrefix = sys.argv[1]
fns          = sys.argv[2:]

def RepeatCounts(df):
    return [int(v.split('x')[-1]) for v in df.columns.values[1:]]

def LogAddExp( row ):
    return reduce(np.logaddexp, [l for l in row[1:] if np.isfinite(l)])

def LogToProb( row ):
    sumLL = row[-1]
    return pd.Series([np.exp(x - sumLL) for x in row[1:-1]])

def MaxRepeatCount( tpls ):
    maxV = 0
    for k, v in tpls.iteritems():
        maxV = max(maxV, max(v))
    maxV = ((maxV / 10) + 2) * 10
    return maxV

def CsvsToDataFrame( csvs ):
    data = None
    for fn in csvs:
        barcode = fn.split('.')[-3].split('-')[0]
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

def PlotRepeats( outputPrefix, name, csvs, tpls ):
    data = CsvsToDataFrame( csvs )
    maxRpt = MaxRepeatCount( tpls )
    plotName = "{0}_{1}_zmws.png".format(outputPrefix.lower(), name.lower())
    g = sns.FacetGrid(data, row="Barcode", size=2.0, aspect=6, xlim=(0,maxRpt))
    g = g.map(plt.plot, "variable", "values", color="darkblue")
    for i, (bc, counts) in enumerate(sorted(tpls.iteritems())):
        ax = g.facet_axis(i, 0)
        for ct in counts:
            ax.axvline(x=ct, color="red", ls='--')
    g.set_axis_labels("Repeat Count", "Density")
    plt.savefig(plotName)

    p = {"caption": "Repeat Analysis Probability Mass Plot",
           "image": plotName,
            "tags": [],
              "id": "{0} - Probability Mass Distribution for {1}".format(outputPrefix, name),
           "title": "{0} - Probability Mass Distribution for {1}".format(outputPrefix, name)}
    return p

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
                    tpl = line.strip().split()[-1]
                    ct = int(tpl.split('x')[-1])
                    if tpl.startswith('CAGx') or tpl.startswith('CTGx'):
                        httSeqs[bc].append( ct )
                    if tpl.startswith('CGGx') or tpl.startswith('GCCx'):
                        fmrSeqs[bc].append( ct )
        elif "_zmws" in fn and fn.endswith('HTT.csv'):
            httCsvs.append(fn)
        elif "_zmws" in fn and fn.endswith('FMR1.csv'):
            fmrCsvs.append(fn)
        else:
            raise SystemExit("Invalid input file! '{0}' is not a FASTQ or ZMW_CSV".format(fn))

    return httSeqs, httCsvs, fmrSeqs, fmrCsvs

def WriteReportJson( plotList=[], tableList=[] ):
    reportDict = {"plots":plotList, "tables":tableList}
    reportStr = json.dumps(reportDict, indent=1)
    with open("report.json", 'w') as handle:
        handle.write(reportStr)


httSeqs, httCsvs, fmrSeqs, fmrCsvs = SortFiles( fns )

plotList = []
if len(httSeqs) > 0 and len(httCsvs) > 0:
    p = PlotRepeats(outputPrefix, "HTT", httCsvs, httSeqs)
    plotList.append( p )
elif len(httSeqs) > 0 or len(httCsvs) > 0:
    raise Warning("Input Error! Recieved HTT sequences but is no ZMW scores, skipping...")

if len(fmrSeqs) > 0 and len(fmrCsvs) > 0:
    p = PlotRepeats(outputPrefix, "FMR1", fmrCsvs, fmrSeqs)
    plotList.append( p )
elif len(fmrSeqs) > 0 or len(fmrCsvs) > 0:
    raise Warning("Input Error! Recieved FMR1 sequences but is no ZMW scores, skipping...")

WriteReportJson( plotList )
