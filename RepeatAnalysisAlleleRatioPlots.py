#! /usr/bin/env python

import sys
import argparse

from collections import Counter, defaultdict

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

NEG_INF = float('-inf')

def Mode( vals ):
    vals = Counter(vals.itervalues())
    return vals.most_common(1)[0][0]

def Median( vals ):
    vals = sorted(vals.itervalues())
    n = len(vals)
    if n % 2 == 1:
        return vals[(n+1)/2]
    if n % 2 == 0:
        low = n / 2
        high = low + 1
        return (vals[low] + vals[high]) / 2.0

def ParseSplitOption( string ):
    elem = [e.split(':') for e in string.split(';')]
    asDict = {e[0]:float(e[1]) for e in elem}
    return asDict

def CleanString( string ):
    return ''.join(c for c in string if (c.isalnum() or c in "_-/."))

def ParseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--prefix", help="Prefix to use for output files", default="output")
    parser.add_argument("-M", "--max", help="Set the maximum value repeat-size to plot", type=int)
    parser.add_argument("-s", "--subreads", nargs='*', help="Subread files")
    parser.add_argument("-m", "--midpoints", help="Size mid-points to split alleles at by barcode, format: '0--0:21;1--1:40;..'")
    parser.add_argument("files", nargs='*', help="input files", metavar="FILE")
    args = parser.parse_args()
    args.midpoints = ParseSplitOption(args.midpoints)
    args.subreads = ReadSubreadCounts( args.subreads ) if args.subreads else None
    return args

def ReadRepeatCounts( fn ):
    bestRpts = {}
    bestScores = {}
    with open(fn) as handle:
        repeats = [int(p.split('x')[-1]) for p in handle.next().split(',')[1:]]
        for line in handle:
            parts = line.strip().split(',')
            zmw = parts[0]
            scores = [float(p) for p in parts[1:]]
            maxLL  = max(scores)
            if maxLL == NEG_INF:
                continue
            maxIdx = scores.index(maxLL)
            maxRpt = repeats[maxIdx]

            bestScores[zmw] = maxLL
            bestRpts[zmw] = maxRpt
    return bestScores, bestRpts

def ReadSubreadCounts( files ):
    data = defaultdict(int)
    for fn in files:
        with open(fn) as handle:
            handle.next()
            for line in handle:
                parts = line.strip().split(',')
                subread = parts[0]
                zmw = '/'.join(parts[0].split('/')[:2])
                scores = [float(p) for p in parts[1:]]
                maxLL  = max(scores)
                if maxLL == NEG_INF:
                    continue
                data[zmw] += 1
    return data

def SummarizeData(files, splits, subreads):
    data = []
    for fn in files:
        cleanFn = CleanString(fn)
        barcode = cleanFn.split('.')[-3]
        bcPart = barcode.split('-')[0]
        if barcode in splits.keys():
            splitVal = splits[barcode]
        elif bcPart in splits.keys():
            splitVal = splits[bcPart]
        else:
            continue
        scores, rpts = ReadRepeatCounts( fn )
        large = {k:v for k,v in rpts.iteritems() if v > splitVal}
        largeMed = Median(large)
        largeMode = Mode(large)
        largeSub = sum(subreads[k] for k in large.keys()) if subreads else "N/A"
        small = {k:v for k,v in rpts.iteritems() if v < splitVal}
        smallMed = Median(small)
        smallMode = Mode(small)
        smallSub = sum(subreads[k] for k in small.keys()) if subreads else "N/A"
        nAssigned = len(large) + len(small)
        frac = round(100 * len(large) / float(nAssigned), 2)
        largeSubFrac = round(100 * largeSub / float(largeSub + smallSub), 2) if largeSub else "N/A"
        largeSubPerZ = round(largeSub / float(len(large)), 2) if largeSub else "N/A"
        smallSubFrac = round(100 * smallSub / float(largeSub + smallSub), 2) if smallSub else "N/A"
        smallSubPerZ = round(smallSub / float(len(small)), 2) if smallSub else "N/A"
        data.append([cleanFn, barcode, "Small", smallMed, smallMode, len(scores), nAssigned, len(small), 100.0-frac, smallSub, smallSubFrac, smallSubPerZ])
        data.append([cleanFn, barcode, "Large", largeMed, largeMode, len(scores), nAssigned, len(large), frac, largeSub, largeSubFrac, largeSubPerZ])
    return data

def WriteCsv(prefix, data):
    with open(prefix + "_summary.csv", 'w') as handle:
        handle.write("File,Barcode,Size,Median,Mode,ZmwCount,AssignedCount,AlleleCount,AlleleFrac,SubreadCount,SubreadFrac,AvgPerZmw\n")
        for row in data:
            row = [str(r) for r in row]
            handle.write("{}\n".format(",".join(row)))

def DataToDf( data ):
    raw = {"Barcode": pd.Series(row[1] for row in data),
           "AlleleClass": pd.Series(row[2] for row in data),
           "Median": pd.Series(row[3] for row in data),
           "Mode": pd.Series(row[4] for row in data),
           "AlleleFraction": pd.Series(row[-4] for row in data),
           "SubreadFraction": pd.Series(row[-2] for row in data),
           "AvgPerZmw": pd.Series(row[-1] for row in data)}
    return pd.DataFrame(raw)

def PlotData(prefix, maxV, doSubreads, df):
    g = sns.lmplot(x="Median", y="AlleleFraction", hue="AlleleClass",
                           truncate=True, size=5, data=df)
    g.set_axis_labels("Estimated Repeat Size (Median)", "Allelic Fraction")
    g.set(xlim=(0, maxV), ylim=(0, 100))
    plt.savefig(prefix + "_frac.png")
    plt.close()

    g = sns.lmplot(x="Median", y="AlleleFraction", truncate=True, size=5, data=df)
    g.set_axis_labels("Estimated Repeat Size (Median)", "Allelic Fraction")
    g.set(xlim=(0, maxV), ylim=(0, 100))
    plt.savefig(prefix + "_all_frac.png")
    plt.close()

    if doSubreads:
        g = sns.lmplot(x="Median", y="SubreadFraction", hue="AlleleClass",
                               truncate=True, size=5, data=df)
        g.set_axis_labels("Estimated Repeat Size (Median)", "Subread Fraction")
        g.set(xlim=(0, maxV), ylim=(0, 100))
        plt.savefig(prefix + "_sub_frac.png")
        plt.close()

        g = sns.lmplot(x="Median", y="SubreadFraction", truncate=True, size=5, data=df)
        g.set_axis_labels("Estimated Repeat Size (Median)", "Subread Fraction")
        g.set(xlim=(0, maxV), ylim=(0, 100))
        plt.savefig(prefix + "_all_sub_frac.png")
        plt.close()

        g = sns.lmplot(x="Median", y="AvgPerZmw", truncate=True, size=5, data=df)
        g.set_axis_labels("Estimated Repeat Size (Median)", "Average Subreads Per Zmw")
        g.set(xlim=(0, maxV))
        plt.savefig(prefix + "_per_zmw.png")
        plt.close()

def Run():
    args = ParseArguments()
    data = SummarizeData(args.files, args.midpoints, args.subreads)
    WriteCsv(args.prefix, data)
    df = DataToDf( data )
    doSubreads = False if args.subreads is None else True
    PlotData(args.prefix, args.max, doSubreads, df)

if __name__ == "__main__":
    Run()
