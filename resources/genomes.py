#! /usr/bin/env python

from collections import defaultdict

import matplotlib.cm as cm

class Genome(object):

    def labels(self):
        return self._labels

    def sizes(self):
        return self._sizes

    def size(self, key):
        return self._sizes[key]

    def targets(self):
        return self._targets

    def targetDictionary(self):
        pass

    def colors(self, N):
        colors = {}
        for i, n in enumerate( sorted(self._sizes.keys()) ):
            color = cm.gist_rainbow(i%(N+1) / float(N))
            colors[i] = color
            colors[n] = color
        return colors


class HG19(Genome):

    _labels = ["chr{0}".format(l) for l in range(1,23) + ['M', 'X', 'Y']]

    _sizes = {"chr1":  249250621, "chr2":  243199373,  "chr3": 198022430, "chr4":  191154276,
              "chr5":  180915260, "chr6":  171115067,  "chr7": 159138663, "chr8":  146364022,
              "chr9":  141213431, "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
              "chr13": 115169878, "chr14": 107349540, "chr15": 102531392, "chr16":  90354753,
              "chr17":  81195210, "chr18":  78077248, "chr19":  59128983, "chr20":  63025520,
              "chr21":  48129895, "chr22":  51304566,
              "chrM":      16571, "chrX":  155270560,  "chrY":  59373566}

    ## Locus,ChrName,ChrIdx,GeneStart,RegionStart,RegionEnd,GeneEnd
    _targets = [["HTT", "chr4", 4, 3075691, 3076603, 3076661, 3076815],
                ["FMR1", "chrX", 23, 146993123, 146993568, 146993629, 146994131],
                ["ALS", "chr9", 9, 27572985, 27573522, 27573541, 27574014],
                ["FUCHS", "chr18", 18, 53251995, 53253386, 53253458, 53253577],
                ["SCA10", "chr22", 22, 46190744, 46191234, 46191305, 46191756],
                ["EWINGS_Chr20", "chr20", 20, 21553989, 21556922, 21557001, 21557036],
                ["EWINGS_ChrX", "chrX", 23, 30325813, 30328875, 30328976, 30329062]]

    def __init__(self):
        # Update the size dictionary so we can index by index as well as name
        for i, c in enumerate(self._labels):
            self._sizes[i] = self._sizes[c]

    def targetDictionary(self):
        tDict = defaultdict(list)
        for t in self.targets():
            if t[2] <= 22:
                target_tId = t[2]-1
            else:
                target_tId = t[2]
            tDict[target_tId].append( t )
        return tDict

class GRC38(Genome):

    _labels = ["chr{0}".format(l) for l in range(1,23) + ['X', 'Y', 'M']]

    _sizes = {"chr1":  248956422, "chr2":  242193529,  "chr3": 198295559, "chr4":  190214555,
              "chr5":  181538259, "chr6":  170805979,  "chr7": 159345973, "chr8":  145138636,
              "chr9":  138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
              "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
              "chr17":  83257441, "chr18":  80373285, "chr19": 58617616 , "chr20": 64444167,
              "chr21":  46709983, "chr22":  50818468,
              "chrX":  156040895,  "chrY":  57227415, "chrM": 16569}

    ## Locus,ChrName,ChrIdx,GeneStart,RegionStart,RegionEnd,GeneEnd
    _targets = [["HTT", "chr4", 3, 3072621, 3074866, 3074949, 3075351],
                ["ALS", "chr9", 8, 27571412, 27573474, 27573556, 27574248],
                ["FUCHS", "chr18", 17, 55584764, 55586145, 55586237, 55586346],
                ["EWINGS_Chr20", "chr20", 19, 21573351, 21576271, 21576374, 21576399],
                ["SCA10", "chr22", 21, 45793649, 45795344, 45795434, 45796093],
                ["EWINGS_ChrX", "chrX", 22, 30307696, 30310741, 30310899, 30310946],
                ["FMR1", "chrX", 22, 147911603, 147912040, 147912120, 147914564]]

    def __init__(self):
        # Update the size dictionary so we can index by index as well as name
        for i, c in enumerate(self._labels):
            self._sizes[i] = self._sizes[c]

    def targetDictionary(self):
        tDict = defaultdict(list)
        for t in self.targets():
            tDict[t[2]].append( t )
        return tDict


def decodeGenome(genome):
    """
    """
    if genome.lower() == "hg19":
        return HG19()
    elif genome.lower() == "grc38":
        return GRC38()
    else:
        raise ValueError("Invalid genome: specified genome must be HG19 or GHC38")
