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
from collections import defaultdict

from pbcore.io import FastqReader, openDataSet

MINIMUM_MAP_PROB = 0.99

dsFile = sys.argv[1]
fns = sys.argv[2:]

def SeparateFiles( fns ):
    seqFiles = []
    mapFiles = []
    for fn in fns:
        if fn.endswith(".fastq"):
            seqFiles.append( fn )
        elif "_mapping" in fn:
            mapFiles.append( fn )
        else:
            print "WARNING: Invalid data file '{0}', only FASTQ and MAPPING files accepted".format(fn)
    return seqFiles, mapFiles

def ReadSeqFiles( fns ):
    seqs = {}
    for fn in fns:
        for record in FastqReader( fn ):
            if record.id in seqs:
                print "ERROR: Duplicate sequence id '{0}'".format(record.id)
                raise SystemExit
            seqs[record.id] = record
    return seqs

def ReadMappingFiles( fns, seqs ):
    mappings = defaultdict(list)
    for fn in fns:
        currMap = defaultdict(list)
        with open(fn) as handle:
            # Parse the current file's header for Column/Consensus mapping
            header = {idx:seqId.strip() for idx, seqId in enumerate(handle.next().split(','))}
            for seqId in header.itervalues():
                if seqId != "SubreadId" and seqId not in seqs:
                    print "ERROR: No matching consensus found for mapped sequence '{0}'".format(seqId)
                    raise SystemExit
            # Parse each row, recording the position of the largest valid mapping probability
            for line in handle:
                parts = line.strip().split(',')
                subreadId = parts[0]
                for idx, part in enumerate(parts):
                    try:
                        prob = float(part)
                    except:
                        continue
                    if prob > MINIMUM_MAP_PROB:
                        currMap[idx].append( subreadId )
                        break
            # Add the current map to our totals while converting Idx->Id
            for idx, subreads in currMap.iteritems():
                seqId = header[idx]
                mappings[seqId] = mappings[seqId] + subreads
    return mappings

def GetPrefix( fn ):
    if fn.lower().endswith('set.xml') or fn.lower().endswith('subreads.bam'):
        return '.'.join(fn.split('.')[:-2])
    else:
        return '.'.join(fn.split('.')[:-1])

def WriteWhitelistedDataSets( dsFile, mappings ):
    prefix = GetPrefix( dsFile )
    for seqId, subreads in mappings.iteritems():
        sset = openDataSet( dsFile )
        sset.filters.addRequirement(qname=[('=', subreadId) for subreadId in subreads])
        #sset.filters.addRequirement(qname=[('=', sorted(subreads))])
        sset.write( prefix + "." + seqId + ".subreadset.xml" )

seqFiles, mapFiles = SeparateFiles( fns )
seqs = ReadSeqFiles( seqFiles )
mappings = ReadMappingFiles( mapFiles, seqs )
WriteWhitelistedDataSets( dsFile, mappings )
