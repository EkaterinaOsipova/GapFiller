#!/usr/bin/env python
#
# This script takes two chains and a gene (coords) as input; extracts sequence of a query
# overlapping this gene, optionally excludes intronic regions;
# outputs percent of identity between sequences.

# ../calculateTwoChainsIdentity.py -c /projects/hillerlab/genome/gbdb-HL/galGal5/lastz/vs_HLcalAnn2B/axtChain/galGal5.HLcalAnn2B.allfilled.chain \
#    -i /projects/hillerlab/genome/gbdb-HL/galGal5/lastz/vs_HLcalAnn2B/axtChain/galGal5.HLcalAnn2B.allfilled.ID.bb \
#    -g ENSGALT00000088977.1 -b ensGene.galGal5.extracted.bed12 -t2 /projects/hillerlab/genome/gbdb-HL/galGal5/galGal5.2bit \
#      -q2 /projects/hillerlab/genome/gbdb-HL/HLcalAnn2B/HLcalAnn2B.2bit --id1 4230  --id2 10643 --introns

import argparse
import subprocess
import sys
from collections import defaultdict
#import Levenshtein
import difflib


def deltaDirac(x):
    return x if x > 0 else 0


def getChainInAxtSeq(id, chainFile, indexFile, ref2bit, query2bit):
    # extract chain by ID
    commandExtractChain = 'chainExtractID -idList={0} {1} {2} stdout'.format(id, indexFile, chainFile)
    chainOut = subprocess.Popen(commandExtractChain.split(), stdout=subprocess.PIPE)

    # convert chain to axt seq
    commandChainToAxt = 'chainToAxt stdin {0} {1} stdout'.format(ref2bit, query2bit)
    chainToAxtOut = subprocess.Popen(commandChainToAxt.split(), stdin=chainOut.stdout, stdout=subprocess.PIPE)
    axtSeq = chainToAxtOut.stdout.read()
    return axtSeq


def findGeneInBed12(gene, bed12File):
    geneLine = "ERROR: no transcript with ID {}".format(gene)
    with open(bed12File, 'r') as infBed:
        for line in infBed.readlines():
            if line.split()[3] == gene:
                geneLine = line
    return geneLine


def calcExonAbsCoords(geneLine):
    # Prepare everything to exclude exonic regions
    # parse gene line, make a list of abs exonic coords
    exonSizeList = map(int, geneLine.split()[10].rstrip(',').split(','))
    exonStartList = map(int, geneLine.split()[11].rstrip(',').split(','))
    exonAbsCoordsList = [(int(geneLine.split()[1]) + exonStartList[i], int(geneLine.split()[1]) + exonStartList[i] + exonSizeList[i]) \
                         for i in range(len(exonSizeList))]
    return exonAbsCoordsList


def extractSubSeq(axtSeq, geneLine, exonAbsCoordsList):
    # parse axt
    # make a dict like {axtHeader1: seqQuery, axtHeader2: seqQuery, ...}
    axtDict = {}
    # print(len(axtSeq1.split('\n\n')))

    for axtBlock in axtSeq.split('\n\n'):
        # add this block to the dict of it's not empty
        if axtBlock != '':
            axtDict[axtBlock.split('\n')[0]] = axtBlock.split('\n')[2]

    # Make an ordered list of dict keys
    # Process this list and exclude blocks that do not overlap a gene under consideration
    overlappingAxtBlocks = []
    geneChrom = geneLine.split()[0]
    geneStart = int(geneLine.split()[1])
    geneEnd = int(geneLine.split()[2])

    for axtBlockHeader in axtDict.keys():
        blockChrom = axtBlockHeader.split()[1]
        blockStart = int(axtBlockHeader.split()[2])
        blockEnd = int(axtBlockHeader.split()[3])

        if (blockChrom == geneChrom) and (((blockStart >= geneStart) and (blockStart <= geneEnd)) or \
                                                  ((blockEnd >= geneStart) and (blockEnd <= geneEnd)) or \
                                                  ((blockEnd >= geneEnd) and (blockStart <= geneStart)) or \
                                                  ((blockStart >= geneStart) and (blockEnd <= geneEnd))):
            print("found an overlapping block:")
            overlappingAxtBlocks.append(axtBlockHeader)
            print(axtBlockHeader)

    if len(overlappingAxtBlocks) == 0:
        print('ERROR: found no overlapping axt blocks!')
        sys.exit(1)

    # sort blocks by order on a chrom
    overlappingAxtBlocks.sort(key=lambda x: x.split()[2])

    # calculate delta for the first and the last blocks
    cutFromFirstBlock = geneStart - int(overlappingAxtBlocks[0].split()[2])
    cutFromLastBlock = int(overlappingAxtBlocks[-1].split()[3]) - geneEnd


    # make a dict like: {axtBlockHeader1: [(), (), ()], ...}
    # where each axtBlock is a key and value is a list of coords (start, end) of exons overlapping this block
    exonsOverlapDict = defaultdict(list)
    for axtBlockHeader in overlappingAxtBlocks:
        blockStart = int(axtBlockHeader.split()[2])
        blockEnd = int(axtBlockHeader.split()[3])

        for exon in exonAbsCoordsList:
            exonStart = exon[0]
            exonEnd = exon[1]
            # if this exon overlaps this block, add to dict tuple with relative exonic coords
            if ((blockStart >= exonStart) and (blockStart <= exonEnd)) or \
            ((blockEnd >= exonStart) and (blockEnd <= exonEnd)) or \
            ((blockEnd >= exonEnd) and (blockStart <= exonStart)) or \
                    ((blockStart >= exonStart) and (blockEnd <= exonEnd)):
                exonsOverlapDict[axtBlockHeader].append((exonStart - blockStart, exonEnd - blockStart))


    # extract sequence for each block. For the first and the last do it manually
    # make a string with a full query sequence
    querySeq = ''

    if args.introns:
        # add the middle: if we exclude exons
        for axtBlockHeader in overlappingAxtBlocks:
            indexExon = [0]
            for exon in exonsOverlapDict[axtBlockHeader]:
                # make a list of indices
                indexExon.append(exon[0])
                indexExon.append(exon[1])
            indexExon.append(len(exonsOverlapDict[axtBlockHeader]))
            for i in range(0, len(indexExon) - 1, 2):
                querySeq += axtDict[axtBlockHeader][indexExon[i]:indexExon[i + 1]]
    else:
        # add the middle: if we don't exclude exons
        for axtBlockHeader in overlappingAxtBlocks:
            querySeq += axtDict[axtBlockHeader]


    # cut the first and the last block if they span more than a gene
    querySeq = querySeq[0 + deltaDirac(cutFromFirstBlock):]
    querySeq = querySeq[:len(querySeq) - deltaDirac(cutFromLastBlock)]

    return querySeq


def calculateIdentity(seq1, seq2):
    seq1 = seq1.upper().replace('-', '')
    seq2 = seq2.upper().replace('-', '')
    matcherSeq1Seq2 = difflib.SequenceMatcher(None, seq1, seq2)
    if len(seq2) < len(seq1):
        return sum(n for i,j,n in matcherSeq1Seq2.get_matching_blocks()) / float(len(seq2))
    else:
        return sum(n for i, j, n in matcherSeq1Seq2.get_matching_blocks()) / float(len(seq1))
    #return Levenshtein.ratio(seq1.upper().replace('-', ''), seq2.upper().replace('-', ''))


###############
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--chains', type=str, help='file with all chains')
parser.add_argument('-i', '--index', type=str, help='index bb file for this chain file')
parser.add_argument('-g', '--geneID', type=str, help='ENS transcript ID')
parser.add_argument('-b', '--bed12', type=str, help='ENS transcripts annotation file in bed12 format')
parser.add_argument('-t2', '--ref2bit', type=str, help='reference 2bit file')
parser.add_argument('-q2', '--query2bit', type=str, help='query 2bit file')
parser.add_argument('--id1', type=str, help='ID of the first chain to consider')
parser.add_argument('--id2', type=str, help='ID of the second chain to consider')
parser.add_argument('--introns', action='store_true', help='specify if you want to consider intronic regions only')

args = parser.parse_args()


###############
# get arguments
chainFile = args.chains
indexFile = args.index
bed12File = args.bed12
geneID = args.geneID
id1 = args.id1
id2 = args.id2
ref2bit = args.ref2bit
query2bit = args.query2bit

# get corresponding transcript line from bed12 file
geneLine = findGeneInBed12(geneID, bed12File)

print('\n....Looking for a gene entry....')
print(geneLine)

# get axt seq for each chain as a string
axtSeq1 = getChainInAxtSeq(id1, chainFile, indexFile, ref2bit, query2bit)
axtSeq2 = getChainInAxtSeq(id2, chainFile, indexFile, ref2bit, query2bit)

exonAbsCoordsList = calcExonAbsCoords(geneLine)
print('...Parsing chain ID {}...'.format(id1))
axtQuerySubSeq1 = extractSubSeq(axtSeq1, geneLine, exonAbsCoordsList)
print('\n...Parsing chain ID {}...'.format(id2))
axtQuerySubSeq2 = extractSubSeq(axtSeq2, geneLine, exonAbsCoordsList)
print('\nsequence query for chain ID {0}: \n{1}'.format(id1, axtQuerySubSeq1))
print('\nsequence query for chain ID {0}: \n{1}'.format(id2, axtQuerySubSeq2))

# calculate % identity
identity = calculateIdentity(axtQuerySubSeq1, axtQuerySubSeq2)
if args.introns:
    print('\nidentity between INTRONIC regions of the two blocks found: {}'.format(identity))
else:
    print('\nidentity between sequence of the two blocks found: {}'.format(identity))


