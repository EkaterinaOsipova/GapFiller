#!/usr/bin/env python
#
# This script takes maf file produced by lastal (wrong maf format)
# and outputs in the correct maf format (ucsc) to stdout.
# Example lastal maf:
# s chr13:118008639-127786733 1486600 4665 + 9778094 AAAGGTGATGAT
# Example ucsc maf:
# s chr19  35156324 490 +  58617616 AAAGGTGATGAT


import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--maf', type=str, help='maf file produced by lastal run')
args = parser.parse_args()

with open(args.maf, 'r') as inf:
    for line in inf.readlines():
        if line.startswith('s '):
            if len(line.split(':')) != 2:
                # different type of line:
                # s chr1 179183808 1519 + 3521390 ATAGTCAGTCA
                newLine = ' '.join(line.split())
            else:
                chrName = line.split(':')[0].split()[-1]
                absStart = line.split()[1].split(':')[-1].split('-')[0]
                relStart = line.split()[2].rstrip()
                restOfLine = line.split()[3:]
                ucscMafLine = 's {0} {1} {2}'.format(chrName, int(absStart) + int(relStart), ' '.join(restOfLine))
                newLine = ucscMafLine
        else:
            newLine = line.rstrip()
        print(newLine)