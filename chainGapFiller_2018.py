#!/usr/bin/env python3

# Ekaterina Osipova, MPI-CBG/MPI-PKS, 2018

##########
# What is does. And how it conceptually works (e.g. input is a chain or chainIDs, for each chain, find a gap of a certain size, run a local lastz job, since the resulting alignments can overlap, chain them. Select the best 'mini-chain' and directly add this into the gap. Then continue iteratiing.
# Note: that 
##########

import sys
import os
import subprocess
import argparse
import logging
import re
import io
import time
import tempfile
############################


############################
# Make a list of jobs to run in temp shell script
# Returns nothing
#
# inChain: string containing all chains
# outf: path to output file; shell commands will be written into this file 
def MakeShellList(inChain, outf):

    # open input file or fail
    try:
        fhout = open(outf, 'w')
    except IOError as e:
        logging.error('Cannot write to shell script', e.errno, e.strerror)
        sys.exit(1)

    # write header for shell script
    fhout.write("#!/usr/bin/env bash\n")
    fhout.write("#set -o pipefail\n")
    fhout.write("#set -e\n")

    # count gaps patched in this file
    gap_count = 0

    # Nmbrs for tracking the line number
    lineNmbr = 0

    # incremeted identifier for mini chains
    #enumerator = 1

    # two bit files
    T2bit = args.T2bit
    Q2bit = args.Q2bit

    lastzVar = args.lastz
    axtChainVar = args.axtChain
    chainSortVar = args.chainSort

    # put chain string into a list
    ChainList = iter([i + '\n' for i in inChain.split('\n')])

    # change to acces by index
    for line in ChainList:

        lineNmbr += 1

        ll = line.split()
        if len(ll) > 0:
            if ll[0] == 'chain':
                # read the chain line:
                # e.g. chain 196228 chr4 62094675 + 12690854 12816143 chr23 24050845 - 20051667 20145391 1252
                score = int(ll[1])
                tName, tStart, tEnd = ll[2], int(ll[5]), int(ll[6])
                qName, qStartx, qEndx, qStrand = ll[7], int(ll[10]), int(ll[11]), ll[9]
                qSize = int(ll[8])
                logging.info('qStrand = {}'.format(qStrand))
                # changing coords for -strand if necessary

                qStart = qStartx
                qEnd = qEndx
                lastzParameters = args.lastzParameters + ' --strand=plus'

                if qStrand == "-":
                    lastzParameters = args.lastzParameters + ' --strand=minus'

                if ll[4] != "+":
                    logging.error("ERROR: target strand is not + for chain:{}".format(line))
                    sys.exit(1)
                # check if we consider this chain
                logging.info("score of this chain = {}".format(score))
                if ((score >= args.chainMinScore) and (tEnd - tStart >= args.chainMinSizeT) and \
                            (qEnd - qStart >= args.chainMinSizeQ)):
                    logging.info("valid chain")
                    curTPos = tStart
                    curQPos = qStart
                    #Chain_string += line

                    line = next(ChainList)
                    lineNmbr += 1

                    while(re.match("^\d+", line) is not None):
                    #while (line != '\n'):

                        a = line.split()
                        if (len(a) == 1):
                            TblockEnd = curTPos + int(a[0])
                            QblockEnd = curQPos + int(a[0])
                            logging.info("it was the last block\n")
                            #Chain_string += line
                            #Chain_string += '\n'
                        else:
                            blockLen = int(a[0])
                            TblockEnd = curTPos + blockLen
                            QblockEnd = curQPos + blockLen
                            TgapEnd = curTPos + blockLen + int(a[1])
                            QgapEnd = curQPos + blockLen + int(a[2])
                            tGapSpan = TgapEnd - TblockEnd
                            qGapSpan = QgapEnd - QblockEnd
                            # check if we want to patch this gap
                            if ((tGapSpan >= args.gapMinSizeT) and (tGapSpan <= args.gapMaxSizeT) and (
                                qGapSpan >= args.gapMinSizeQ) and (qGapSpan <= args.gapMaxSizeQ)):
                                logging.info("yes, this gap will be patched: {}".format(line.strip()))
                                TblockEnd += 1
                                QblockEnd += 1

                                # replace the content of the unmask by '[unmask]' if the user sets this flag, otherwise ''
                                if args.unmask:
                                    unmask = '[unmask]'
                                else:
                                    unmask = ''

                                if qStrand == "-":
                                    realQblockEnd = qSize - QgapEnd + 1
                                    realQgapEnd = qSize - QblockEnd + 1
                                else:
                                    realQblockEnd = QblockEnd
                                    realQgapEnd = QgapEnd

                                logging.info("running lastz on the block:")
                                regionToBePatched = [tName, str(TblockEnd), str(TgapEnd), qName, str(realQblockEnd), str(realQgapEnd)]
                                logging.info(' '.join(regionToBePatched))

                                # making lastz command for this region

                                #command_lastz = lastzVar + ' {0}/{1}[{2}..{3}]{9} {4}/{5}[{6}..{7}]{9} --format=axt {8} | ' + axtChainVar + ' -linearGap=loose stdin {0} {4} stdout 2> /dev/null | ' + chainSortVar + ' stdin stdout'.format(T2bit, tName, TblockEnd, TgapEnd, Q2bit, qName, realQblockEnd, realQgapEnd, lastzParameters, unmask)
                                command1 = ' {0}/{1}[{2}..{3}]{9} {4}/{5}[{6}..{7}]{9} --format=axt {8} | '.format(T2bit, tName, TblockEnd, TgapEnd, Q2bit, qName, realQblockEnd, realQgapEnd, lastzParameters, unmask)
                                command2 = ' -linearGap=loose stdin {0} {1} stdout 2> /dev/null | '.format(T2bit, Q2bit)
                                command3 = ' stdin stdout'
                                command_lastz = lastzVar + command1 + axtChainVar + command2 + chainSortVar + command3
                                ### adding this lastz run to a shell command list; lineNmbr - 1 because we start with 1 and later with 0
                                shellCommand = 'echo -e "LINE{0}\\n{1}\\n{2}\\n{3}\\n{4}\\n{5}\\n"; {6}; echo -e "LINE{0}\\n"\n'.\
                                    format((lineNmbr-1), blockLen, TblockEnd, TgapEnd, realQblockEnd, realQgapEnd, command_lastz)
                                fhout.write(shellCommand)
                                #logging.error("Edit line number:" + str(lineNmbr))

                            curQPos = QgapEnd
                            curTPos = TgapEnd


                        # get next line, break if no more line in string
                        try:
                             line = next(ChainList)
                             lineNmbr += 1
                        except Exception:
                             break
                else:
                    logging.info("invalid chain\n")

                    # save chain header line; get next line
                    line = next(ChainList)
                    lineNmbr += 1

                    # read the rest of the chain store blocks
                    while(re.match("^\d+", line) is not None):                        
                        try:
                            line = next(ChainList)
                            lineNmbr += 1
                        except Exception:
                            break

    logging.info("Done with reading gaps")
    logging.info("Gaps patched in this chain = {}".format(gap_count))
    logging.info('\n')
    logging.info('\n')
    #MergedChainString = MergeZeros(Chain_string)
    
    # close file handler
    fhout.close()

    return()

######################
# Takes temp file with all shell commands to run and returns lastz output in a single string

def RunAllShell(shellfile):

    AllShellCommand = 'bash ' + shellfile
    #print(AllShellCommand)

    try:
        allMiniChains = subprocess.check_output(AllShellCommand, shell=True)

        '''
        # for debugging write all mini chains to a file
        with open("allminifile", 'w') as f:
            allMiniChainsStr = allMiniChains.decode()
            for el in allMiniChainsStr.split('\n'):
                f.write(el)
                f.write('\n')
        '''

    except subprocess.CalledProcessError as shellrun:
        #logging.error('shell command failed', shellrun.returncode, shellrun.output)
        sys.exit(1)

    allMiniChains = allMiniChains.decode()

    return(allMiniChains)



######################
# Takes the whole lastz output chain list and return list containing chain block starting at  cur_position
# 
# (lineNmbr: [blockLen, TblockEnd, TgapEnd, realQblockEnd, realQgapEnd, 'allChainsStrings']) from LINE# to LINE#
# returns a dictionary
#
# allMiniChainsSplit: list containing line-wise lastz output file includung LINE### statements as block separator
# positions: starting position of current block (should start with LINE###
#
# returns list of line split strings: one output block starting with LINE### and ending with LINE#
# raises ValueError if block not properly separated by LINE#
#
def GetChainBlockFromLastzOutput(allMiniChainsSplit, cur_position):
    
    position = cur_position
    start = position
    end = None
    line = allMiniChainsSplit[position]
    reLine = re.compile(r'LINE\d+')

    # check whether initial line start with LINE#, raise error otherwise
    if reLine.match(line) is not None:
        
        position += 1 #get next line 
        line = allMiniChainsSplit[position]
        
        # process block until LINE# as end separator is encountered
        while reLine.match(line) is None:
              position += 1 #get next line 
              line = allMiniChainsSplit[position]

        # check that last line contains LINE#
        if reLine.match(line) is not None: 
            end = position
        else:
            raise ValueError('ERROR! allMiniChainsSplit end separator line at position ' + str(position) + ' does not start with LINE...')

    else:
        raise ValueError('ERROR! allMiniChainsSplit start separator line at position ' + str(position) + ' does not start with LINE...')

    curBlockList = allMiniChainsSplit[start:(end+1)]

    return(curBlockList)


######################
# Takes first chain from a chain list
# returns the header: "chain 52633 chr..." and a list of lines of this chain
# returns twice None if no chains are present
#
def TakeFirstChainFromList(chainList):
    
    headLine = None
    chainContent = None
    chainStart = None
    chainEnd = None

    for pos in range(0, len(chainList)):

        line = chainList[pos]

        # check if chain line
        m = re.match(r'chain', line)
        if m is not None:
            headLine = line.strip('\n')
            
            # process and store end position
            pos +=1
            line = chainList[pos]
            chainStart = pos
            
            while(re.match("^\d+", line) is not None):                
                pos +=1
                line = chainList[pos]            
            chainEnd = pos # actually position after chain            

            # don't process lower scoring chains
            break 

    # extract chain
    if chainStart is not None:
        chainContent = chainList[chainStart:chainEnd]

    return(headLine, chainContent)




######################
# Enumerates all mini chains and writes them to a file

def WriteMiniChainsFile(s, outfile, enum):
    list = [l + '\n' for l in s.split('\n') if l]
    ouf = open(outfile, 'a')
    for element in list:
        if element.startswith('chain'):
            element = ' '.join(element.split()[:-1]) + '\t{}\n'.format(enum)
            enum += 1
            ouf.write(element)
        else:
            ouf.write(element)
    ouf.close()
    return(enum)


###################
# After patching chain we need to insert it back on the right place
# InsertChainContent calculates new coordinates for a chain to be inserted
# and returns a list of lines, that were changed in comparision with an old chain file

def InsertChainContent(ChainContent, BestChain, BlockLenA, TblockEnd, TgapEnd, LoQblockEnd, LoQgapEnd):
    TlastzStart = int(BestChain.split()[5]) + 1
    TlastzEnd = int(BestChain.split()[6])

    if BestChain.split()[9] == '+':
        QlastzStart = int(BestChain.split()[10]) + 1
        QlastzEnd = int(BestChain.split()[11])
    else:
        # recalculate -strand coords to +strand:
        QlastzStart = int(BestChain.split()[8]) - int(BestChain.split()[10])
        QlastzEnd = int(BestChain.split()[8]) - int(BestChain.split()[11]) + 1

        tempQ = LoQgapEnd
        LoQgapEnd = LoQblockEnd
        LoQblockEnd = tempQ

    BlockToAdd = []

    if BestChain.split()[9] == '+':
        FirstQgap = abs(QlastzStart - int(LoQblockEnd))
        LastQGap = abs(int(LoQgapEnd) - QlastzEnd)
    else:
        FirstQgap = abs(QlastzStart - int(LoQblockEnd))
        LastQGap = abs(int(LoQgapEnd) - QlastzEnd)

    FirstLine = str(BlockLenA) + '\t' + str(TlastzStart - int(TblockEnd)) + '\t' + str(FirstQgap) + '\n'

    BlockToAdd.append(FirstLine)
    for i in range(0, len(ChainContent)-1):
        BlockToAdd.append(ChainContent[i])

    LastLine = ChainContent[len(ChainContent)-1].strip() + '\t' + str(int(TgapEnd)-TlastzEnd) +\
                   '\t' + str(LastQGap) + '\n'
    BlockToAdd.append(LastLine)
    return(BlockToAdd)



######## main ###########

# tracking runtime
start_time = time.time()

# initialize variables
chainIDs = ""

# use list containing batches of max. maxChainIDs chainIDs per chainExtract call
maxChainIDs = 5000
listChainID = []


###### Optional arguments #######
parser = argparse.ArgumentParser(description='This script extracts a chain from all.chain file by ID, finds gaps and using lastz patches these gaps, then inserts new blocks to a chain', \
                                 epilog='Example of use: ChainHelper.py -c hg38.speTri2.all.chain -ix hg38.speTri2.all.bb -T2 hg38.2bit -Q2 speTri2.2bit -um -m mini.458.165.chains -o out.458.165.chain')
parser.add_argument('-v', '--verbose', action='store_true', \
                    help="if -v is not specified, only ERROR messages will be shown")
#parser.add_argument('-m', '--mini', type=str, help="name of a file to put all mini chains from lastz")
parser.add_argument('--idList', type=str, help="idList=X,Y,Z a list of IDs of chains that have to be patched")
parser.add_argument('--idListFile', type=str, help="File containing list of IDs of chains that have to be patched, one per row")
parser.add_argument('--lastz', '-l', type=str, default='lastz', help="path to lastz executable, default = lastz")
parser.add_argument('--axtChain', '-x', type=str, default='axtChain', help="path to axtChain executable, default = axtChain")
parser.add_argument('--chainExtractID', '-cid', type=str, default='chainExtractID', help="path to chainExtractID executable, default = chainExtractID")
parser.add_argument('--chainSort', '-s', type=str, default='chainSort', help="path to chainSort executable, default = chainSort")
parser.add_argument('--output', '-o', type=str, help="name of output chain file. If not specified chains go to stdout")
parser.add_argument('--workdir', '-w', type=str, default='./', help="working directory for temp files, default = ./")

# initial parameters
parser.add_argument('--chainMinScore', '-mscore', type=int, default=0, help="consider only chains with a chainMinScore, default consider all")
parser.add_argument('--chainMinSizeT', '-mst', type=int, default=0, help="consider only chains with a chainMinSizeT, default consider all")
parser.add_argument('--chainMinSizeQ', '-msq', type=int, default=0, help="consider only chains with a chainMinSizeQ, default consider all")
parser.add_argument('--gapMinSizeT', '-gmint', type=int, default=10, help="patch only gaps that are at least that long on the target side, default gmint = 10")
parser.add_argument('--gapMinSizeQ', '-gminq', type=int, default=10, help="patch only gaps that are at least that long on the query side, default gminq = 10")
parser.add_argument('--gapMaxSizeT', '-gmaxt', type=int, default=100000, help="patch only gaps that are at most that long on the target side, default gmaxt = 100000")
parser.add_argument('--gapMaxSizeQ', '-gmaxq', type=int, default=100000, help="patch only gaps that are at most that long on the query side, default gmaxq = 100000")
parser.add_argument('--lastzParameters', '-lparam', type=str, default=' K=1500 L=2000 M=0 T=0 W=6 ', help="line with lastz parameters, default 'K=1500 L=2000 M=0 T=0 W=6' ")
parser.add_argument('--unmask', '-um', action='store_true', help="unmasking (lower case to upper case) characters from the 2bit files")
parser.add_argument('--scoreThreshold', '-st', type=int, default=2000, help="insert only chains that have at leaset this score, default st = 2000")
parser.add_argument("--index", "-ix", type=str, help="index.bb file for chains")

###### Required arguments #######
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument("--chain", "-c", type=str, help="all.chain file", required=True)
requiredNamed.add_argument("--T2bit", '-T2', type=str, help="reference 2bit file", required=True)
requiredNamed.add_argument("--Q2bit", '-Q2', type=str, help="query 2bit file", required=True)
args = parser.parse_args()



# check conflicting parameters
if args.idList and args.idListFile:
    logging.error("ERROR! Choose either idListFile or idList. Both is not supported.")
    sys.exit(1)

if args.idList or args.idListFile:
    if not args.index:
        logging.error("ERROR! index most be specified if idListFile or idList is used.")
        sys.exit(1)

if args.verbose:
    logging.basicConfig(level=logging.INFO)

if args.idList:
    chainIDs = args.idList

if args.idListFile:
    try:
        # read IDList file
        with open(args.idListFile, 'r') as fin:
            chainIDLines = fin.read().split('\n')
    except IOError as e:
        logging.error('Cannot read to idListFile', e.errno, e.strerror)
        sys.exit(1)

    # we might exceed shell argument maximum byte size, therefore we use smaller batches of chains
    # for chainExtractID; most simple way of getting batches of chains
    nchain = 0
    tchainstr = "" #add IDs to string until maxChainIDs chains reached
    for line in iter(chainIDLines):
        tchainstr = tchainstr + "," + line
        nchain += 1

        if nchain == maxChainIDs:
            listChainID.append(tchainstr)
            nchain = 0
            tchainstr = ""

    # get final batch of chainIDs
    listChainID.append(tchainstr)
                                    

# extract chains with these ids from all.chain file

# get chains; either iterate over batches of chains or extract single chain string
CurrentChainString = ""
if listChainID:

    # loop over chains and concatenate string
    for chainIDs in iter(listChainID):

        extract_command = '{0} {1} {2} stdout -idList={3}'.format(args.chainExtractID, args.index, args.chain, chainIDs) 

        try:
            TCurrentChainString = subprocess.check_output(extract_command, shell=True)
            TCurrentChainString = TCurrentChainString.decode()
            logging.info('Running patching on the chain ID = {}'.format(id))        
            CurrentChainString = CurrentChainString + "\n" + TCurrentChainString
        except subprocess.CalledProcessError as extractrun:
            logging.error('extract chain command failed', extractrun.returncode, extractrun.output)
            sys.exit(1)
else:
    
    if chainIDs:
        extract_command = '{0} {1} {2} stdout -idList={3}'.format(args.chainExtractID, args.index, args.chain, chainIDs)
        try:
            CurrentChainString = subprocess.check_output(extract_command, shell=True)
            CurrentChainString = CurrentChainString.decode()
            logging.info('Running patching on the chain ID = {}'.format(id))        
        except subprocess.CalledProcessError as extractrun:
            logging.error('extract chain command failed', extractrun.returncode, extractrun.output)
            sys.exit(1)    
    else:
        try:
            with open(args.chain, 'r') as content_file:
                CurrentChainString = content_file.read()
        except IOError as e:
            logging.error('Cannot read chain file {} {}'.format(e.errno, e.strerror) )
            sys.exit(1)

#if args.verbose:
#    logging.basicConfig(level=logging.INFO)


# 1) Loop through .all.chain and make a jobList in a shell script
if not os.path.isdir(args.workdir):
    logging.error("ERROR! Working directory '" + args.workdir + "' does not exist.")
    sys.exit(1)

# create 
try:
    temp = tempfile.NamedTemporaryFile(prefix = "tempCGFjobList", dir = args.workdir, delete = False )
    temp.close()
except PermissionError as e:
    logging.error('ERROR! Failed to create temporary file inside \'{}\'. {} {}'.format(args.workdir, e.errno, e.strerror) )
    sys.exit(1)

# Find gaps and write corresponding jobs to a shell script
MakeShellList(CurrentChainString, temp.name)


# 2) Run the prepared shell script
allMiniChains = RunAllShell(temp.name)

# remove temp file
os.unlink(temp.name)

# string will be used for final output
OutputChain = ''

# position of next block in miniChain output 
next_pos = 0
# next line number where chain is filled
next_lineNmbr = None

# regexp for getting lineNumber
relineNmbr = re.compile(r'LINE(\d+)')

# list of initial chains
currentChainLines = [i + '\n' for i in CurrentChainString.split('\n')]

# list of mini Chain Blocks
allMiniChainLines = [i + '\n' for i in allMiniChains.split('\n')]
lenAllMiniChainLines = len(allMiniChainLines)

# get current mini chain
curMiniBlockLines = GetChainBlockFromLastzOutput(allMiniChainLines, next_pos)
# get next lineNmbr where intial chain will be filled with gaps
m = relineNmbr.match(curMiniBlockLines[0])
if m is not None:
    next_lineNmbr = int(m.group(1))
else:
    raise ValueError('ERROR! Could not extract line number from separator current miniChain block')



### process initial chain and fill gaps with miniChainBlocks

# write output; open file handle
if args.output:
    try:
        ouf = open(args.output, 'w')
    except IOError as e:
        logging.error('Cannot write to outputfile', e.errno, e.strerror)
        sys.exit(1)
    
for lineNmbr in range(0, len(currentChainLines)):

    # get current initial chain line
    line = currentChainLines[lineNmbr]
    
    # empty output string
    OutputChain = ""

    # update chain
    if lineNmbr == next_lineNmbr:
        
        # strip first and last line containing LINE# from block
        valueList = curMiniBlockLines[1:(len(curMiniBlockLines)-1)] 
        Coords = valueList[:5]
        # remove new lines from Coords elements
        Coords = map(lambda s: s.strip(), Coords)

        # update next_pos and get next mini chain block; +1 since we have new line after each block in the output
        next_pos = next_pos + len(curMiniBlockLines) + 1 
        # test that we are not out of bounds, i.e. last entry, -1 since last line is new line
        if next_pos < lenAllMiniChainLines-1:
            curMiniBlockLines = GetChainBlockFromLastzOutput(allMiniChainLines, next_pos)            
            # get next lineNmbr
            m = relineNmbr.match(curMiniBlockLines[0])
            if m is not None:
                next_lineNmbr = int(m.group(1))
            else:
                raise ValueError('ERROR! Could not extract line number from separator current miniChain block')
        
        # get chain to be inserted
        BestChain, ChainContent = TakeFirstChainFromList(valueList[5:])                

        # insert nothing if no chain in block
        if BestChain is not None:                        
            if int(BestChain.split()[1]) >= args.scoreThreshold:
                logging.info('Best lastz output chain = {}'.format(BestChain))

                InsertBlock = InsertChainContent(ChainContent, BestChain, *Coords)
                OutputChain = ''.join(InsertBlock)
                logging.info('--- %s seconds ---' % (time.time() - start_time))
            else:
                logging.info("lastz output chains have low score\n")
                OutputChain = line                
        else:
            logging.info("lastz changed nothing in this block\n")
            OutputChain = line
    else:
        # Just add this line to the chain string and go further
        OutputChain = line        

    # print output
    if args.output:    
        ouf.write(OutputChain)
    else:
        print(OutputChain)

# close output file handle
if args.output:    
    ouf.close()

logging.info('--- Final runtime: %s seconds ---' % (time.time() - start_time))
