#!/sw/bin/python3
#

""" Chain_lib-based script. Merges chains with the same ID to one chain with this ID """
import argparse
from chain_lib import Chain
#from collections import Counter
from collections import defaultdict
import subprocess
import tempfile
import sys
import os


def main():
    """Parse args, split the chain file."""
    parser = argparse.ArgumentParser()
    requiredNamed = parser.add_argument_group('required named arguments')

    # Required arguments
    requiredNamed.add_argument("--chainFile", '-c', type=str, help="chain file to extract", required=True)
    requiredNamed.add_argument('--output', '-o', type=str, help="name of output chain file", required=True)
    requiredNamed.add_argument("--T2bit", '-T2', type=str, help="reference 2bit file", required=True)
    requiredNamed.add_argument("--Q2bit", '-Q2', type=str, help="query 2bit file", required=True)
    #requiredNamed.add_argument("--Tchroms", '-ts', type=str, help="reference chrom.sizes file", required=True)
    #requiredNamed.add_argument("--Qchroms", '-qs', type=str, help="query chrom.sizes file", required=True)
    requiredNamed.add_argument("--linearGapVar", '-lg', type=str, help="linearGap file", required=True)
    requiredNamed.add_argument('--tempStore', '-t', type=str, help="fake temporary file to store chains with same IDs", required=True)
    args = parser.parse_args()

    # Get a list of chain-related strings
    allChains = readChainFile(args.chainFile)

    # Get list of non unique IDs
    nonUniqIDs, uniqIDs, dictChains = getListNonUniqIDs(allChains)
    # print("number of non uniq = ", len(nonUniqIDs))
    # print("number of uniq = ", len(uniqIDs))
    # print("intersections: ", set(nonUniqIDs).intersection(uniqIDs))

    with open(args.output, 'w') as ouf:

        # Write to the output file combined chains that had non unique IDs
        for nid in nonUniqIDs:
            print("non unique ID = ", nid) # for verbosity
            nChains = dictChains.get(nid)
            temp = open(args.tempStore, 'w')
            temp.write("\n".join(nChains) + "\n")
            # Close temp file!
            temp.close()

            # Make merged chain string using chainToAxt -> axtChain
            mergedChainStr = runAxtConversions(args.tempStore, args.T2bit, args.Q2bit, args.linearGapVar)
            #print(mergedChainStr.decode('utf-8'))


            # Write to a new chain file changing ID field
            for line in mergedChainStr.decode('utf-8').split('\n')[:-1]:
                if line.startswith('chain'):
                    ouf.write(' '.join(line.split()[:-1])+' {}\n'.format(nid))
                elif line != '\n':
                    ouf.write(line)
                    ouf.write('\n')

        # And write all the rest that was unique to the same file
        for uid in uniqIDs:
            uChains = dictChains.get(uid)
            ouf.write("\n".join(uChains) + "\n")


###############################


def readChainFile(chain_file):
    """Read chain file, return list of chains."""
    with open(chain_file, "r") as f:
        content = f.read().split("\n\n")[:-1]

    return(content)

###############################


def getListNonUniqIDs(allChains):
    # Create a list of all IDs
    allIDs, task_size = [], len(allChains)
    #id_num = Counter()
    dictChains = defaultdict(list)
    for num, chain in enumerate(allChains):
        chain_instance = Chain(chain, parse_blocks=False)
        # allIDs.append(chain_instance.attrs["chain_id"])
        chid = chain_instance.attrs["chain_id"]
        # id_num[chain_instance.attrs["chain_id"]] += 1
        dictChains[chid].append(chain_instance.restore_string())
        #print("{0}/{1}".format(num, task_size))

    # Make a list of non unique IDs
    # nonUniqIDs = set([i for i in allIDs if allIDs.count(i) > 1])
    nonUniqIDs = [k for k, v in dictChains.items() if len(v) > 1]
    uniqIDs = [k for k, v in dictChains.items() if len(v) == 1]
    return nonUniqIDs, uniqIDs, dictChains

###############################


def extractNonUniqByID(id, allChains):
    oneIDchainStr = ''
    # check each of them one-by-one
    for chain in allChains:
        chain_instance = Chain(chain)
        chid = chain_instance.attrs["chain_id"]
        if chid == id:
            chain += '\n\n'
            oneIDchainStr += chain

    return(oneIDchainStr)

###############################


def runAxtConversions(oneIDfile, T2, Q2, linGap):
    commandChainToAxt = "chainToAxt {0} {1} {2} stdout".format(oneIDfile, T2, Q2)
    commandAxtChain = "axtChain -linearGap={2} stdin {0} {1} stdout ".format(T2, Q2, linGap)
    combinedCommand = "{0} | {1}".format(commandChainToAxt, commandAxtChain)
    FNULL = open(os.devnull, 'w')
    print(combinedCommand)

    mergedChainStr = subprocess.check_output(combinedCommand, shell=True, stderr=FNULL)

    return(mergedChainStr)

###############################
# chainToAxt hg38.rheMac8.1942.chain /projects/hillerlab/genome/gbdb-HL/hg38/hg38.2bit \
#  /projects/hillerlab/genome/gbdb-HL/rheMac8/rheMac8.2bit hg38.rheMac8.1942.axt

# axtChain -linearGap=loose  hg38.rheMac8.1942.axt /projects/hillerlab/genome/gbdb-HL/hg38/hg38.2bit \
# /projects/hillerlab/genome/gbdb-HL/rheMac8/rheMac8.2bit hg38.rheMac8.1942.axt.chain


if __name__ == "__main__":
    main()
