#!/sw/bin/python3
"""Python class for operating chains."""
import argparse
import sys
import os

__author__ = 'Bogdan Kirilenko, 2018'
__version__ = 0.02


class Chain():
    """
    Toolset for chain processing with python.

    The class may be initiated with a chain file or a string representing content of a chain file.
    It is assumed that one instance -> one chain, that means it should be only one chain in the
    input file.

    After initiating of the class you get an object with attributes corresponding with a chain
    header like "chain_score" or "qStrand". Moreover, it is assumed as default, you get a set
    of blocks-related attributes. It might be called like Chain.block_info[block_num].
    Block-related attributes represent such features as "relative and absolute coordinates" or
    "size". For example, c.blocks[4]["Size"] returns size of the block with index 4.

    By the way, you can skip parsing of the blocks. Just set parse_blocks=False when you call
    the class. In this case all the block-dependend methods shall be unavailable!
    If you wish to parse this information later (by runtime purposes, for example) call
    Chain.block_info().
    """

    def __init__(self, chain_string, parse_blocks=True):
        """Read chain, extract main attrs."""
        # check if it is a file
        if os.path.isfile(chain_string):
            with open(chain_string, 'r') as f:
                self._raw_lines = f.read().split('\n')
        else:  # it is already a content of some file
            self._raw_lines = chain_string.split('\n')

        # get header info
        self.header = self._raw_lines[0].split()
        # check if header is ok:
        if len(self.header) != 13 or self.header[0] != 'chain':
            err_msg = "Error! There are data fields of incorrect type in chain header!\n" \
                      "Should starts with 'chain' and contain 12 params\n" \
                      "The header is: {0}\n".format(' '.join(self.header))
            raise self.__die(err_msg)
        # header is ok: parse it
        self.attrs = self._parse_header()

        # parse_blocks=True as default
        # should be false if user is interested in headers only for time economy
        # if you need it later just call Chain.add_blocks()
        if parse_blocks:
            self.blocks_info()
            self.blocks_parsed = True
        else:
            self.blocks_parsed = False
            sys.stderr.write('Warning! Blocks info is not available in this mode.\n')
            sys.stderr.write('Call Chain.blocks_info() if you need this info.\n')

    def blocks_info(self):
        """Parse blocks separately from __init__."""
        # parse blocks
        self.raw_blocks = [x for x in self._raw_lines[1:] if x != '']
        # get attrs for each block and interval
        self.blocks, self.blocks_iters_num, self.t_intervals, self.q_intervals = self._parse_blocks()
        # the last check (tEnd must be equal with tEnd of the last block)
        self.__check()
        self.blocks_parsed = True

    def _parse_header(self):
        """Parse chain header."""
        try:
            chain_score = int(self.header[1])
            tName = self.header[2]
            tSize = int(self.header[3])
            tStrand = self.header[4]
            tStrand_b = True if self.header[4] == "+" else False
            tStart = int(self.header[5])
            tEnd = int(self.header[6])
            qName = self.header[7]
            qSize = int(self.header[8])
            qStrand = self.header[9]
            qStrand_b = True if self.header[9] == "+" else False
            qStart = int(self.header[10])
            qEnd = int(self.header[11])
            chain_id = int(self.header[12])
            # save and return attributes
            attrs = {'chain_score': chain_score, 'tName': tName, 'tSize': tSize, 'tStrand': tStrand, 'tStart': tStart,
                     'tEnd': tEnd, 'qName': qName, 'qSize': qSize, 'qStrand': qStrand, 'qStart': qStart, 'qEnd': qEnd,
                     'chain_id': chain_id, 'qStrand_b': qStrand_b, 'tStrand_b': tStrand_b}
            return attrs
        except ValueError:  # it was impossible to read the header properly
            err_msg = "Error! There are data fields of incorrect type in chain header!\n" \
                      "The header is: {0}\n".format(' '.join(self.header))
            raise self.__die(err_msg)

    def _parse_blocks(self):
        """Parse blocks."""
        # the last block is so special
        last_block_check_1 = len(self.raw_blocks[-1].split()) != 1
        last_block_check_2 = not self.raw_blocks[-1].isdigit()
        if last_block_check_1 or last_block_check_2:
            # it should be only one number in the last blocks
            err_msg = "Error! The last block is broken!\n" \
                      "expected single number, got: {0}\n".format(self.raw_blocks[-1])
            raise self.__die(err_msg)  # kill if it is wrong

        tail_size = int(self.raw_blocks[-1])

        # so I can iterate by this blocks
        # get initial values for coordinates
        t_flag, q_flag = self.attrs['tStart'], self.attrs['qStart']  # start in target and query genomes
        block_attrs = {}  # put the info for blocks here
        t_intervals, q_intervals = {}, {}  # save intervals info here
        i = -1  # if self.raw_blocks[:-1] == [] --> 'i' referenced before assignment
        for i in range(len(self.raw_blocks[:-1])):
            # for each block get start end in both target and query genomes
            block_info = self.raw_blocks[i].split('\t')  # ['13', '4', '0']
            if len(block_info) != 3:  # check if block is correct
                err_msg = "Error! Wrong number of elements in a block!\n" \
                          "Block num. {0}, expected 3 elems, got: {1}\n".format(i, self.raw_blocks[i])
                raise self.__die(err_msg)

            try:  # check if type of data is correct
                size = int(block_info[0])
                dt = int(block_info[1])
                dq = int(block_info[2])
            except ValueError:
                err_msg = "Error! Data of wrong types in the block!\n" \
                          "Block num. {0}, expected 3 ints, got: {1}\n".format(i, self.raw_blocks[i])
                raise self.__die(err_msg)

            # get start and end positions
            block_start_t, block_start_q = t_flag, q_flag  # where it starts
            block_end_t, block_end_q = t_flag + size, q_flag + size  # where it ends in T and
            # increase flags for the next iteration
            t_flag += size + dt
            q_flag += size + dq
            # save the results for blocks
            block_attrs[i] = {'size': size, 'dt': dt, 'dq': dq,
                              'btstart': block_start_t, 'bqstart': block_start_q,
                              'btend': block_end_t, 'bqend': block_end_q}
            # save interval-related info, at first for T
            t_intervals[i] = (block_start_t, block_end_t)  # for the block as is
            t_intervals[(i, i + 1)] = (block_end_t, t_flag)  # for the gap between
            # the same for Q
            q_intervals[i] = (block_start_q, block_end_q)  # block
            q_intervals[(i, i + 1)] = (block_end_q, q_flag)  # gap
        # save for the tail
        block_attrs[i + 1] = {'size': tail_size, 'dt': 0, 'dq': 0,
                              'btstart': t_flag, 'bqstart': q_flag,
                              'btend': t_flag + tail_size, 'bqend': q_flag + tail_size}
        t_intervals[i + 1] = (t_flag, t_flag + tail_size)
        q_intervals[i + 1] = (q_flag, q_flag + tail_size)
        blocks_iters_num = i + 2  # range(self.blocks_iters_num) - iter by all the blocks
        return block_attrs, blocks_iters_num, t_intervals, q_intervals

    # TODO merge into one func
    def __die(self, err_msg):
        """Die in case of error."""
        sys.stderr.write(err_msg)
        sys.stderr.write("Program finished with exit code -1\n")
        sys.exit(1)

    def __incorrect_order_error(self, err_msg):
        """Def raised in case if the range is wrong in some sence."""
        sys.stderr.write("The order of blocks {0} is wrong!\n".format(err_msg))
        sys.stderr.write("Left block must be <= right block\n")
        sys.stderr.write("Program finished with exit code -1\n")
        sys.exit(1)  # send signal that everything is wrong

    def __blocks_are_not_parsed(self):
        """Raise if blocks were not parsed but called."""
        sys.stderr.write("Error: blocks-dependend method was called.\n")
        sys.stderr.write("Your possible options: \n")
        sys.stderr.write("1) Do not call this class with parse_blocks=False attribute\n")
        sys.stderr.write("2) Call .blocks_info() method before using block-dependend methods.\n")
        sys.exit(1)

    def __check(self):
        """Check that all the blocks and intervals are correct."""
        # get ends and starts according the blocks
        blocks_tStart = self.blocks[0]['btstart']
        blocks_qStart = self.blocks[0]['bqstart']
        max_block = self.blocks_iters_num - 1
        blocks_tEnd = self.blocks[max_block]['btend']
        blocks_qEnd = self.blocks[max_block]['bqend']

        # compare with the values from header, die if something is wrong
        # for tStart
        if self.attrs['tStart'] != blocks_tStart:
            sys.stderr.write("tStart's according chain header and blocks are not equal!\n")
            sys.stderr.write("tStart in header: {0}, tStart in blocks: {1}\n".format(self.attrs['tStart'],
                                                                                     blocks_tStart))
            sys.exit(1)
        # for qStart:
        if self.attrs['qStart'] != blocks_qStart:
            sys.stderr.write("qStart's according chain header and blocks are not equal!\n")
            sys.stderr.write("qStart in header: {0}, qStart in blocks: {1}\n".format(self.attrs['qStart'],
                                                                                     blocks_qStart))
            sys.exit(1)
        # for tEnd:
        if self.attrs['tEnd'] != blocks_tEnd:
            sys.stderr.write("tEnd's according chain header and blocks are not equal!\n")
            sys.stderr.write("tEnd in header: {0}, tEnd in blocks: {1}\n".format(self.attrs['tEnd'],
                                                                                 blocks_tEnd))
            sys.exit(1)
        # for qEnd:
        if self.attrs['qEnd'] != blocks_qEnd:
            sys.stderr.write("qEnd's according chain header and blocks are not equal!\n")
            sys.stderr.write("qEnd in header: {0}, qEnd in blocks: {1}\n".format(self.attrs['qEnd'],
                                                                                 blocks_qEnd))
            sys.exit(1)
        # Nothing? So, 99.99999% everything is fine with input data

    def __apply_shift(self, block_num, upstr, shift):
        """Apply shift for the block number given."""
        if upstr:  # add shith for the block_number estimated
            block_with_shift = block_num + shift
            # check if it is in range
            if block_with_shift > self.blocks_iters_num - 1:
                return "upstream"  # the shift is too big
        else:  # substract shift in this case
            block_with_shift = block_num - shift
            if block_with_shift < 0:
                return "downstream"
        return block_with_shift

    def __locked_method(self):
        """Use to protect unavailable functions."""
        if not self.blocks_parsed:
            raise self.__blocks_are_not_parsed()

    def __block_number(self, coordinate, upstr=True, T=True, shift=0):
        """Return the number of interval corresponding to the certain coordinate."""
        # check if blocks are exist
        self.__locked_method()
        coordinate = int(coordinate)  # just in case if it is a string
        block_num = None  # default value, return if nothing found
        # get the interval for check
        if T:  # coordinate is for target genome
            current_interval = self.t_intervals
            # check if coordinate is out of chain scope
            if coordinate < self.attrs['tStart']:  #
                sys.stderr.write("WARNING: coord {0} is less than {1}\n".format(coordinate, self.attrs['tStart']))
                return 'downstream'  # stands for "lower"
            elif coordinate > self.attrs['tEnd']:
                sys.stderr.write("WARNING: coord {0} is more than {1}\n".format(coordinate, self.attrs['tEnd']))
                return 'upstream'
            # and some very simple cases
            elif coordinate == self.attrs['tStart']:
                return 0
            elif coordinate == self.attrs['tEnd']:
                return self.blocks_iters_num - 1

        else:  # otherwise, for the query
            current_interval = self.q_intervals
            # check the same there
            if coordinate < self.attrs['qStart']:  #
                sys.stderr.write("WARNING: coord {0} is less than {1}\n".format(coordinate, self.attrs['qStart']))
                return 'downstream'  # stands for "lower"
            elif coordinate > self.attrs['qEnd']:
                sys.stderr.write("WARNING: coord {0} is more than {1}\n".format(coordinate, self.attrs['qEnd']))
                return 'upstream'
            # and some very simple cases
            elif coordinate == self.attrs['qStart']:
                return 0
            elif coordinate == self.attrs['qEnd']:
                return self.blocks_iters_num - 1

        # now go interval-by-interval
        for num, coords in current_interval.items():
            # set the possible values
            block_range = range(coords[0], coords[1])
            # check if it is here
            if coordinate in block_range:
                # num - exactly what do we need
                if type(num) == tuple:  # it is not a block, but a gap in between
                    if upstr:  # we should choose between left and right blocks
                        block_num = num[1]  # if we are looking for right block
                    else:
                        block_num = num[0]
                else:  # it is block
                    block_num = num
                # we don't need to iterate it
                break
        # apply shift
        block_num = self.__apply_shift(block_num, upstr, shift)
        return block_num

    def __grange_blocks(self, genomic_range, T=True, shift=0):
        """Get block numbers corresponding a genome range."""
        # check if blocks were parsed
        self.__locked_method()
        try:  # check if it is a proper genomic range
            grange_chrom, grange_crd = genomic_range.split(':')
            assert len(grange_crd.split('-')) == 2
            grange_start, grange_end = int(grange_crd.split('-')[0]), int(grange_crd.split('-')[1])
            assert grange_start <= grange_end
            # print(grange_chrom, grange_start, grange_end)
        except Exception:
            err_msg = "The genomic range {0} is wrong!\n" \
                      "It should be 'chrom:start-end', where start <= end.\n" \
                      "Please check the source of the input data.\n".format(genomic_range)
            raise self.__die(err_msg)

        try:
            if T:  # check if the chrom is the same
                assert grange_chrom == self.attrs['tName']
            else:  # and also if we are looking in Q genone
                assert grange_chrom == self.attrs['qName']
        except AssertionError:
            err_msg = "Error! Genomic range {0} contains wrong {1} chromosome.\n" \
                      "It is {2} in target genome and {3} in the query.\n" \
                      "".format(genomic_range, grange_chrom, self.attrs['tName'], self.attrs['qName'])
            raise self.__die(err_msg)

        # get left and right block numbers
        left_block = self.__block_number(grange_start, upstr=False, T=T, shift=shift)
        right_block = self.__block_number(grange_end, upstr=True, T=T, shift=shift)

        # and return it finally
        return left_block, right_block

    def __arrange_subchain(self, left_block, right_block):
        """Make subchain using left and right blocks."""
        # check if blocks were parsed
        if not self.blocks_parsed:
            raise self.__blocks_are_not_parsed()
        # check if there are "down/upstream" blocks
        if left_block == right_block == "downstream":
            # this case looks like:
            # -----l---r----cccccccccccc--------
            # l is left side, r - right one and c represents chain
            error_data = '"left block: {0}; right block: {1}"'.format(left_block, right_block)
            raise self.__incorrect_order_error(error_data)

        elif left_block == right_block == "upstream":
            # there it is:
            # ------ccccccccccccc----l------r---
            # actually we should return nothing
            error_data = '"left block: {0}; right block: {1}"'.format(left_block, right_block)
            raise self.__incorrect_order_error(error_data)

        # check for some weird cases
        if left_block == "upstream" or right_block == "downstream":
            # if they are not equal but the condition is true -->
            # --> it was wrong range actually
            error_data = '"left block: {0}; right block: {1}"'.format(left_block, right_block)
            raise self.__incorrect_order_error(error_data)

        # there are no dangerous cases anymore
        # but left and tight blocks should be numbers
        if left_block == "downstream":
            # possibilities here are:
            # ---l------cccccc(r)ccc------------
            # or
            # ---l---ccccccccccccc-----r--------
            left_block = 0
        # the same for right side
        if right_block == "upstream":
            # looks like:
            # -----cccccc(l)ccccc-----r---------
            # or:
            # --l----cccccccccccc-----r---------
            right_block = self.blocks_iters_num - 1  # the last block

        # it must be respected
        if left_block > right_block:
            error_data = '"left block: {0}; right block: {1}"'.format(left_block, right_block)
            raise self.__incorrect_order_error(error_data)

        subchain_blocks = []
        # we do need a header for this chain
        header_template = 'chain 0 {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10}\n'
        # extract starts and ends
        # starts of the first blocks == starts for all the chain
        new_tStart, new_qStart = self.blocks[left_block]['btstart'], self.blocks[left_block]['bqstart']
        # ends of the last blocks == ends for all the chain
        new_tEnd, new_qEnd = self.blocks[right_block]['btend'], self.blocks[right_block]['bqend']
        # arrange the header
        subchain_header = header_template.format(self.attrs['tName'],
                                                 self.attrs['tSize'],
                                                 self.attrs['tStrand'],
                                                 new_tStart,
                                                 new_tEnd,
                                                 self.attrs['qName'],
                                                 self.attrs['qSize'],
                                                 self.attrs['qStrand'],
                                                 new_qStart,
                                                 new_qEnd,
                                                 self.attrs["chain_id"])
        # let's collect all the blocks
        for bnum in range(left_block, right_block):
            # parse info
            current_block = self.blocks[bnum]
            curr_size = current_block['size']
            curr_dt = current_block['dt']
            curr_dq = current_block['dq']
            # add into an array
            new_block = '{0}\t{1}\t{2}'.format(curr_size, curr_dt, curr_dq)
            subchain_blocks.append(new_block)
        # add the tail
        tail_size = self.blocks[right_block]['size']
        subchain_blocks.append(str(tail_size))
        # convert to string
        blocks_lines = '\n'.join(subchain_blocks) + '\n'
        # arrange the output finally
        subchain = subchain_header + blocks_lines
        return subchain

    def make_subchain(self, grange, shift=0, output=''):
        """Get genomic range - return or save to file a subchain."""
        # check if blocks were parsed
        self.__locked_method()
        # define chain_id == 0
        self.attrs["chain_id"] = 0
        # if user needs to write it in file check if it is possible
        if output not in ('', 'stdout', 'object'):
            if not os.path.isfile(output):  # this is wrong output file location
                err_msg = "Error! Wrong output! Location {0} is unreachable.\n".format(output)
                raise self.__die(err_msg)

        # get the numbers of blocks
        lb, rb = self.__grange_blocks(grange, shift=shift)
        # get the corresponding subchain
        subchain = self.__arrange_subchain(lb, rb)
        # save the result
        if output == 'stdout':  # just print, for > redirection
            sys.stdout.write(subchain)
            return None
        elif output == 'object':  # in case if you need a python object
            return subchain
        elif output != '':  # means, it should be a file
            with open(output, 'w') as f:
                f.write(subchain)
                return None
        else:  # if you need a python object
            return subchain

    def change_id(self, new_id):
        """Change chain_id to user-defined."""
        # TODO extend to another attrs
        # check if new_id might be an integer
        try:
            self.attrs["chain_id"] = int(new_id)
        except TypeError:
            self.__die("Error! New if ")
        # arrange new chain
        new_chain = self.__arrange_subchain(left_block="downstream", right_block="upstream")
        return new_chain

    def restore_string(self):
        """Return the string representation of chain as it was."""
        # str_chain = self.__arrange_subchain(left_block="downstream", right_block="upstream")
        str_chain = "\n".join(self._raw_lines) + "\n"  # just make this chain again
        return str_chain


def main():
    """Chain_lib-based program for subchain-cutting purposes."""
    app = argparse.ArgumentParser(description=main.__doc__)
    app.add_argument("chain_file", type=str, help="Chain file. Take into accout it should be one chain in the file.")
    app.add_argument("range", type=str, help="Genomic range for the subchain. Like chr1:1-100.")
    app.add_argument("--output", type=str, default="stdout", help="Output, stdout as default.")
    app.add_argument("--shift", type=int, default=0, help="Number of blocks in up/downstream region.")
    # print help if there are no arguments
    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    # parse args
    args = app.parse_args()
    # call the methods
    chain = Chain(args.chain_file)
    subchain = chain.make_subchain(grange=args.range, output=args.output, shift=args.shift)


if __name__ == '__main__':
    main()
    # Tests:
    # test_file = '5857.chain'
    # test_range = 'chr11:63183500-63184500'
    # original_range = 'chr11:63162697-63239687'
    #
    # c = Chain(test_file, parse_blocks=True)
    # print(c.attrs['chain_score'])
    # c.make_subchain(test_range, output='stdout', shift=1)
