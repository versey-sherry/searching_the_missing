#!/usr/bin/env python3
# Name: Sherry Lin

import sys

"""
Homework 1: finding the missing motif

Input/Output: STDIN / STDOUT
"""


# TODO: Read in fa or fna file
# TODO: Produce count for particular k-mers
# TODO: Compute expected Pr(K)
# TODO: Compute mu and sd
# TODO: Output the results
# print('{0:8}:{1:8}\t{2:0d}\t{3:0.2f}\t{4:0.2f}'.format(
# seq, rSeq, count,E,pVal))
# TODO: Add p value by scipy.stats.norm.cdf(z)

class CommandLine():
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.

    '''

    def __init__(self, inOpts=None):
        '''
        CommandLine constructor.

        Implements a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(
            description='Parse arguments for search for the missing',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
        )
        self.parser.add_argument('--minMotif', type=int, help='minimum motif size to evaluate (int>=3)')
        self.parser.add_argument('--maxMotif', type=int, help='maximum motif size to evaluate (int<=8)')
        self.parser.add_argument('--cutoff', type=int, help='Z-score cutoff (negative int)')
        self.parser.add_argument('--kScoring', action='store_true', help='displaying ')
        self.parser.add_argument('--version', action='version', version='%(prog)s 0.1')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


class FastAreader():
    """
    Helper function that returns objects of header and the DNA sequences separately

    """

    def __init__(self, fname=''):
        self.fname = fname

    def doOpen(self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):

        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence


class SearchMissing():
    """
    Algorithms that takes single strand of DNA, count the number of k-mers, output the sequences count and Z-score
    Null model is binomial
    Probability is computed with Markovian(2)
    """

    def __init__(self, sequence, min, max, cutoff, pvalflag):
        """

        Args:
            sequence: DNA sequence
            min: min for the ker to consider min is 3
            max: max for the ker to consider max is 8
            cutoff: cut off for z-score, negative value less than
            pvalflag: if present, compute the normal distributed group
        """
        self.sequence = sequence
        self.min = min
        self.max = max
        self.cutoff = cutoff
        self.pValFlag = pvalflag

    def countSeqRseq(self, sequence, k):
        # Ignore the seq and reverse seq relations and record k-mer counts
        seqDict = {}
        for n in range(len(sequence) + 1 - k):
            tempSeq = sequence[n:n + k]
            if seqDict.get(tempSeq):
                seqDict[tempSeq] += 1
            else:
                seqDict[tempSeq] = 1

        pairDict = {}


        return pairDict


class Usage(Exception):
    '''
    Used to signal a Usage error, evoking a usage statement and eventual exit when raised.
    '''

    def __init__(self, msg):
        self.msg = msg


def main(myCommandLine=None):
    '''
    Implement finding the missing sequence

    '''

    myCommandLine = CommandLine()  # read options from the command line

    try:
        print(myCommandLine.args)  # print the parsed argument string .. as there is nothing better to do
    except Usage as err:
        print(err.msg)

    # Get the commandline arguments
    min = myCommandLine.args.minMotif
    max = myCommandLine.args.maxMotif
    cutoff = myCommandLine.args.cutoff
    pValFlag = myCommandLine.args.kScoring

    # print(min, max, cutoff, pValFlag)

    fastaFile = FastAreader().readFasta()
    for header, sequence in fastaFile:
        # print('header is', header)
        # print('seq is', sequence[0:10])
        print(SearchMissing(sequence, min, max, cutoff, pValFlag))


if __name__ == "__main__":
    main()
