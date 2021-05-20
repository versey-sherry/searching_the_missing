#!/usr/bin/env python3
# Name: Sherry Lin
#

import sys
from math import sqrt
from collections import defaultdict
from scipy.stats import norm
import time


"""
Homework 1: finding the missing motif

Input/Output: STDIN / STDOUT

This script aims at finding the underrepresented motif by the z-score

Examples:
Order sequence by the z-score with respect to the whole sequence
python missingMotif.py --minMotif 3 --maxMotif 8 --cutoff -5  < xx.fna  > output.out

Order sequence by the z-score with respect to the k mer group and the p value
python missingMotif.py --minMotif 3 --maxMotif 8 --cutoff 0 --kScoring  < xx.fna  > output.out

"""

# print('{0:8}:{1:8}\t{2:0d}\t{3:0.2f}\t{4:0.2f}'.format(
# seq, rSeq, count,E,pVal))

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
        self.parser.add_argument('--kScoring', action='store_true', help='Using p value to score the motif and display p value')
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
    The output is ordered by z score of a DNA sequence by k-mer group

    To estimate the z score for a particular motif, we assume the distribution is binomial
    and the probability and expected value for the motif is estimated by Markovian(2)

    To estimate the z score for a particular motif in a k-mer group, we assume that within the K-mer group,
    the motif count, normalized by the expected value should be normally distributed.
    expected value and standard deviation within k-mer group is approximated by sample mean and sample deviation

    For my computation, I assume that the motif must appear in the DNA sequence
    at least once to be considered as a relevant motif
    Therefore, I didn't permutate all possible k-mer and my results don't include motif that has 0 count
    My mean of normalized values would be slightly greater than the mean of normalized values that include 0 counts.
    However, the order of non-zero count motifs would be the same as the results that includes 0 counts

    """

    def __init__(self, sequences, min, max, cutoff, pValFlag):
        """
        Initialize the objects and add all necessary attributes, including the dictionary of all needed k-mer
        Including such a dictionary will avoid computing the dictionary multiple times

        Args:
            sequence: DNA sequences, a list containing all fasta sequences
            min: min for the ker to consider min is 3
            max: max for the ker to consider max is 8
            cutoff: cut off for z-score, negative value less than
            pvalflag: if present, compute the normal distributed group
        """
        self.sequences = sequences
        self.min = min
        self.max = max
        self.cutoff = cutoff
        self.pValFlag = pValFlag
        self.kDict = self.genAllk()

    def countkSeqRseq(self, k):
        """
        This function counts the motif of length k
        The algorithm counts the sequence and the reverse complement sequence equivalently
        and return the total count of sequence and its reverse sequence as the count for the pair

        Args:
            k: the length of the motif we are counting

        Returns:
            a dictionary contains the Seq:rSeq pair count for all the k-mer

        """
        # Ignore the seq and reverse seq relations and record k-mer counts
        seqDict = {}

        for sequence in self.sequences:
            for n in range(len(sequence) + 1 - k):
                tempSeq = sequence[n:n + k]
                if seqDict.get(tempSeq):
                    seqDict[tempSeq] += 1
                else:
                    seqDict[tempSeq] = 1

        # Pairwise dictionary with keys seq:rSeq
        pairDict = {}
        # sort the list so the pairs will be in alpha order
        seqList = sorted(seqDict.keys())
        # add the seq to counted list to avoid double count
        counted = []

        # Using string.translate() method to generate complement sequence
        # https://www.programiz.com/python-programming/methods/string/translate
        # ascii table https://www.ascii-code.com/
        translation = {65: 84, 84: 65, 67: 71, 71: 67}
        # sum = 0
        for seq in seqList:
            rSeq = seq.translate(translation)[::-1]
            if seq not in counted:
                value = seqDict[seq]
                seqDict[seq] = 0
                counted.append(seq)
                # Reverse complement sequence may not be in the dictionary
                if seqDict.get(rSeq):
                    value += seqDict[rSeq]
                    seqDict[rSeq] = 0
                counted.append(rSeq)
                key = ':'.join(sorted([seq, rSeq]))
                pairDict[key] = value
                # print(seq, ':', rSeq, value)
                # sum += value
        # print('total', len(self.sequence)+1-k, 'sum is', sum)
        # sum = 0
        # for item, value in pairDict.items():
        #     sum+=value
        # print('sanity check sum', sum)
        return pairDict

    def genAllk(self):
        """
        Use the input arguments to generate all k-mers from min to max count in the sequence
        Returns:
            A dictionary with k being the key, and a dictionary that contains the seq:rSeq pairs counts as the value
            This method generate all needed k-mers from the sequence. Since all k-mer with this method will appear at least 1
            in the DNA sequence, there is no divide by 0 issue.
        """
        kDict = {}
        for k in range(self.min - 2, self.max + 1):
            # print('generating k-mer dictionary for length', k)
            kDict[k] = self.countkSeqRseq(k)
        return kDict

    def zScore(self, targetSeq):
        """
        Function that computes expected value of the motif by Markovian(2) and assume the distribution to be binomial
        Then the function compute the z score for the motif

        Args:
            targetSeq: the k-mer that we are looking for zScore

        Returns:
            targetSeq: target sequence reverse complement sequence pair
            countK: the count for the target sequence
            mu: expected count for the sequence computed by Markovian(2)
            zScore: zScore for the sequence
        """
        # length of DNA sequence is significantly large so n = N-k+1 can be approximated length of n
        # n = len(self.sequence)
        n = sum([len(sequence) for sequence in self.sequences])

        prefix = targetSeq[0:-1]
        suffix = targetSeq[1:]
        mid = targetSeq[1:-1]

        # generate reverse complement sequence and convert to alpha order pairs
        translation = {65: 84, 84: 65, 67: 71, 71: 67}

        rTarget = targetSeq.translate(translation)[::-1]
        targetSeq = ':'.join(sorted([targetSeq, rTarget]))
        countK = self.kDict[len(rTarget)].get(targetSeq)

        rPrefix = prefix.translate(translation)[::-1]
        prefix = ':'.join(sorted([prefix, rPrefix]))
        countPrefix = self.kDict[len(rPrefix)].get(prefix)

        rSuffix = suffix.translate(translation)[::-1]
        suffix = ':'.join(sorted([suffix, rSuffix]))
        countSuffix = self.kDict[len(rSuffix)].get(suffix)

        rMid = mid.translate(translation)[::-1]
        mid = ':'.join(sorted([mid, rMid]))
        countMid = self.kDict[len(rMid)].get(mid)

        # print('target',targetSeq, countK, 'prefix', prefix, countPrefix, 'suffix', suffix, countSuffix, 'mid', mid, countMid)
        mu = (countPrefix * countSuffix) / countMid
        prK = mu / n
        sd = sqrt(mu * (1 - prK))

        if sd == 0:
            zScore = 0
        else:
            zScore = (countK - mu) / sd

        return targetSeq, countK, mu, zScore

    def genzScore(self):
        """
        loop through the DNA sequence and find the Z Score for all the k-mers and output the zScore
        Returns:
            a dictionary with count of the target sequence, expected and zScore
        """
        resultDict = {}
        for k in range(self.min, self.max + 1):
            for sequence in self.sequences:
                for n in range(len(sequence) + 1 - k):
                    tempSeq = sequence[n:n + k]
                    targetSeq, countK, mu, zScore = self.zScore(tempSeq)
                    resultDict[targetSeq] = [countK, mu, zScore]
        return resultDict

    def pVal(self):
        """
        Use the Z score generated from the previous function, normalized the sequence count by the expected value
        compute z scores for the normalized values by kmer groups
        Returns:
            a dictionary with count of the target sequence, expected and z score, and normalized score,
            z score for normalized values, p value
        """
        zScoreResults = self.genzScore()
        # Normalize Score
        for key, values in zScoreResults.items():
            # generate normal score by computing count/expected
            # 0th is count, 1st is mu, 2nd is z score, 3th is normalized score
            values.append(values[0]/values[1])
        # print('z score results', zScoreResults)

        # Extract the normalized values
        groupNorm = defaultdict(list)
        # print(kStats)
        for key, values in zScoreResults.items():
            k = key.find(':')
            groupNorm[k].append(values[3])
        # print('Group Norm', groupNorm)

        # compute mu and sd
        kStats = defaultdict(list)
        for key, values in groupNorm.items():
            mu = sum(values)/len(values)
            sd = sqrt(sum([value**2 for value in values])/len(values) - mu**2)
            kStats[key] = [mu, sd]
        # print('kStats', kStats)

        for key, values in zScoreResults.items():
            k = key.find(':')
            mu, sd = kStats[k]
            #compute z-score for normalized values
            # print(values)
            if sd ==0:
                normalizedz = 0
            else:
                normalizedz = (values[3]-mu)/sd
            # 0th is count, 1st is mu, 2nd is z score, 3th is normalized score, 4th is the z score for normalized values
            values.append(normalizedz)
            # 0th is count, 1st is mu, 2nd is z score, 3th is normalized score
            # 4th is the z score for normalized values, 5th is the p value
            values.append(norm.cdf(normalizedz))
        return zScoreResults

    def printReuslts(self):
        """
        pretty printing the results
        Returns:

        """
        results = self.pVal()
        # print(results.items())
        # https://stackoverflow.com/questions/613183/how-do-i-sort-a-dictionary-by-value
        # https://stackoverflow.com/questions/4233476/sort-a-list-by-multiple-attributes/4233482
        if self.pValFlag:
            results = {key: values for key, values in sorted(results.items(), key=lambda item: (-len(item[0]), item[1][4]))}

        else:
            results = {key: values for key, values in sorted(results.items(), key=lambda item: (-len(item[0]), item[1][2]))}

        # print(results)
        n = sum([len(sequence) - self.max for sequence in self.sequences])
        print('N =', n)
        for key, values in results.items():
            index = key.find(':')
            if self.pValFlag:
                if values[4] < self.cutoff:
                    print('{0:8}:{1:8}\t{2:0d}\t{3:0.2f}\t{4:0.2f}\t{5:0.2f}'.format(
                        key[0:index], key[index+1:], values[0], values[1], values[4], values[5]))
            else:
                if values[2] < self.cutoff:
                    print('{0:8}:{1:8}\t{2:0d}\t{3:0.2f}\t{4:0.2f}'.format(
                        key[0:index], key[index+1:], values[0], values[1], values[2]))


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

    try:
        myCommandLine = CommandLine()  # read options from the command line
        #print(myCommandLine.args)  # print the parsed argument string .. as there is nothing better to do
    except Usage as err:
        print(err.msg)

    # Get the commandline arguments
    min = myCommandLine.args.minMotif
    max = myCommandLine.args.maxMotif
    cutoff = myCommandLine.args.cutoff
    pValFlag = myCommandLine.args.kScoring

    # print(min, max, cutoff, pValFlag)

    fastaFile = FastAreader().readFasta()
    # store all sequence in a list
    sequences = []
    for header, sequence in fastaFile:
        # print('header is', header)
        # print('seq is', sequence)
        # print(len(sequence))
        sequences.append(sequence)

    searchSequence = SearchMissing(sequences, min, max, cutoff, pValFlag)
    searchSequence.printReuslts()
    #print(searchSequence.kDict)

if __name__ == "__main__":
    #start = time.time()
    main()
    #print('time consumed is', time.time() - start)
