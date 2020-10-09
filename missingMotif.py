#!/usr/bin/env python3
#Name: Sherry Lin

"""
Homework 1: finding the missing motif

Input/Output: STDIN / STDOUT
"""
#TODO: Read in fa or fna file
#TODO: Produce count for particular k-mers
#TODO: Compute expected Pr(K)
#TODO: Compute mu and sd
#TODO: Output the results
# print('{0:8}:{1:8}\t{2:0d}\t{3:0.2f}\t{4:0.2f}'.format(
# seq, rSeq, count,E,pVal))
#TODO: Add p value by scipy.stats.norm.cdf(z)

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
        self.parser.add_argument('--option', action='store', help='foo help')
        self.parser.add_argument('--minMotif', type=int, help='minimum motif size to evaluate (int>=3)')
        self.parser.add_argument('--maxMotif', type=int, help='maximum motif size to evaluate (int<=8)')
        self.parser.add_argument('--cutoff', type=int, help='Z-score cutoff (negative int)')
        self.parser.add_argument('--kScoring', action='store_true', )
        self.parser.add_argument('-l', '--list', action='append', nargs='?', help='list help')  # allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


class Usage(Exception):
    '''
    Used to signal a Usage error, evoking a usage statement and eventual exit when raised.
    '''

    def __init__(self, msg):
        self.msg = msg


def main(myCommandLine=None):
    '''
    Implement the Usage exception handler that can be raised from anywhere in process.

    '''
    if myCommandLine is None:
        myCommandLine = CommandLine()  # read options from the command line
    else:
        myCommandLine = CommandLine(
            myCommandLine)  # interpret the list passed from the caller of main as the commandline.

    try:

        print(myCommandLine.args)  # print the parsed argument string .. as there is nothing better to do

        if myCommandLine.args.requiredBool:
            print('requiredBool is', str(myCommandLine.args.requiredBool))  ## this is just an example
        else:
            pass
        raise Usage(
            'testing')  # this is an example of how to raise a Usage exception and you can include some text that will get printed. Delete this is you dont need it

    except Usage as err:
        print(err.msg)


if __name__ == "__main__":
#   main(['-r'])  # this would make this program.py behave as written

# if you want to make use of a test commandLine, you could do it like this ( notice that it is a list of strings ):
#   main([ '--bool',
#          '--character=b',
#          '-i=10'])