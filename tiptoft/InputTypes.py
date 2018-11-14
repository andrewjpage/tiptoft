'''check the input types from the command line script'''
import os
import argparse


class InputTypes:

    '''The input file should exist'''
    def is_fastq_file_valid(filename):
        if not os.path.exists(filename) and filename != "-":
            raise argparse.ArgumentTypeError('Cannot access input file')
        return filename

    '''The input file should exist'''
    def is_plasmid_file_valid(filename):
        if not os.path.exists(filename):
            raise argparse.ArgumentTypeError('Cannot access input file')
        return filename

    def is_kmer_valid(value_str):
        if value_str.isdigit():
            kmer = int(value_str)
            if kmer >= 7 and kmer <= 31:
                return kmer
        raise argparse.ArgumentTypeError(
            """Invalid Kmer value, it must be an integer between 7 and 31.
 If you need values outside of this range, your probably
 trying to use this software for something it was never
 designed to do.""")
