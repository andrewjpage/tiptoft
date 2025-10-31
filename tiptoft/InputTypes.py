"""
Input type validation utilities for command-line arguments.

This module provides validation functions for command-line inputs to ensure
they meet the requirements before processing begins.
"""
import os
import argparse


class InputTypes:
    """
    Utility class containing static methods for validating command-line input types.
    
    This class is used by argparse to validate inputs before the main program runs,
    providing early error detection and clear error messages to users.
    """

    def is_fastq_file_valid(filename):
        """
        Validate that a FASTQ input file exists or is stdin.
        
        Args:
            filename (str): Path to the FASTQ file, or "-" for stdin
            
        Returns:
            str: The filename if valid
            
        Raises:
            argparse.ArgumentTypeError: If file doesn't exist and isn't stdin
            
        Note:
            The special filename "-" is accepted to indicate stdin input.
        """
        if not os.path.exists(filename) and filename != "-":
            raise argparse.ArgumentTypeError('Cannot access input file')
        return filename

    def is_plasmid_file_valid(filename):
        """
        Validate that a plasmid database file exists.
        
        Args:
            filename (str): Path to the plasmid database FASTA file
            
        Returns:
            str: The filename if valid
            
        Raises:
            argparse.ArgumentTypeError: If file doesn't exist
        """
        if not os.path.exists(filename):
            raise argparse.ArgumentTypeError('Cannot access input file')
        return filename

    def is_kmer_valid(value_str):
        """
        Validate that k-mer size is within acceptable range.
        
        The k-mer size must be between 7 and 31 inclusive. These limits are based
        on the biological constraints of long-read sequencing error rates and the
        computational efficiency of the k-mer matching algorithm.
        
        Args:
            value_str (str): String representation of the k-mer size
            
        Returns:
            int: The k-mer size as an integer if valid
            
        Raises:
            argparse.ArgumentTypeError: If k-mer is not a valid integer between 7 and 31
            
        Note:
            - Smaller k-mers (< 7) will match too frequently due to sequencing errors
            - Larger k-mers (> 31) may miss matches due to high error rates in long reads
        """
        if value_str.isdigit():
            kmer = int(value_str)
            if kmer >= 7 and kmer <= 31:
                return kmer
        raise argparse.ArgumentTypeError(
            """Invalid Kmer value, it must be an integer between 7 and 31.
 If you need values outside of this range, your probably
 trying to use this software for something it was never
 designed to do.""")
