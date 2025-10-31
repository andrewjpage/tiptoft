"""
FASTQ read representation and manipulation.

This module contains the Read class for handling FASTQ sequence reads,
including operations like reverse complementing and subsequence extraction.

Derived from https://github.com/sanger-pathogens/Fastaq
Author: Martin Hunt
Copyright (c) 2013 - 2017 by Genome Research Ltd.
GNU GPL version 3
"""


class Read:
    """
    Represents a single FASTQ sequence read.
    
    A FASTQ read consists of an identifier, the nucleotide sequence,
    and quality scores for each base.
    
    Attributes:
        id (str): Read identifier (without @ prefix)
        seq (str): Nucleotide sequence
        qual (str): Quality scores (same length as seq)
    """
    
    def __init__(self, id=None, seq=None, qual=None):
        """
        Initialize a Read object.
        
        Args:
            id (str, optional): Read identifier
            seq (str, optional): Nucleotide sequence
            qual (str, optional): Quality scores string
        """
        self.id = id
        self.seq = seq
        self.qual = qual

    def subsequence(self, start, end):
        """
        Extract a subsequence from the read.
        
        Creates a new Read object containing only the specified range
        of the sequence and quality scores.
        
        Args:
            start (int): Start position (0-based, inclusive)
            end (int): End position (0-based, exclusive)
            
        Returns:
            Read: New Read object with subsequence and modified ID
            
        Note:
            The new read's ID will be appended with "_start_end"
        """
        return Read(
            id=self.id+"_"+str(start)+"_"+str(end),
            seq=self.seq[start:end],
            qual=self.qual[start:end])

    def __str__(self):
        """
        Format read as FASTQ text.
        
        Returns:
            str: Four-line FASTQ format string with newlines
        """
        return '@' + self.id + '\n' + self.seq + '\n+\n' + self.qual + '\n'

    def reverse_complement_sequence(self):
        """
        Calculate reverse complement of the sequence.
        
        Returns:
            str: Reverse complement of the sequence
            
        Note:
            Handles both uppercase and lowercase nucleotides
        """
        return self.seq.translate(str.maketrans("ATCGatcg", "TAGCtagc"))[::-1]

    def reverse_read(self):
        """
        Create a reverse-complemented copy of the read.
        
        Generates a new Read with the reverse complement sequence and
        reversed quality scores. The ID is appended with "_reverse".
        
        Returns:
            Read: New Read object with reverse-complemented sequence
            
        Note:
            Quality scores are reversed but not recalibrated
        """
        return Read(
            id=self.id+"_reverse",
            seq=self.reverse_complement_sequence(),
            qual=self.qual)

    def get_next_from_file(self, f):
        """
        Read the next FASTQ record from an open file handle.
        
        This method populates the current Read object with data from
        the file handle. It expects standard 4-line FASTQ format.
        
        Args:
            f (file): Open file handle positioned at a FASTQ record
            
        Returns:
            Read|bool: Self if read successfully, False if EOF reached
            
        Note:
            - Skips empty lines
            - Strips whitespace from ID, sequence, and quality lines
        """
        line = f.readline()

        # Skip empty lines
        while line == '\n':
            line = f.readline()

        # Check for end of file
        if not line:
            return False

        # Parse the four FASTQ lines
        self.id = line.rstrip()[1:]  # Remove @ prefix
        line = f.readline()

        self.seq = line.strip()
        line = f.readline()  # Skip '+' line
        line = f.readline()

        self.qual = line.rstrip()
        return self
