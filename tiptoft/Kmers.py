"""
K-mer extraction and management from DNA sequences.

This module provides the Kmers class for extracting all k-mers from a DNA sequence,
with optional homopolymer compression to handle sequencing errors in long reads.

The homopolymer compression is quite intensive so it can be run with Cython
for performance, with a pure Python fallback.
"""
from homopolymer_compression import homopolymer_compression_of_sequence


class Kmers:
    """
    Extracts and manages k-mers from a DNA sequence.
    
    K-mers are subsequences of length k. This class can optionally apply
    homopolymer compression before k-mer extraction, which reduces consecutive
    identical bases to single bases (e.g., "AAA" -> "A"). This helps match
    k-mers across long reads with homopolymer length errors.
    
    Attributes:
        sequence (str): The DNA sequence (possibly compressed)
        k (int): Length of k-mers to extract
        end (int): Last valid starting position for a k-mer
    """
    
    def __init__(self, sequence, k, homopolyer_compression):
        """
        Initialize Kmers extractor for a sequence.
        
        Args:
            sequence (str): DNA sequence to process
            k (int): K-mer size
            homopolyer_compression (bool): Whether to apply homopolymer compression
                                          before k-mer extraction
        """
        # Apply homopolymer compression if requested
        self.sequence = homopolymer_compression_of_sequence(
            sequence) if homopolyer_compression else sequence
        self.k = k
        self.end = len(self.sequence) - self.k + 1

    def get_all_kmers_counter(self, max_kmer_count=10):
        """
        Get all k-mers as a counter dictionary initialized to 0.
        
        This returns a dictionary where keys are k-mer sequences and values
        are 0. It's used as a starting point for counting k-mer occurrences.
        
        Args:
            max_kmer_count (int): Maximum number of times a k-mer can occur
                                 to be included (filters repetitive k-mers)
            
        Returns:
            dict: Dictionary mapping k-mer sequences to 0
        """
        kmers = self.get_all_kmers_filtered(max_kmer_count)
        return {x: 0 for x in kmers.keys()}

    def get_all_kmers_freq(self, max_kmer_count=10):
        """
        Get all k-mers with their frequency counts.
        
        Args:
            max_kmer_count (int): Maximum number of times a k-mer can occur
                                 to be included (filters repetitive k-mers)
            
        Returns:
            dict: Dictionary mapping k-mer sequences to their frequency counts
        """
        kmers = self.get_all_kmers_filtered(max_kmer_count)
        return {k: len(v) for k, v in kmers.items()}

    def get_all_kmers(self, max_kmer_count=10):
        """
        Get all k-mers with their positions (legacy API).
        
        Args:
            max_kmer_count (int): Maximum number of times a k-mer can occur
                                 to be included
            
        Returns:
            dict: Dictionary mapping k-mer sequences to lists of positions
        """
        return self.get_all_kmers_filtered(max_kmer_count)

    def get_all_kmers_filtered(self, max_kmer_count=10):
        """
        Extract all k-mers and filter out overly repetitive ones.
        
        Repetitive k-mers (those occurring more than max_kmer_count times)
        are excluded because they are likely to cause false positive matches.
        
        Args:
            max_kmer_count (int): Maximum times a k-mer can occur to be kept
            
        Returns:
            dict: Dictionary mapping k-mer sequences to lists of their
                  starting positions in the sequence
                  
        Note:
            Positions are 0-based indices
        """
        kmers = {}

        # Extract all k-mers with their positions
        kmer_sequences = [self.sequence[i:i+self.k]
                          for i in range(0, self.end)]

        # Build dictionary of k-mer positions
        for i, k in enumerate(kmer_sequences):
            if k in kmers:
                kmers[k].append(i)
            else:
                kmers[k] = [i]

        # Filter out overly repetitive k-mers
        filtered_kmers = {k: v for k,
                          v in kmers.items() if len(v) <= max_kmer_count}
        return filtered_kmers

    def get_one_x_coverage_of_kmers(self):
        """
        Get a 1x coverage sampling of k-mers from the sequence.
        
        Instead of extracting overlapping k-mers (sliding window), this
        extracts non-overlapping k-mers at intervals of k bases. This is
        useful for quickly checking if a sequence is present without
        examining every k-mer.
        
        Returns:
            list: List of k-mer sequences sampled at k-base intervals
            
        Example:
            For sequence "ATCGATCG" with k=3:
            Returns ["ATC", "GAT"] instead of ["ATC", "TCG", "CGA", "GAT", "ATC"]
        """
        return [self.sequence[i:i+self.k] for i in range(0, self.end, self.k)]
