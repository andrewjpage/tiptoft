"""
FASTA database loading and k-mer indexing.

This module loads plasmid typing sequences from a FASTA file and builds
k-mer indices for efficient matching against sequencing reads.
"""
from Bio import SeqIO
from tiptoft.Kmers import Kmers


class Fasta:
    """
    Loads plasmid database and builds k-mer indices.
    
    This class reads a FASTA file containing plasmid typing sequences,
    extracts all k-mers from each sequence, and creates lookup structures
    for efficient matching during read processing.
    
    Attributes:
        logger: Python logger for diagnostic messages
        filename (str): Path to FASTA database file
        k (int): K-mer size
        homopolyer_compression (bool): Whether to compress homopolymers
        max_kmer_count (int): Max occurrences of a k-mer to include
        sequences_to_kmers (dict): Map of sequence ID to k-mer counter
        sequences_to_kmers_count (dict): Map of sequence ID to k-mer frequencies
        all_kmers (dict): All k-mers across database with occurrence counts
        kmers_to_genes (dict): Map of k-mer to list of sequences containing it
        kmer_keys_set (set): Set of all k-mers for fast lookup
    """
    
    def __init__(
            self,
            logger,
            filename,
            k,
            homopolyer_compression,
            max_kmer_count=10):
        """
        Initialize Fasta database loader.
        
        Args:
            logger: Python logger instance
            filename (str): Path to FASTA database file
            k (int): K-mer size for indexing
            homopolyer_compression (bool): Apply homopolymer compression
            max_kmer_count (int): Maximum times a k-mer can occur to be included
        """
        self.logger = logger
        self.filename = filename
        self.k = k
        self.homopolyer_compression = homopolyer_compression
        self.max_kmer_count = max_kmer_count

        # Build all k-mer indices during initialization
        self.sequences_to_kmers = self.sequence_kmers('get_all_kmers_counter')
        self.sequences_to_kmers_count =\
            self.sequence_kmers('get_all_kmers_freq')
        self.all_kmers = self.all_kmers_in_file()
        self.kmers_to_genes = self.all_kmers_to_seq_in_file()
        self.kmer_keys_set = set(self.all_kmers.keys())

    def sequence_kmers(self, kmer_action='get_all_kmers_counter'):
        """
        Extract k-mers from all sequences in the FASTA file.
        
        Processes each sequence and extracts k-mers using the specified
        method (counter or frequency).
        
        Args:
            kmer_action (str): Method name to call on Kmers object
                              ('get_all_kmers_counter' or 'get_all_kmers_freq')
            
        Returns:
            dict: Map of sequence ID to k-mer dictionary
            
        Note:
            Assumes sequence IDs are unique in the FASTA file
        """
        seq_counter = 0
        kmer_to_sequences = {}
        
        # Process each sequence in the FASTA file
        for record in SeqIO.parse(self.filename, "fasta"):
            kmers = Kmers(str(record.seq), self.k, self.homopolyer_compression)
            # Extract k-mers using the specified method
            kmer_to_sequences[record.id] = getattr(kmers, kmer_action)(
                max_kmer_count=self.max_kmer_count)

            seq_counter += 1

        return kmer_to_sequences

    def all_kmers_to_seq_in_file(self):
        """
        Create reverse index from k-mers to sequences.
        
        Builds a dictionary mapping each k-mer to the list of sequences
        (genes) that contain it. This enables finding which genes might
        be present based on k-mer matches.
        
        Returns:
            dict: Map of k-mer to list of sequence IDs containing that k-mer
        """
        kmers_to_genes = {}
        
        for seq_name, kmer_counts in self.sequences_to_kmers.items():
            for kmer, count in kmer_counts.items():
                # Add this sequence to the list for this k-mer
                if kmer not in kmers_to_genes:
                    kmers_to_genes[kmer] = []
                kmers_to_genes[kmer].append(seq_name)

        return kmers_to_genes

    def all_kmers_in_file(self):
        """
        Extract all unique k-mers from the database with occurrence counts.
        
        Creates a dictionary of all k-mers found across all sequences,
        counting how many sequences each k-mer appears in. This is used
        for filtering common k-mers during read matching.
        
        Returns:
            dict: Map of k-mer to number of sequences containing it
        """
        self.logger.info("Finding all k-mers in plasmid FASTA file")
        all_kmers = {}
        
        for seq_name, kmer_counts in self.sequences_to_kmers.items():
            for kmer, count in kmer_counts.items():
                if kmer in all_kmers:
                    all_kmers[kmer] += 1
                else:
                    all_kmers[kmer] = 1

        return all_kmers
