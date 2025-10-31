"""
Block identification and merging for k-mer matches.

This module contains the Blocks class which identifies contiguous regions
of k-mer matches within reads and merges nearby blocks to account for
sequencing errors in long reads.
"""


class Blocks:
    """
    Identifies and manages contiguous blocks of k-mer matches.
    
    In long read sequencing, k-mer matches may be interrupted by sequencing
    errors. This class finds blocks of consecutive k-mer matches and can merge
    nearby blocks that are likely part of the same feature.
    
    Attributes:
        k (int): K-mer size used for matching
        min_block_size (int): Minimum size in bases for a block to be considered
        max_gap (int): Maximum gap between blocks (in k-mer units) for merging
        margin (int): Flanking bases to add around blocks
    """
    
    def __init__(self, k, min_block_size, max_gap, margin):
        """
        Initialize Blocks processor with parameters.
        
        Args:
            k (int): K-mer size
            min_block_size (int): Minimum block size in bases
            max_gap (int): Maximum gap in k-mer units for merging blocks
            margin (int): Margin in bases to add around blocks
        """
        self.k = k
        self.min_block_size = min_block_size
        self.max_gap = max_gap  # multiples of the kmer
        self.margin = margin

    def merge_blocks(self, blocks):
        """
        Merge nearby blocks that are likely part of the same feature.
        
        Blocks separated by less than max_gap k-mers are merged together.
        This accounts for small regions of elevated sequencing errors that
        interrupt k-mer matching.
        
        Args:
            blocks (list): List of [start, end] coordinate pairs
            
        Returns:
            list: List of merged [start, end] coordinate pairs
            
        Note:
            Blocks are assumed to be in order by start position
        """
        for i in range(0, len(blocks)-1):
            # Check if blocks are close enough to merge
            if blocks[i][1] + self.max_gap > blocks[i+1][0]:
                # Extend the first block to encompass the second
                if blocks[i][1] < blocks[i+1][1]:
                    blocks[i][1] = blocks[i+1][1]
                # Copy merged block forward
                blocks[i+1][0] = blocks[i][0]
                blocks[i+1][1] = blocks[i][1]
        return blocks

    def find_largest_block(self, sequence_hits):
        """
        Find the largest contiguous block of k-mer matches.
        
        This identifies all blocks, merges nearby ones, and returns the
        coordinates of the largest resulting block.
        
        Args:
            sequence_hits (list): Array where each element is the count of
                                 k-mer matches at that position
            
        Returns:
            tuple: (start, end) coordinates of largest block, or (0, 0) if
                   no block meets minimum size requirement
        """
        blocks = self.find_all_blocks(sequence_hits)
        merged_blocks = self.merge_blocks(blocks)

        largest_block = 0
        lbi = 0  # largest block index

        # Find the largest block
        for i, block in enumerate(merged_blocks):
            block_size = block[1] - block[0]
            if block_size > largest_block:
                lbi = i
                largest_block = block_size

        # Check if largest block meets minimum size
        if largest_block < (self.min_block_size/self.k):
            return 0, 0

        return merged_blocks[lbi][0], merged_blocks[lbi][1]

    def find_all_blocks(self, sequence_hits):
        """
        Identify all contiguous blocks of k-mer matches in a sequence.
        
        A block is a contiguous region where k-mer matches occur (value > 0).
        
        Args:
            sequence_hits (list): Array where each element is the count of
                                 k-mer matches at that position
            
        Returns:
            list: List of [start, end] coordinate pairs for each block
            
        Note:
            Coordinates are in k-mer units, not bases
        """
        blocks = []
        in_block = False
        current_block_start = 0
        
        for i, val_count in enumerate(sequence_hits):
            if not in_block and val_count > 0:
                # Start of a new block
                in_block = True
                current_block_start = i
            elif in_block and val_count == 0:
                # End of current block
                in_block = False
                blocks.append([current_block_start, i])

        # Handle case where sequence ends while in a block
        if in_block:
            blocks.append([current_block_start, len(sequence_hits)])

        return blocks

    def adjust_block_start(self, block_start):
        """
        Adjust block start coordinate with margin and bounds checking.
        
        Converts from k-mer units to bases, subtracts margin, and ensures
        the result is not negative.
        
        Args:
            block_start (int): Block start position in k-mer units
            
        Returns:
            int: Adjusted start position in bases (>= 0)
        """
        block_start *= self.k
        if block_start - self.margin < 0:
            block_start = 0
        else:
            block_start -= self.margin
        return block_start

    def adjust_block_end(self, block_end, seq_length):
        """
        Adjust block end coordinate with margin and bounds checking.
        
        Converts from k-mer units to bases, adds margin, and ensures
        the result does not exceed the sequence length.
        
        Args:
            block_end (int): Block end position in k-mer units
            seq_length (int): Length of the sequence in bases
            
        Returns:
            int: Adjusted end position in bases (<= seq_length)
        """
        block_end *= self.k
        if block_end + self.margin > seq_length:
            block_end = seq_length
        else:
            block_end += self.margin
        return block_end
