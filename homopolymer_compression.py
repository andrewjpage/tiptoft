"""Pure Python fallback for homopolymer compression.
This module provides a Python implementation when the C extension is not available."""

def homopolymer_compression_of_sequence(sequence):
    """
    Compress homopolymer runs in a DNA sequence.
    
    A homopolymer run is a consecutive stretch of the same nucleotide.
    This function reduces each run to a single base.
    
    Args:
        sequence (str): DNA sequence to compress
        
    Returns:
        str: Compressed sequence with homopolymer runs reduced to single bases
        
    Example:
        >>> homopolymer_compression_of_sequence('AAATTTGGG')
        'ATG'
    """
    previous_base = ''
    compressed_sequence = []
    for base in sequence:
        if base != previous_base:
            previous_base = base
            compressed_sequence.append(base)
    return ''.join(compressed_sequence)
