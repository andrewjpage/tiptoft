"""
Database downloader driver for TipToft.

This module provides the TipToftDatabaseDownloader class which coordinates
downloading the latest plasmid typing database from the PlasmidFinder repository.
"""
import logging
from tiptoft.RefGenesGetter import RefGenesGetter


class TipToftDatabaseDownloader:
    """
    Coordinates downloading plasmid database from PlasmidFinder.
    
    This class manages the process of downloading the latest plasmid
    incompatibility group and replicon typing sequences from the
    PlasmidFinder database repository.
    
    Attributes:
        logger: Python logger for diagnostic messages
        output_prefix (str): Prefix for output files (.fa and .tsv)
        verbose (bool): Whether to enable debug logging
    """
    
    def __init__(self, options):
        """
        Initialize database downloader with command-line options.
        
        Args:
            options: Parsed command-line arguments with attributes:
                    - output_prefix: Prefix for output files
                    - verbose: Enable debug logging
        """
        self.logger = logging.getLogger(__name__)
        self.output_prefix = options.output_prefix
        self.verbose = options.verbose

        # Configure logging level
        if self.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.ERROR)

    def run(self):
        """
        Execute the database download process.
        
        Downloads plasmid typing sequences from PlasmidFinder and saves
        them to files with the specified output prefix.
        
        The method will create two files:
        - <output_prefix>.fa: FASTA file with sequences
        - <output_prefix>.tsv: TSV file with metadata
        
        If no output prefix is provided, logs an error and returns.
        """
        refgenes = RefGenesGetter(verbose=self.verbose)

        if self.output_prefix:
            refgenes.run(self.output_prefix)
        else:
            self.logger.error("Please check the input parameters")
