"""
Main driver class for the TipToft plasmid detection workflow.

This module coordinates the overall plasmid incompatibility group prediction
process by managing the Fasta database loading and Fastq read processing.
"""
import logging
import os
import sys
import time
from importlib.resources import files
from tiptoft.Fasta import Fasta
from tiptoft.Fastq import Fastq


class TipToft:
    """
    Main coordinator for plasmid typing from long reads.
    
    This class manages the workflow of:
    1. Loading the plasmid database (FASTA)
    2. Processing input reads (FASTQ)
    3. Matching k-mers between reads and database
    4. Reporting detected plasmid replicons
    
    Attributes:
        logger: Python logger for diagnostic messages
        plasmid_data (str): Path to plasmid database FASTA file
        input_fastq (str): Path to input FASTQ file (or "-" for stdin)
        kmer (int): K-mer size for matching
        verbose (bool): Whether to enable debug logging
        min_fasta_hits (int): Minimum k-mer matches to consider a read
        print_interval (int): Print results every N reads (None = end only)
        output_file (str): Output file path (None = stdout)
        filtered_reads_file (str): Path to save matching reads
        max_gap (int): Max gap between blocks for merging (in k-mer units)
        min_block_size (int): Minimum block size in bases
        margin (int): Flanking bases around blocks
        start_time (int): Unix timestamp when analysis started
        min_kmers_for_onex_pass (int): Min k-mers for initial pass
        min_perc_coverage (int): Minimum coverage % to report
        homopolyer_compression (bool): Whether to compress homopolymers
        max_kmer_count (int): Max times a k-mer can occur to be used
        no_gene_filter (bool): Whether to disable gene filtering
    """
    
    def __init__(self, options):
        """
        Initialize TipToft with command-line options.
        
        Args:
            options: Parsed command-line arguments object with attributes
                    matching the parameters documented above
                    
        Raises:
            SystemExit: If filtered_reads_file already exists
        """
        self.logger = logging.getLogger(__name__)
        
        # Basic input/output settings
        self.plasmid_data = options.plasmid_data
        self.input_fastq = options.input_fastq
        self.kmer = options.kmer
        self.verbose = options.verbose
        self.min_fasta_hits = options.min_fasta_hits
        self.print_interval = options.print_interval
        self.output_file = options.output_file
        self.filtered_reads_file = options.filtered_reads_file
        
        # Advanced algorithm parameters
        self.max_gap = options.max_gap
        self.min_block_size = options.min_block_size
        self.margin = options.margin
        self.start_time = int(time.time())
        self.min_kmers_for_onex_pass = options.min_kmers_for_onex_pass
        self.min_perc_coverage = options.min_perc_coverage
        
        # K-mer processing options
        if options.no_hc_compression:
            self.homopolyer_compression = False
        else:
            self.homopolyer_compression = True
        self.max_kmer_count = options.max_kmer_count
        self.no_gene_filter = options.no_gene_filter

        # Use bundled database if none specified
        if self.plasmid_data is None:
            self.plasmid_data = str(files('tiptoft').joinpath('data', 'plasmid_data.fa'))

        # Validate output files don't already exist
        # Note: Output file check is commented out to allow overwriting
        #if self.output_file and os.path.exists(self.output_file):
        #    self.logger.error(
        #        "The output file already exists, "
        #        "please choose another filename: "
        #        + self.output_file)
        #    sys.exit(1)

        if(self.filtered_reads_file and
                os.path.exists(self.filtered_reads_file)):
            self.logger.error(
                "The output filtered reads file already exists,"
                " please choose another filename: " + self.filtered_reads_file)
            sys.exit(1)

        # Configure logging level
        if self.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.ERROR)

    def run(self):
        """
        Execute the plasmid detection workflow.
        
        This method:
        1. Loads the plasmid database and extracts k-mers
        2. Processes the input FASTQ file
        3. Matches k-mers between reads and database
        4. Reports detected plasmid incompatibility groups
        
        The results are written to stdout or the specified output file.
        """
        self.logger.info("Starting analysis")
        
        # Load plasmid database and build k-mer index
        fasta = Fasta(self.logger,
                      self.plasmid_data,
                      self.kmer,
                      self.homopolyer_compression,
                      max_kmer_count=self.max_kmer_count)
        
        # Process FASTQ reads and match against database
        fastq = Fastq(self.logger,
                      self.input_fastq,
                      self.kmer,
                      fasta.all_kmers_in_file(),
                      self.min_fasta_hits,
                      self.print_interval,
                      self.output_file,
                      self.filtered_reads_file,
                      fasta,
                      self.homopolyer_compression,
                      max_gap=self.max_gap,
                      min_block_size=self.min_block_size,
                      margin=self.margin,
                      start_time=self.start_time,
                      min_kmers_for_onex_pass=self.min_kmers_for_onex_pass,
                      min_perc_coverage=self.min_perc_coverage,
                      max_kmer_count=self.max_kmer_count,
                      no_gene_filter=self.no_gene_filter)
        fastq.read_filter_and_map()
