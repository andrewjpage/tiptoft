'''Driver class for the tiptoft script'''
import logging
import os
import sys
import time
import pkg_resources
from tiptoft.Fasta import Fasta
from tiptoft.Fastq import Fastq


class TipToft:
    def __init__(self, options):
        self.logger = logging.getLogger(__name__)
        self.plasmid_data = options.plasmid_data
        self.input_fastq = options.input_fastq
        self.kmer = options.kmer
        self.verbose = options.verbose
        self.min_fasta_hits = options.min_fasta_hits
        self.print_interval = options.print_interval
        self.output_file = options.output_file
        self.filtered_reads_file = options.filtered_reads_file
        self.max_gap = options.max_gap
        self.min_block_size = options.min_block_size
        self.margin = options.margin
        self.start_time = int(time.time())
        self.min_kmers_for_onex_pass = options.min_kmers_for_onex_pass
        self.min_perc_coverage = options.min_perc_coverage
        if options.no_hc_compression:
            self.homopolyer_compression = False
        else:
            self.homopolyer_compression = True
        self.max_kmer_count = options.max_kmer_count
        self.no_gene_filter = options.no_gene_filter

        if self.plasmid_data is None:
            self.plasmid_data = str(pkg_resources.resource_filename(
                __name__, 'data/plasmid_data.fa'))

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

        if self.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.ERROR)

    '''pass everything over to other classes'''

    def run(self):
        self.logger.info("Starting analysis")
        fasta = Fasta(self.logger,
                      self.plasmid_data,
                      self.kmer,
                      self.homopolyer_compression,
                      max_kmer_count=self.max_kmer_count)
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
