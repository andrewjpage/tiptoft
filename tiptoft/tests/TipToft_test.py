"""
Unit tests for TipToft main driver class.

This module tests the TipToft class which coordinates the plasmid detection workflow.
"""
import unittest
import tempfile
import os
from unittest.mock import Mock, MagicMock, patch
from tiptoft.TipToft import TipToft


class TestTipToft(unittest.TestCase):
    """Test cases for TipToft driver class."""

    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary FASTQ file
        self.temp_fastq = tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False)
        self.temp_fastq.write('@read1\nATGC\n+\nIIII\n')
        self.temp_fastq.close()

    def tearDown(self):
        """Clean up test fixtures."""
        if os.path.exists(self.temp_fastq.name):
            os.unlink(self.temp_fastq.name)

    def create_options(self, **kwargs):
        """Helper to create mock options object."""
        defaults = {
            'plasmid_data': None,
            'input_fastq': self.temp_fastq.name,
            'kmer': 13,
            'verbose': False,
            'min_fasta_hits': 10,
            'print_interval': None,
            'output_file': None,
            'filtered_reads_file': None,
            'max_gap': 3,
            'min_block_size': 130,
            'margin': 10,
            'min_kmers_for_onex_pass': 10,
            'min_perc_coverage': 85,
            'no_hc_compression': False,
            'max_kmer_count': 10,
            'no_gene_filter': False
        }
        defaults.update(kwargs)
        return Mock(**defaults)

    def test_initialization_with_defaults(self):
        """Test TipToft initialization with default parameters."""
        options = self.create_options()
        tiptoft = TipToft(options)
        
        self.assertEqual(tiptoft.kmer, 13)
        self.assertEqual(tiptoft.min_fasta_hits, 10)
        self.assertTrue(tiptoft.homopolyer_compression)
        self.assertIsNotNone(tiptoft.plasmid_data)

    def test_initialization_with_hc_disabled(self):
        """Test initialization with homopolymer compression disabled."""
        options = self.create_options(no_hc_compression=True)
        tiptoft = TipToft(options)
        
        self.assertFalse(tiptoft.homopolyer_compression)

    def test_initialization_with_custom_plasmid_data(self):
        """Test initialization with custom plasmid data file."""
        temp_plasmid = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False)
        temp_plasmid.write('>test\nATGC\n')
        temp_plasmid.close()
        
        try:
            options = self.create_options(plasmid_data=temp_plasmid.name)
            tiptoft = TipToft(options)
            self.assertEqual(tiptoft.plasmid_data, temp_plasmid.name)
        finally:
            os.unlink(temp_plasmid.name)

    def test_initialization_with_existing_filtered_reads_file(self):
        """Test that existing filtered reads file causes error."""
        temp_filtered = tempfile.NamedTemporaryFile(delete=False)
        temp_filtered.close()
        
        try:
            options = self.create_options(filtered_reads_file=temp_filtered.name)
            with self.assertRaises(SystemExit):
                TipToft(options)
        finally:
            os.unlink(temp_filtered.name)

    def test_verbose_mode_sets_debug_logging(self):
        """Test that verbose mode enables debug logging."""
        options = self.create_options(verbose=True)
        tiptoft = TipToft(options)
        
        import logging
        self.assertEqual(tiptoft.logger.level, logging.DEBUG)

    def test_non_verbose_mode_sets_error_logging(self):
        """Test that non-verbose mode sets error logging."""
        options = self.create_options(verbose=False)
        tiptoft = TipToft(options)
        
        import logging
        self.assertEqual(tiptoft.logger.level, logging.ERROR)
