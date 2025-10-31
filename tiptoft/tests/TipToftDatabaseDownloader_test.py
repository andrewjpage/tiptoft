"""
Unit tests for TipToftDatabaseDownloader class.

This module tests the database downloader functionality.
"""
import unittest
from unittest.mock import Mock, patch
from tiptoft.TipToftDatabaseDownloader import TipToftDatabaseDownloader


class TestTipToftDatabaseDownloader(unittest.TestCase):
    """Test cases for TipToftDatabaseDownloader class."""

    def create_options(self, **kwargs):
        """Helper to create mock options object."""
        defaults = {
            'output_prefix': 'test_output',
            'verbose': False
        }
        defaults.update(kwargs)
        return Mock(**defaults)

    def test_initialization_with_defaults(self):
        """Test initialization with default parameters."""
        options = self.create_options()
        downloader = TipToftDatabaseDownloader(options)
        
        self.assertEqual(downloader.output_prefix, 'test_output')
        self.assertFalse(downloader.verbose)

    def test_initialization_verbose_mode(self):
        """Test initialization with verbose mode enabled."""
        options = self.create_options(verbose=True)
        downloader = TipToftDatabaseDownloader(options)
        
        import logging
        self.assertEqual(downloader.logger.level, logging.DEBUG)

    def test_initialization_non_verbose_mode(self):
        """Test initialization with verbose mode disabled."""
        options = self.create_options(verbose=False)
        downloader = TipToftDatabaseDownloader(options)
        
        import logging
        self.assertEqual(downloader.logger.level, logging.ERROR)

    def test_run_with_no_output_prefix(self):
        """Test that run without output prefix logs error."""
        options = self.create_options(output_prefix=None)
        downloader = TipToftDatabaseDownloader(options)
        
        # Should not raise an exception, just log error
        downloader.run()

    @patch('tiptoft.TipToftDatabaseDownloader.RefGenesGetter')
    def test_run_with_output_prefix(self, mock_refgenes_class):
        """Test that run with output prefix calls RefGenesGetter."""
        mock_refgenes = Mock()
        mock_refgenes_class.return_value = mock_refgenes
        
        options = self.create_options(output_prefix='test_db')
        downloader = TipToftDatabaseDownloader(options)
        downloader.run()
        
        # Verify RefGenesGetter was instantiated
        mock_refgenes_class.assert_called_once_with(verbose=False)
        # Verify run was called with output prefix
        mock_refgenes.run.assert_called_once_with('test_db')
