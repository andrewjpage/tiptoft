"""
Unit tests for InputTypes module.

This module tests the validation functions for command-line inputs.
"""
import unittest
import tempfile
import os
import argparse
from tiptoft.InputTypes import InputTypes


class TestInputTypes(unittest.TestCase):
    """Test cases for InputTypes validation functions."""

    def test_valid_fastq_file(self):
        """Test validation of existing FASTQ file."""
        with tempfile.NamedTemporaryFile(delete=False) as tmp:
            tmp_path = tmp.name
        try:
            result = InputTypes.is_fastq_file_valid(tmp_path)
            self.assertEqual(result, tmp_path)
        finally:
            os.unlink(tmp_path)

    def test_stdin_fastq_file(self):
        """Test that stdin dash character is accepted."""
        result = InputTypes.is_fastq_file_valid("-")
        self.assertEqual(result, "-")

    def test_invalid_fastq_file(self):
        """Test that non-existent FASTQ file raises error."""
        with self.assertRaises(argparse.ArgumentTypeError):
            InputTypes.is_fastq_file_valid("/nonexistent/file.fastq")

    def test_valid_plasmid_file(self):
        """Test validation of existing plasmid file."""
        with tempfile.NamedTemporaryFile(delete=False) as tmp:
            tmp_path = tmp.name
        try:
            result = InputTypes.is_plasmid_file_valid(tmp_path)
            self.assertEqual(result, tmp_path)
        finally:
            os.unlink(tmp_path)

    def test_invalid_plasmid_file(self):
        """Test that non-existent plasmid file raises error."""
        with self.assertRaises(argparse.ArgumentTypeError):
            InputTypes.is_plasmid_file_valid("/nonexistent/plasmid.fa")

    def test_valid_kmer_minimum(self):
        """Test validation of minimum valid k-mer size (7)."""
        result = InputTypes.is_kmer_valid("7")
        self.assertEqual(result, 7)

    def test_valid_kmer_maximum(self):
        """Test validation of maximum valid k-mer size (31)."""
        result = InputTypes.is_kmer_valid("31")
        self.assertEqual(result, 31)

    def test_valid_kmer_middle(self):
        """Test validation of middle-range k-mer size (13)."""
        result = InputTypes.is_kmer_valid("13")
        self.assertEqual(result, 13)

    def test_invalid_kmer_too_small(self):
        """Test that k-mer size below minimum (6) raises error."""
        with self.assertRaises(argparse.ArgumentTypeError):
            InputTypes.is_kmer_valid("6")

    def test_invalid_kmer_too_large(self):
        """Test that k-mer size above maximum (32) raises error."""
        with self.assertRaises(argparse.ArgumentTypeError):
            InputTypes.is_kmer_valid("32")

    def test_invalid_kmer_non_numeric(self):
        """Test that non-numeric k-mer value raises error."""
        with self.assertRaises(argparse.ArgumentTypeError):
            InputTypes.is_kmer_valid("abc")

    def test_invalid_kmer_negative(self):
        """Test that negative k-mer value raises error."""
        with self.assertRaises(argparse.ArgumentTypeError):
            InputTypes.is_kmer_valid("-5")
