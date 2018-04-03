import unittest
import os
from plasmidpredictor.Gene import Gene

class TestGene(unittest.TestCase):

	def test_full_coverage(self):
		g = Gene('ABC_123', 10, 0)
		self.assertTrue(g.is_full_coverage())
		self.assertEqual(str(g), 'ABC_123')

	def test_no_coverage(self):
		g = Gene('ABC_123', 0, 10)
		self.assertFalse(g.is_full_coverage())
		self.assertEqual(str(g), 'ABC_123*')
		
	def test_medium_coverage(self):
		g = Gene('ABC_123', 5, 5)
		self.assertFalse(g.is_full_coverage())
		self.assertEqual(str(g), 'ABC_123*')
