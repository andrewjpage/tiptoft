import unittest
import os
import logging
from plasmidpredictor.Kmers import Kmers

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','kmers')

class TestKmers(unittest.TestCase):

	def test_four_kmers(self):
		k = Kmers('AAAAATTTTT',4)
		self.assertEqual(k.get_all_kmers_counter(), {'AAAA': 0, 'AAAT': 0, 'AATT': 0, 'ATTT': 0, 'TTTT': 0})
		
	def test_short_sequence(self):
		k = Kmers('A',10)
		self.assertEqual(k.get_all_kmers_counter(),{})
		