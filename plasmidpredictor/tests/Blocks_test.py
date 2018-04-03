import unittest
import os
from plasmidpredictor.Blocks import Blocks

class TestBlocks(unittest.TestCase):

	def test_two_blocks(self):
		b = Blocks(7, 7,2, 5)
		hits = [0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1]
		self.assertEqual(b.find_all_blocks(hits), [[3, 14], [27, 36]])
		
	def test_merging_blocks(self):
		b = Blocks(7, 7,3, 5)
		blocks = [[4,8],[10,14],[20,24],[26,30]]
		self.assertEqual(b.merge_blocks(blocks), [[4, 14], [4, 14], [20, 30], [20, 30]])
		
	def test_largest_block(self):
		b = Blocks(7, 7,2, 5)
		hits = [0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1]
		self.assertEqual(b.find_largest_block(hits), (3, 14))	
		
	def test_adjust_block_start(self):
		b = Blocks(7, 0, 0, 10)
		self.assertEqual(b.adjust_block_start(10),60)
		self.assertEqual(b.adjust_block_start(1),0)
	
	def test_adjust_block_end(self):
		b = Blocks(7, 0, 0, 10)
		self.assertEqual(b.adjust_block_end(10, 90),80)
		self.assertEqual(b.adjust_block_end(10, 80),80)
		self.assertEqual(b.adjust_block_end(10, 50),50)
