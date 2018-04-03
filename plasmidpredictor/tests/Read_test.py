import unittest
import os
from plasmidpredictor.Read import Read

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','read')


class TestRead(unittest.TestCase):

	def test_initialise(self):
		read = Read()
		self.assertEqual(read.id, None)
		self.assertEqual(read.seq, None)
		self.assertEqual(read.qual, None)
		
	def test_one_read(self):
		f = open(os.path.join(data_dir, 'sample.fastq'))
		read = Read()
		read.get_next_from_file(f)
		self.assertEqual(read.id, 'read1')
		self.assertEqual(read.seq, 'AAAAAAAAAAAAGGGGGGGGGGGGGGAAAAAA')
		self.assertEqual(read.qual, '77777777777788888888888888777777')
		f.close()
		
	def test_subsequence(self):
		f = open(os.path.join(data_dir, 'sample.fastq'))
		read = Read()
		read.get_next_from_file(f)
		self.assertEqual(read.id, 'read1')
		self.assertEqual(read.seq, 'AAAAAAAAAAAAGGGGGGGGGGGGGGAAAAAA')
		self.assertEqual(read.qual, '77777777777788888888888888777777')
		
		sub_read  = read.subsequence(12,26)
		self.assertEqual(sub_read.id, 'read1_12_26')
		self.assertEqual(sub_read.seq, 'GGGGGGGGGGGGGG')
		self.assertEqual(sub_read.qual, '88888888888888')
		f.close()
		