import unittest
import os
import logging
import filecmp
from plasmidpredictor.Fasta import Fasta
from plasmidpredictor.Fastq import Fastq

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','fastq')

#import cProfile, pstats, io

class TestFastq(unittest.TestCase):

	def test_four_kmers(self):
		logger = logging.getLogger(__name__)
		logger.setLevel(logging.ERROR)
		fasta = Fasta(logger, os.path.join(data_dir,'plasmid_data.fa'),4, True)
		
		fastq = Fastq(logger, os.path.join(data_dir,'query.fastq'), 4 , fasta.all_kmers_in_file(), 1, 50, None, None, fasta, True)
				
		self.assertTrue(fastq.read_filter_and_map())
		
	def test_with_nonmatching_read(self):
		logger = logging.getLogger(__name__)
		logger.setLevel(logging.ERROR)
		fasta = Fasta(logger, os.path.join(data_dir,'plasmid_data.fa'),4, True)
		
		fastq = Fastq(logger, os.path.join(data_dir,'query.fastq'), 4 , fasta.all_kmers_in_file(), 1, 50, None, None, fasta, True)
		
		self.assertFalse(fastq.does_read_contain_quick_pass_kmers("AAAAAAAAAAAAAAAA"))
		
	def test_with_matching_read(self):
		logger = logging.getLogger(__name__)
		logger.setLevel(logging.ERROR)
		fasta = Fasta(logger, os.path.join(data_dir,'plasmid_data.fa'),11, True)
	
		fastq = Fastq(logger, os.path.join(data_dir,'query.fastq'), 11 , fasta.all_kmers_in_file(), 1, 50, None, None, fasta, True)
	
		self.assertTrue(fastq.does_read_contain_quick_pass_kmers("ATCAATACCTTCTTTATTGATTTTGATATTCACACGGCAAAAGAAACTATTTCAGCAAGCGATATTTTAACAACCGCTATTGATTTAGGTTTTATGCCTACTATGATTATCAAATCTGATAAAGGTTATCAAGCATATTTTGTTTTAGAAACGCCAGTCTATGTGACTTCAAAATCAGAATTTAAATCTGTCAAAGCAGCCAAAATAATTTCGCAAAATATCCGAGAATATTTTGGAAAGTCTTTGCCAGTTGATCTAACGTGTAATCATTTTGGTATTGCTCGCATACCAAGAACGGACAATGTAGAATTTTTTGATCCTAATTACCGTTATTCTTTCAAAGAATGGCAAGATTGGTCTTTCAAACAAACAGATAATAAGGGCTTTACTCGTTCAAGTCTAACGGTTTTAAGCGGTACAGAAGGCAAAAAACAAGTAGATGAACCCTGGTTTAATCTCTTATTGCACGAAACGAAATTTTCAGGAGAAAAGGGTTTAATAGGGCGTAATAACGTCATGTTTACCCTCTCTTTAGCCTACTTTAGTTCAGGCTATTCAATCGAAACGTGCGAATATAATATGTTTGAGTTTAATAATCGATTAGATCAACCCTTAGAAGAAAAAGAAGTAATCAAAATTGTTAGAAGTGCCTATTCAGAAAACTATCAAGGGGCTAATAGGGAATACATTACCATTCTTTGCAAAGCTTGGGTATCAAGTGATTTAACCAGTAAAGATTTATTTGTCCGTCAAGGGTGGTTTAAATTCAAGAAAAAAAGAAGCGAACGTCAACGTGTTCATTTGTCAGAATGGAAAGAAGATTTAATGGCTTATATTAGCGAAAAAAGCGATGTATACAAGCCTTATTTAGTGACGACCAAAAAAGAGATTAGAGAAGTG"))


#	def test_stap_aureus_pacbio(self):
#		logger = logging.getLogger(__name__)
#		kmer = 11
#		
#		pr = cProfile.Profile()
#		pr.enable()      
#		
#		fastas = Fastas(logger, databases+'/Staphylococcus_aureus/',kmer, False)
#		mlst_profile = MlstProfile(databases+'/Staphylococcus_aureus/profile.txt')
#		
#		fastq = Fastq(logger, os.path.join(data_dir,'NCTC11150.fastq.gz'), kmer , fastas.get_fastas_to_kmers(), 10,  mlst_profile, 50, None, None)
#		self.assertTrue(fastq.read_filter_and_map())
#
#		pr.disable()
#		s = io.StringIO()
#		sortby = 'cumulative'
#		ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
#		ps.print_stats()
#		print(s.getvalue())


