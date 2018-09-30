import unittest
import os
import logging
import filecmp
from tiptoft.Fasta import Fasta
from tiptoft.Fastq import Fastq
from tiptoft.Fastq import Gene

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
	
	def test_gzipped_input(self):
		logger = logging.getLogger(__name__)
		logger.setLevel(logging.ERROR)
		fasta = Fasta(logger, os.path.join(data_dir,'plasmid_data.fa'),4, True)
		
		fastq = Fastq(logger, os.path.join(data_dir,'query_gz.fastq.gz'), 4 , fasta.all_kmers_in_file(), 1, 50, None, None, fasta, True)
		self.assertTrue(fastq.read_filter_and_map())
		
	def test_writting_to_output_file(self):
		logger = logging.getLogger(__name__)
		logger.setLevel(logging.ERROR)
		fasta = Fasta(logger, os.path.join(data_dir,'plasmid_data.fa'),4, True)
		
		fastq = Fastq(logger, os.path.join(data_dir,'query_gz.fastq.gz'), 4 , fasta.all_kmers_in_file(), 1, 50, 'outputfile', None, fasta, True)
		fastq.read_filter_and_map()
		self.assertTrue(os.path.exists('outputfile'))
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_outputfile'), 'outputfile'))
		os.remove('outputfile')
		
		
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


	def test_filtering_alleles_one_complete(self):
		logger = logging.getLogger(__name__)
		logger.setLevel(logging.ERROR)
		fastq = Fastq(logger, os.path.join(data_dir,'query.fastq'), 11 , None, 1, 50, None, None, None, True)

		input_alleles = [ Gene('rep7.1_repC(Cassette)_AB037671', 9, 1), Gene('rep7.5_CDS1(pKC5b)_AF378372', 8, 2), Gene('rep7.6_ORF(pKH1)_SAU38656', 10, 0), Gene('repUS14.1_repA(VRSAp)_AP003367', 10, 0)]
		expected_allele_names = ['rep7.6','repUS14.1']
		filtered_alleles = fastq.filter_contained_alleles(input_alleles)
		self.assertEqual(expected_allele_names, sorted(list(map(lambda x: x.short_name(), filtered_alleles))))

	def test_filtering_alleles_all_partial(self):
		logger = logging.getLogger(__name__)
		logger.setLevel(logging.ERROR)
		fastq = Fastq(logger, os.path.join(data_dir,'query.fastq'), 11 , None, 1, 50, None, None, None, True)

		input_alleles = [ Gene('rep7.1_repC(Cassette)_AB037671', 9, 1), Gene('rep7.5_CDS1(pKC5b)_AF378372', 8, 2), Gene('rep7.6_ORF(pKH1)_SAU38656', 7, 3), Gene('repUS14.1_repA(VRSAp)_AP003367', 10, 0)]
		expected_allele_names = ['rep7.1','repUS14.1']
		filtered_alleles = fastq.filter_contained_alleles(input_alleles)
		self.assertEqual(expected_allele_names, sorted(list(map(lambda x: x.short_name(), filtered_alleles))))

	def test_filtering_alleles_partial_equal_values(self):
		logger = logging.getLogger(__name__)
		logger.setLevel(logging.ERROR)
		fastq = Fastq(logger, os.path.join(data_dir,'query.fastq'), 11 , None, 1, 50, None, None, None, True)

		input_alleles = [ Gene('rep7.1_repC(Cassette)_AB037671', 9, 1), Gene('rep7.5_CDS1(pKC5b)_AF378372', 9, 1), Gene('rep7.6_ORF(pKH1)_SAU38656', 9, 1), Gene('repUS14.1_repA(VRSAp)_AP003367', 10, 0)]
		expected_allele_names = ['rep7.1','repUS14.1']
		filtered_alleles = fastq.filter_contained_alleles(input_alleles)
		self.assertEqual(expected_allele_names, sorted(list(map(lambda x: x.short_name(), filtered_alleles))))

	def test_filtering_alleles_all_complete(self):
		logger = logging.getLogger(__name__)
		logger.setLevel(logging.ERROR)
		fastq = Fastq(logger, os.path.join(data_dir,'query.fastq'), 11 , None, 1, 50, None, None, None, True)

		input_alleles = [ Gene('rep7.1_repC(Cassette)_AB037671', 10, 0), Gene('rep7.5_CDS1(pKC5b)_AF378372', 10, 0), Gene('rep7.6_ORF(pKH1)_SAU38656', 10, 0), Gene('repUS14.1_repA(VRSAp)_AP003367', 10, 0)]
		expected_allele_names = ['rep7.1', 'rep7.5', 'rep7.6', 'repUS14.1']
		filtered_alleles = fastq.filter_contained_alleles(input_alleles)
		self.assertEqual(expected_allele_names, sorted(list(map(lambda x: x.short_name(), filtered_alleles))))

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


