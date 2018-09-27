'''Given a string of nucleotides and k, return all kmers'''

# The homopolymer compression is quite intensive so run it with cython
import pyximport
pyximport.install()

from homopolymer_compression import homopolymer_compression_of_sequence

class Kmers:
	def __init__(self, sequence, k, homopolyer_compression):
		self.sequence = homopolymer_compression_of_sequence(sequence) if homopolyer_compression else sequence
		self.k = k
		self.end = len(self.sequence) - self.k + 1
		
	def get_all_kmers_counter(self, max_kmer_count = 10):
		kmers = self.get_all_kmers_filtered(max_kmer_count)
		return { x:0 for x in kmers.keys()}
		
	def get_all_kmers_freq(self, max_kmer_count = 10):
		kmers = self.get_all_kmers_filtered(max_kmer_count)
		return { k:len(v) for k,v in kmers.items()}
	
	def get_all_kmers(self, max_kmer_count = 10):
		return self.get_all_kmers_filtered(max_kmer_count)
	     
	def get_all_kmers_filtered(self, max_kmer_count = 10):
		kmers = {}

		kmer_sequences = [ self.sequence[i:i+self.k] for i in range(0,self.end)]
		
		for i,k in enumerate(kmer_sequences):
		    if k in kmers:
		        kmers[k].append(i)
		    else:
		        kmers[k] = [i]
		     
		filtered_kmers = {k:v for k,v in kmers.items() if len(v) <= max_kmer_count }
		return filtered_kmers
	     
	def get_one_x_coverage_of_kmers(self):
		return [ self.sequence[i:i+self.k] for i in range(0,self.end, self.k)]
		
