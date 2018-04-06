'''Given a string of nucleotides and k, return all kmers'''
import re

class KmerHit():
	def __init__(self,count, coord):
		self.count = 0
		self.coordinates = [coord]

class Kmers:
	def __init__(self, sequence, k, homopolyer_compression):
		self.sequence = self.homopolymer_compression_of_sequence(sequence) if homopolyer_compression else sequence
		self.k = k
		self.end = len(self.sequence) - self.k + 1
		
	def homopolymer_compression_of_sequence(self,sequence):
		p = re.compile(r'(.)\1*')
		compressed_sequence = p.sub(r'\1', sequence)
		return compressed_sequence
	
	def get_all_kmers_counter(self, max_kmer_count = 10):
		kmers = self.get_all_kmers_filtered(max_kmer_count)
		return { x:0 for x in kmers.keys()}
	
	def get_all_kmers(self, max_kmer_count = 10):
		return self.get_all_kmers_filtered(max_kmer_count)
	     
	def get_all_kmers_filtered(self, max_kmer_count = 10):
		kmers = {}

		kmer_sequences = [ self.sequence[i:i+self.k] for i in range(0,self.end)]
		
		for i,k in enumerate(kmer_sequences):
		    if k in kmers:
		        kmers[k].count += 1
		        kmers[k].coordinates.append(i)
		    else:
		        kmers[k] = KmerHit(1,i)
		     
		filtered_kmers = {k:v for k,v in kmers.items() if int(v.count) <= max_kmer_count }
		return filtered_kmers
	     
	def get_one_x_coverage_of_kmers(self):
		kmers = [ self.sequence[i:i+self.k] for i in range(0,self.end, self.k)]
		return kmers
		