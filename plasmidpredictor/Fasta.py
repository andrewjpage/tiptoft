'''Read in a FASTA file and extract all the k-mers'''
import operator
from Bio import SeqIO
from plasmidpredictor.Kmers import Kmers

class Fasta:
	def __init__(self,logger, filename, k):
		self.logger = logger
		self.filename = filename
		self.k = k

		self.sequences_to_kmers = self.sequence_kmers()
		self.all_kmers = self.all_kmers_in_file()
		self.kmer_keys_set = set(self.all_kmers.keys())

	def sequence_kmers(self):
		seq_counter = 0
		
		kmer_to_sequences = {}
		for record in SeqIO.parse(self.filename, "fasta"):
			sequence_length  = len(record.seq)
			
			kmers = Kmers(str(record.seq), self.k)
			# We assume here that the sequence name is unique in the FASTA file
			kmer_to_sequences[record.id] = kmers.get_all_kmers_counter()
			
			seq_counter += 1
			
		return kmer_to_sequences
		
	def all_kmers_in_file(self):
		self.logger.info("Finding all k-mers in plasmid FASTA file")
		all_kmers = {}
		for seq_name, kmer_counts in self.sequences_to_kmers.items():
			for kmer, count in kmer_counts.items():
				if kmer in all_kmers:
					all_kmers[kmer] += 1
				else:
					all_kmers[kmer] = 1
		return all_kmers

		