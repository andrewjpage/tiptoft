'''Read in a FASTQ file and identify matching plasmids'''
from plasmidpredictor.Kmers import Kmers
from plasmidpredictor.Read import Read
from plasmidpredictor.Gene import Gene
from plasmidpredictor.Blocks import Blocks
import subprocess
import os
import numpy
import time
import sys

class Error (Exception): pass

class Fastq:
	def __init__(self,logger, filename, k, fasta_kmers, min_fasta_hits, print_interval, output_file, filtered_reads_file, fasta_obj, max_gap = 4, min_block_size = 150, margin = 100, start_time = 0, min_kmers_for_onex_pass = 10, max_kmers = 5, min_perc_coverage = 95 ):
		self.logger = logger
		self.filename = filename
		self.k = k
		self.fasta_kmers = fasta_kmers
		self.min_fasta_hits = min_fasta_hits
		self.print_interval = print_interval
		self.output_file = output_file
		self.filtered_reads_file = filtered_reads_file
		self.max_gap = max_gap # multiples of the kmer
		self.min_block_size = min_block_size
		self.margin = margin
		self.start_time = start_time
		self.min_kmers_for_onex_pass = min_kmers_for_onex_pass
		self.max_kmers = max_kmers
		self.fasta_obj = fasta_obj
		self.min_perc_coverage = min_perc_coverage

	def read_filter_and_map(self):
		counter = 0 
		match_counter = 0
		
		self.logger.info("Reading in FASTQ file")
		fh = self.open_file_read()
		read = Read()

		while read.get_next_from_file(fh):
			counter += 1
			if counter % self.print_interval == 0:
				self.full_gene_coverage(counter)

			if not self.map_read(read):
				self.map_read(read.reverse_read())
				
		self.full_gene_coverage(counter)
								
		return self
		
	def map_read(self, read):
		if self.does_read_contain_quick_pass_kmers(read.seq):
			self.map_kmers_to_read(read.seq, read)
			return True
		else:
			return False
		
	def does_read_contain_quick_pass_kmers(self, sequence):
		self.logger.info("Perform quick pass k-mer check on read")
		seq_length = len(sequence)
		if seq_length < self.min_block_size:
			self.logger.info("Read below minimum size")
			return False
		
		kmers_obj = Kmers(sequence, self.k)
		read_onex_kmers = kmers_obj.get_one_x_coverage_of_kmers()
		
		fasta_kmers_seqs = self.fasta_kmers.keys()
		hit_counter = 0
		for r in read_onex_kmers:
			if r in fasta_kmers_seqs:
				hit_counter += 1
		
				if hit_counter > self.min_kmers_for_onex_pass:
					return True

		return False
		
		
	def put_kmers_in_read_bins(self, seq_length, end, fasta_kmers, read_kmers):
		self.logger.info("Put k-mers in read bins")
		sequence_hits = numpy.zeros(int(seq_length/self.k)+1, dtype=int)
		hit_counter = 0
		
		for read_kmer, read_kmer_hit in read_kmers.items():
			if read_kmer in fasta_kmers:
				for coordinate in read_kmer_hit.coordinates:
					hit_counter += 1
					sequence_hits[int(coordinate/self.k)] += 1
		return sequence_hits, hit_counter

		
	def map_kmers_to_read(self, sequence, read):	
		self.logger.info("Map k-mers to read")	

		seq_length = len(sequence)
		end = seq_length - self.k
		
		kmers_obj = Kmers(sequence, self.k)
		read_kmers = kmers_obj.get_all_kmers_filtered(self.max_kmers)
		is_read_matching = False
		
		sequence_hits, hit_counter = self.put_kmers_in_read_bins( seq_length, end, self.fasta_kmers, read_kmers)
			
		blocks_obj = Blocks(self.k, self.min_block_size, self.max_gap, self.margin)
		block_start, block_end = blocks_obj.find_largest_block(sequence_hits)
			
		block_start = blocks_obj.adjust_block_start(block_start)
		block_end = blocks_obj.adjust_block_end(block_end, seq_length)

		block_kmers = self.create_kmers_for_block(block_start, block_end, sequence)
		self.apply_kmers_to_genes(self.fasta_obj,block_kmers)
		is_read_matching = True
		
		if self.filtered_reads_file:
			self.append_subread_to_fastq_file(read, block_start, block_end)
				
		return is_read_matching
			
	def append_subread_to_fastq_file(self, read, block_start, block_end):
			with open(self.filtered_reads_file, 'a+') as output_fh:
				output_fh.write(str(read.subsequence(block_start, block_end)))
			
	def create_kmers_for_block(self, block_start, block_end, sequence):
		if block_end ==  0:
			return {}
		
		block_seq = sequence[block_start:block_end]
			
		return Kmers(block_seq, self.k).get_all_kmers(self.max_kmers)
	
	def apply_kmers_to_genes(self,fasta_obj,hit_kmers):
		for (gene_name, kmers_dict) in fasta_obj.sequences_to_kmers.items():
			for kmer in kmers_dict.keys():
				if kmer in hit_kmers:
					fasta_obj.sequences_to_kmers[gene_name][kmer] += 1 		
		
	def full_gene_coverage(self, counter):
		self.logger.info("Check the coverage of a sequence")
		alleles = []
		
		for (gene_name, kmers_dict) in self.fasta_obj.sequences_to_kmers.items():
			kv = kmers_dict.values()
			
			sum_kmer_coverage = sum(kv)

			kv_coverage = [x for x in kv if x >= 1]
			kv_zero = [x for x in kv if x == 0]
			kl  = len(kv_coverage)
			kz = len(kv_zero)
			if kl > kz:
				alleles.append(Gene(gene_name, kl, kz))
				
		self.print_out_alleles(alleles)
		print("****")
		return alleles
		
	def print_out_alleles(self,alleles):
		for g in alleles:
			if g.percentage_coverage() >= self.min_perc_coverage:
				print(g)
		
	# Derived from https://github.com/sanger-pathogens/Fastaq
	# Author: Martin Hunt	
	def open_file_read(self):
		if self.filename == '-':
			f = sys.stdin
		elif self.filename.endswith('.gz'):
			# first check that the file is OK according to gunzip
			retcode = subprocess.call('gunzip -t ' + self.filename, shell=True)
			
			# now open the file
			f = os.popen('gunzip -c ' + self.filename)
		else:
			try:
				f = open(self.filename)
			except:
				raise Error("Error opening for reading file '" + self.filename + "'")
		
		return f
	
