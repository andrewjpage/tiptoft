'''Given an array of integers, identify blocks which are likely to be genes of interest'''

class Blocks:
	def __init__(self, k, min_block_size,max_gap, margin):
		self.k = k
		self.min_block_size = min_block_size
		self.max_gap = max_gap # multiples of the kmer
		self.margin = margin
		
	'''An array of pairs, with start and end coordinate of each block is passed in as input. If these blocks are close together, merge them.'''
	def merge_blocks(self, blocks):
		for i in range(0,len(blocks)-1 ):
			if blocks[i][1] + self.max_gap > blocks[i+1][0]:
				if blocks[i][1]  < blocks[i+1][1]:
					blocks[i][1] = 	blocks[i+1][1]
				blocks[i+1][0] = blocks[i][0]
				blocks[i+1][1] = blocks[i][1]
		return blocks
		
	'''Given a set of blocks (start and end coordinates), return the largest block coordinates'''
	def find_largest_block(self, sequence_hits):
		blocks = self.find_all_blocks(sequence_hits)
		merged_blocks = self.merge_blocks(blocks)
		
		largest_block = 0
		largest_block_index = 0
		
		for i, block in enumerate(merged_blocks):
			block_size = block[1] - block[0]
			if block_size > largest_block:
				largest_block_index = i
				largest_block = block_size
				
		if largest_block < (self.min_block_size/self.k):
			return 0,0
				
		return merged_blocks[largest_block_index][0], merged_blocks[largest_block_index][1]
		
	'''Given an array with kmer matches, find the coordinates of each contiguous segment'''
	def find_all_blocks(self, sequence_hits):
		blocks = []
		in_block = False
		current_block_start = 0
		for i,val_count in enumerate(sequence_hits):
			
			if not in_block and val_count > 0:
				in_block = True
				current_block_start = i
			elif in_block and val_count == 0: 
				in_block = False
				blocks.append([current_block_start, i])
				
		if in_block:
			blocks.append([current_block_start, len(sequence_hits)])
				
		return blocks	
		
	'''Rescale the blocks, add a flanking margin and make sure the start is in bounds'''
	def adjust_block_start(self,block_start):
		block_start *= self.k
		if block_start - self.margin < 0:
			block_start = 0
		else:
			block_start -= self.margin
		return block_start
		
	'''Rescale the blocks, add a flanking margin and make sure the end is in bounds'''
	def adjust_block_end(self,block_end,seq_length):
		block_end *= self.k
		if block_end + self.margin > seq_length:
			block_end = seq_length
		else:
			block_end +=  self.margin
		return block_end
		