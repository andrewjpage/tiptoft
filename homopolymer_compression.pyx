# Run it with Cython
def homopolymer_compression_of_sequence(sequence):
	previous_base = ''
	compressed_sequence = []
	for base in sequence:
		if base != previous_base:
			previous_base = base
			compressed_sequence.append(base)
	return ''.join(compressed_sequence)
	