'''Read in a FASTA file and extract all the k-mers'''
from Bio import SeqIO
from tiptoft.Kmers import Kmers


class Fasta:
    def __init__(
            self,
            logger,
            filename,
            k,
            homopolyer_compression,
            max_kmer_count=10):
        self.logger = logger
        self.filename = filename
        self.k = k
        self.homopolyer_compression = homopolyer_compression
        self.max_kmer_count = max_kmer_count

        self.sequences_to_kmers = self.sequence_kmers('get_all_kmers_counter')
        self.sequences_to_kmers_count =\
            self.sequence_kmers('get_all_kmers_freq')
        self.all_kmers = self.all_kmers_in_file()
        self.kmers_to_genes = self.all_kmers_to_seq_in_file()
        self.kmer_keys_set = set(self.all_kmers.keys())

    '''Count the kmers in a sequence'''

    def sequence_kmers(self, kmer_action='get_all_kmers_counter'):
        seq_counter = 0
        kmer_to_sequences = {}
        for record in SeqIO.parse(self.filename, "fasta"):
            kmers = Kmers(str(record.seq), self.k, self.homopolyer_compression)
            # We assume here that the sequence name is unique in the FASTA file
            kmer_to_sequences[record.id] = getattr(kmers, kmer_action)(
                max_kmer_count=self.max_kmer_count)

            seq_counter += 1

        return kmer_to_sequences

    '''create a dictionary of kmers to underlying genes'''

    def all_kmers_to_seq_in_file(self):
        kmers_to_genes = {}
        for seq_name, kmer_counts in self.sequences_to_kmers.items():
            for kmer, count in kmer_counts.items():

                if kmer not in kmers_to_genes:
                    kmers_to_genes[kmer] = []
                kmers_to_genes[kmer].append(seq_name)

        return kmers_to_genes

    '''given a fasta file, extract all kmers'''

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
