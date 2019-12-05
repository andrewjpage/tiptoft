'''Read in a FASTQ file and identify matching plasmids'''
from tiptoft.Kmers import Kmers
from tiptoft.Read import Read
from tiptoft.Gene import Gene
from tiptoft.Blocks import Blocks
import subprocess
import os
import numpy
import sys


class Error (Exception):
    pass


class Fastq:
    def __init__(
            self,
            logger,
            filename,
            k,
            fasta_kmers,
            min_fasta_hits,
            print_interval,
            output_file,
            filtered_reads_file,
            fasta_obj,
            homopolyer_compression,
            max_gap=4,
            min_block_size=150,
            margin=100,
            start_time=0,
            min_kmers_for_onex_pass=10,
            min_perc_coverage=95,
            max_kmer_count=10,
            no_gene_filter=False):
        self.logger = logger
        self.filename = filename
        self.k = k
        self.fasta_kmers = fasta_kmers
        self.min_fasta_hits = min_fasta_hits
        self.print_interval = print_interval
        self.output_file = output_file
        self.filtered_reads_file = filtered_reads_file
        self.max_gap = max_gap  # multiples of the kmer
        self.min_block_size = min_block_size
        self.margin = margin
        self.start_time = start_time
        self.min_kmers_for_onex_pass = min_kmers_for_onex_pass
        self.fasta_obj = fasta_obj
        self.min_perc_coverage = min_perc_coverage
        self.genes_with_100_percent = {}
        self.homopolyer_compression = homopolyer_compression
        self.max_kmer_count = max_kmer_count
        self.no_gene_filter = no_gene_filter

    '''Read in a whole FASTQ file, do a quick pass filter, then map more
         sensitively'''

    def read_filter_and_map(self):
        counter = 0
        match_counter = 0
        f = 0
        r = 0

        self.logger.info("Reading in FASTQ file")
        self.print_out_header()
        fh = self.open_file_read()
        read = Read()

        while read.get_next_from_file(fh):
            counter += 1

            if(self.print_interval is not None and
                    counter % self.print_interval == 0):
                self.full_gene_coverage(counter)

            if self.map_read(read):
                match_counter += 1
                f += 1
            elif self.map_read(read.reverse_read()):
                match_counter += 1
                r += 1
        fh.close()

        alleles = self.full_gene_coverage(counter)
        self.logger.info("Number of reads: "+str(counter))
        self.logger.info("Number of matching reads: " +
                         str(match_counter)+"\t"+str(f)+"\t"+str(r))
		
        print("\t".join([self.filename, str(len(alleles)), str(match_counter), str(counter), str(match_counter*100/counter)]))

        return self

    '''Take a single read, do a quick kmer check to see if it matches any gene,
         then map more sensitively'''

    def map_read(self, read):
        candidate_gene_names = self.does_read_contain_quick_pass_kmers(
            read.seq)
        if len(candidate_gene_names) > 0:
            self.logger.info("Read passes 1X check")
            return self.map_kmers_to_read(read.seq, read, candidate_gene_names)
        else:
            return False

    '''Taking kmers at 1x, see do any match the gene kmers'''

    def does_read_contain_quick_pass_kmers(self, sequence):
        self.logger.info("Perform quick pass k-mer check on read")
        seq_length = len(sequence)
        if seq_length < self.min_block_size:
            self.logger.info("Read below minimum size")
            return {}

        kmers_obj = Kmers(sequence, self.k, self.homopolyer_compression)
        read_onex_kmers = kmers_obj.get_one_x_coverage_of_kmers()

        intersect_read_fasta_kmers = self.fasta_obj.kmer_keys_set & set(
            read_onex_kmers)

        if len(intersect_read_fasta_kmers) > self.min_kmers_for_onex_pass:
            gene_names = self.genes_containing_first_pass_kmers(
                self.fasta_obj, intersect_read_fasta_kmers)

            mk1x = self.min_kmers_for_onex_pass
            candidate_gene_names = {
                k: v for k, v in gene_names.items() if v > mk1x}

            return candidate_gene_names

        return {}

    '''place kmers into read bins'''

    def put_kmers_in_read_bins(self, seq_length, end, fasta_kmers, read_kmers):
        self.logger.info("Put k-mers in read bins")
        sequence_hits = numpy.zeros(int(seq_length/self.k)+1, dtype=int)
        hit_counter = 0

        hit_kmers = {}
        for read_kmer, read_kmer_hit in read_kmers.items():
            if read_kmer in fasta_kmers:
                for coordinate in read_kmer_hit:
                    hit_counter += 1
                    sequence_hits[int(coordinate/self.k)] += 1
                    hit_kmers[read_kmer] = read_kmer_hit
        return sequence_hits, hit_counter, hit_kmers

    '''do a fine grained mapping of kmers to a read'''

    def map_kmers_to_read(self, sequence, read, candidate_gene_names):
        self.logger.info("Map k-mers to read")

        seq_length = len(sequence)
        end = seq_length - self.k

        kmers_obj = Kmers(sequence, self.k, self.homopolyer_compression)
        read_kmers = kmers_obj.get_all_kmers(
            max_kmer_count=self.max_kmer_count)
        is_read_matching = False

        sequence_hits, hit_counter, read_kmer_hits =\
            self.put_kmers_in_read_bins(
                seq_length,
                end,
                self.fasta_kmers,
                read_kmers)

        blocks_obj = Blocks(self.k, self.min_block_size,
                            self.max_gap, self.margin)
        block_start, block_end = blocks_obj.find_largest_block(sequence_hits)

        block_start = blocks_obj.adjust_block_start(block_start)
        block_end = blocks_obj.adjust_block_end(block_end, seq_length)

        block_kmers = self.create_kmers_for_block(
            block_start, block_end, read_kmer_hits)
        is_read_matching = self.apply_kmers_to_genes(
            self.fasta_obj, block_kmers, candidate_gene_names)

        if self.filtered_reads_file:
            self.append_read_to_fastq_file(read, block_start, block_end)

        return is_read_matching

    '''optional method to output only reads matching the genes allowing for
         offline assembly'''

    def append_read_to_fastq_file(self, read, block_start, block_end):
        with open(self.filtered_reads_file, 'a+') as output_fh:
            output_fh.write(str(read))

    '''given a putative block where a gene may match, get all the kmers'''

    def create_kmers_for_block(self, block_start, block_end, read_kmer_hits):
        if block_end == 0:
            return {}

        def get_one_x_coverage_of_kmers(self, sequence, k, end):
            return [sequence[i:i+k] for i in range(0, end, k)]
        block_kmers = {}

        for read_kmer, read_kmer_hit in read_kmer_hits.items():
            found_match_in_block = False
            for coordinate in read_kmer_hit:
                if coordinate >= block_start and coordinate <= block_end:
                    found_match_in_block = True
                    continue
            if found_match_in_block:
                block_kmers[read_kmer] = 1

        return block_kmers

    '''Get a list of genes which may be in the reads'''

    def genes_containing_first_pass_kmers(self, fasta_obj, first_pass_kmers):
        genes = {}
        for current_kmer in first_pass_kmers:
            if current_kmer in fasta_obj.kmers_to_genes:
                for gene_name in fasta_obj.kmers_to_genes[current_kmer]:
                    if gene_name in self.genes_with_100_percent:
                        continue

                    if gene_name in genes:
                        genes[gene_name] += 1
                    else:
                        genes[gene_name] = 1
        return genes

    '''given some kmer hits, apply the kmers to the genes so theres a count'''

    def apply_kmers_to_genes(self, fasta_obj, hit_kmers, gene_names):
        num_genes_applied = 0
        min_kmers = self.min_kmers_for_onex_pass * self.k
        hit_kmers_set = set(hit_kmers)

        for gene_name in gene_names.keys():
            if gene_names[gene_name] < self.min_kmers_for_onex_pass:
                continue
            kmers_dict = fasta_obj.sequences_to_kmers[gene_name]
            intersection_hit_keys = set(kmers_dict) & hit_kmers_set

            if len(intersection_hit_keys) > min_kmers:
                num_genes_applied += 1
                for kmer in intersection_hit_keys:

                    if fasta_obj.sequences_to_kmers_count[gene_name][kmer] > 0:
                        fasta_obj.sequences_to_kmers[gene_name][kmer] += 1 / \
                            fasta_obj.sequences_to_kmers_count[gene_name][kmer]

        if num_genes_applied > 0:
            return True
        else:
            return False

    '''calculate the coverage of the genes'''

    def full_gene_coverage(self, counter):
        self.logger.info("Check the coverage of a sequence")
        alleles = []
        s2k = self.fasta_obj.sequences_to_kmers
        for (gene_name, kmers_dict) in s2k.items():
            kv = kmers_dict.values()
            kv_total_length = len(kmers_dict)
            kz = len([x for x in kv if x == 0])
            kl = kv_total_length - kz

            if kl > kz:
                alleles.append(Gene(gene_name, kl, kz))

        self.print_out_alleles(self.filter_contained_alleles(alleles))
        self.identify_alleles_with_100_percent(alleles)

        return alleles

    '''Identify genes which match 100 percent'''

    def identify_alleles_with_100_percent(self, alleles):
        for g in alleles:
            if(g.percentage_coverage() == 100
                    and g.name not in self.genes_with_100_percent):
                self.genes_with_100_percent[g.name] = 1

    '''Some genes are subsets of other genes, so optionally filter out'''

    def filter_contained_alleles(self, alleles):
        if self.no_gene_filter:
            return alleles

        prefix_to_coverage = {}
        filtered_alleles = []
        for g in alleles:
            if g.prefix_short_name() in prefix_to_coverage:
                prefix_to_coverage[g.prefix_short_name()].append(g)
            else:
                prefix_to_coverage[g.prefix_short_name()] = [g]

        for genes in prefix_to_coverage.values():
            genes.sort(key=lambda x: x.percentage_coverage(), reverse=True)
            for index, gene in enumerate(genes):
                if gene.percentage_coverage() == 100 or index == 0:
                    filtered_alleles.append(gene)

        return filtered_alleles

    '''print alleles to stdout'''

    def print_out_alleles(self, alleles):
        found_alleles = False

        for g in alleles:
            if g.percentage_coverage() >= self.min_perc_coverage:
                found_alleles = True
                if self.output_file:
                    with open(self.output_file, 'a+') as output_fh:
                        output_fh.write(str(g) + "\n")
                else:
                    print(g)

        if found_alleles and self.print_interval is not None:
            print("****")

    '''print out the header'''

    def print_out_header(self):
        header = "GENE\tCOMPLETENESS\t%COVERAGE\tACCESSION\tDATABASE\tPRODUCT"
        if self.output_file:
            with open(self.output_file, 'w+') as output_fh:
                output_fh.write(header + "\n")
        else:
            print(header)

    # Derived from https://github.com/sanger-pathogens/Fastaq
    # Author: Martin Hunt
    def open_file_read(self):
        if self.filename == '-':
            f = sys.stdin
        elif self.filename.endswith('.gz'):
            # first check that the file is OK according to gunzip
            retcode = subprocess.call('gunzip -t ' + self.filename, shell=True)
            if retcode > 1:
                raise Error("Error file may not be gzipped correctly")

            # now open the file
            f = os.popen('gunzip -c ' + self.filename)
        else:
            try:
                f = open(self.filename)
            except IOError:
                raise Error("Error opening for reading file '" +
                            self.filename + "'")

        return f
