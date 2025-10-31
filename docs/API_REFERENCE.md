# TipToft API Reference

Complete API documentation for all TipToft modules and classes.

## Core Modules

### tiptoft.TipToft

Main workflow coordinator for plasmid detection.

#### Class: `TipToft`

```python
class TipToft(options)
```

**Description**: Coordinates the overall plasmid typing workflow.

**Parameters**:
- `options`: Parsed command-line arguments object

**Attributes**:
- `logger`: Logging instance
- `plasmid_data` (str): Path to plasmid database FASTA
- `input_fastq` (str): Path to input FASTQ file
- `kmer` (int): K-mer size for matching
- `verbose` (bool): Enable debug logging
- `min_fasta_hits` (int): Minimum k-mer matches threshold
- `print_interval` (int): Print results every N reads
- `output_file` (str): Output file path
- `filtered_reads_file` (str): Path to save matching reads
- `max_gap` (int): Maximum gap for block merging
- `min_block_size` (int): Minimum block size in bases
- `margin` (int): Flanking bases around blocks
- `homopolyer_compression` (bool): Use homopolymer compression
- `max_kmer_count` (int): Max k-mer repetitions to include
- `no_gene_filter` (bool): Disable gene filtering

**Methods**:

##### `run()`

Execute the plasmid detection workflow.

```python
def run() -> None
```

**Returns**: None (results written to output)

**Example**:
```python
from tiptoft.TipToft import TipToft

options = parse_args()  # Your argument parser
tiptoft = TipToft(options)
tiptoft.run()
```

---

### tiptoft.Fasta

Database loading and k-mer indexing.

#### Class: `Fasta`

```python
class Fasta(logger, filename, k, homopolyer_compression, max_kmer_count=10)
```

**Description**: Loads plasmid database and builds k-mer indices.

**Parameters**:
- `logger`: Python logging instance
- `filename` (str): Path to FASTA database
- `k` (int): K-mer size
- `homopolyer_compression` (bool): Apply homopolymer compression
- `max_kmer_count` (int, optional): Max k-mer occurrences (default: 10)

**Attributes**:
- `sequences_to_kmers` (dict): Sequence ID → k-mer counter
- `sequences_to_kmers_count` (dict): Sequence ID → k-mer frequencies
- `all_kmers` (dict): All k-mers → occurrence count
- `kmers_to_genes` (dict): K-mer → list of sequence IDs
- `kmer_keys_set` (set): Set of all k-mers

**Methods**:

##### `sequence_kmers(kmer_action='get_all_kmers_counter')`

Extract k-mers from all sequences.

```python
def sequence_kmers(kmer_action: str = 'get_all_kmers_counter') -> dict
```

**Parameters**:
- `kmer_action` (str): Method name ('get_all_kmers_counter' or 'get_all_kmers_freq')

**Returns**: dict mapping sequence ID to k-mer dictionary

##### `all_kmers_in_file()`

Get all unique k-mers with occurrence counts.

```python
def all_kmers_in_file() -> dict
```

**Returns**: dict mapping k-mer to occurrence count

##### `all_kmers_to_seq_in_file()`

Create k-mer to sequence mapping.

```python
def all_kmers_to_seq_in_file() -> dict
```

**Returns**: dict mapping k-mer to list of sequence IDs

---

### tiptoft.Fastq

FASTQ read processing and matching.

#### Class: `Fastq`

```python
class Fastq(logger, filename, k, all_kmers, min_fasta_hits, print_interval,
            output_file, filtered_reads_file, fasta, homopolyer_compression,
            max_gap=3, min_block_size=130, margin=10, start_time=0,
            min_kmers_for_onex_pass=10, min_perc_coverage=85,
            max_kmer_count=10, no_gene_filter=False)
```

**Description**: Processes reads and matches against database.

**Key Methods**:

##### `read_filter_and_map()`

Main processing loop for all reads.

```python
def read_filter_and_map() -> None
```

**Process**:
1. Parse FASTQ reads
2. Filter reads with quick k-mer pass
3. Detailed analysis of passing reads
4. Output results

---

### tiptoft.Kmers

K-mer extraction and management.

#### Class: `Kmers`

```python
class Kmers(sequence, k, homopolyer_compression)
```

**Description**: Extracts k-mers from DNA sequences.

**Parameters**:
- `sequence` (str): DNA sequence
- `k` (int): K-mer size
- `homopolyer_compression` (bool): Apply homopolymer compression

**Attributes**:
- `sequence` (str): Processed sequence
- `k` (int): K-mer size
- `end` (int): Last valid k-mer start position

**Methods**:

##### `get_all_kmers_counter(max_kmer_count=10)`

Get k-mers with zero counts.

```python
def get_all_kmers_counter(max_kmer_count: int = 10) -> dict
```

**Parameters**:
- `max_kmer_count` (int): Filter k-mers occurring more than this

**Returns**: dict mapping k-mer → 0

**Example**:
```python
from tiptoft.Kmers import Kmers

kmers = Kmers("ATCGATCG", k=3, homopolyer_compression=False)
result = kmers.get_all_kmers_counter()
# {'ATC': 0, 'TCG': 0, 'CGA': 0, 'GAT': 0}
```

##### `get_all_kmers_freq(max_kmer_count=10)`

Get k-mers with frequencies.

```python
def get_all_kmers_freq(max_kmer_count: int = 10) -> dict
```

**Returns**: dict mapping k-mer → frequency count

##### `get_all_kmers_filtered(max_kmer_count=10)`

Get k-mers with positions.

```python
def get_all_kmers_filtered(max_kmer_count: int = 10) -> dict
```

**Returns**: dict mapping k-mer → list of positions

##### `get_one_x_coverage_of_kmers()`

Get non-overlapping k-mers.

```python
def get_one_x_coverage_of_kmers() -> list
```

**Returns**: list of k-mer sequences at k-base intervals

---

### tiptoft.Gene

Gene representation with coverage.

#### Class: `Gene`

```python
class Gene(name, kmers_with_coverage, kmers_without_coverage)
```

**Description**: Represents a detected gene with coverage information.

**Parameters**:
- `name` (str): Full gene name from database
- `kmers_with_coverage` (int): K-mers detected
- `kmers_without_coverage` (int): K-mers not detected

**Methods**:

##### `percentage_coverage()`

Calculate coverage percentage.

```python
def percentage_coverage() -> int
```

**Returns**: int (0-100)

##### `is_full_coverage()`

Check if fully covered.

```python
def is_full_coverage() -> bool
```

**Returns**: True if all k-mers detected

##### `completeness()`

Get completeness status.

```python
def completeness() -> str
```

**Returns**: "Full" or "Partial"

##### `short_name()`

Extract short gene name.

```python
def short_name() -> str
```

**Returns**: Human-readable short name

##### `accession()`

Extract accession number.

```python
def accession() -> str
```

**Returns**: GenBank/EMBL accession

##### `__str__()`

Format for output.

```python
def __str__() -> str
```

**Returns**: Tab-delimited string for output

**Example**:
```python
from tiptoft.Gene import Gene

gene = Gene("rep7.1_repC_AB037671", 100, 0)
print(gene.percentage_coverage())  # 100
print(gene.is_full_coverage())     # True
print(gene.short_name())            # "rep7.1"
print(gene.accession())             # "AB037671"
```

---

### tiptoft.Blocks

Block identification and merging.

#### Class: `Blocks`

```python
class Blocks(k, min_block_size, max_gap, margin)
```

**Description**: Identifies contiguous k-mer match blocks.

**Parameters**:
- `k` (int): K-mer size
- `min_block_size` (int): Minimum block size in bases
- `max_gap` (int): Max gap for merging (k-mer units)
- `margin` (int): Flanking bases to add

**Methods**:

##### `find_all_blocks(sequence_hits)`

Identify all contiguous blocks.

```python
def find_all_blocks(sequence_hits: list) -> list
```

**Parameters**:
- `sequence_hits` (list): Array of k-mer match counts per position

**Returns**: list of [start, end] coordinate pairs

##### `merge_blocks(blocks)`

Merge nearby blocks.

```python
def merge_blocks(blocks: list) -> list
```

**Parameters**:
- `blocks` (list): List of [start, end] pairs

**Returns**: list of merged [start, end] pairs

##### `find_largest_block(sequence_hits)`

Find largest block.

```python
def find_largest_block(sequence_hits: list) -> tuple
```

**Returns**: (start, end) of largest block, or (0, 0) if none meet criteria

##### `adjust_block_start(block_start)`

Adjust start with margin.

```python
def adjust_block_start(block_start: int) -> int
```

**Returns**: Adjusted start position in bases

##### `adjust_block_end(block_end, seq_length)`

Adjust end with margin.

```python
def adjust_block_end(block_end: int, seq_length: int) -> int
```

**Returns**: Adjusted end position in bases

---

### tiptoft.Read

FASTQ read representation.

#### Class: `Read`

```python
class Read(id=None, seq=None, qual=None)
```

**Description**: Represents a FASTQ sequence read.

**Parameters**:
- `id` (str, optional): Read identifier
- `seq` (str, optional): Nucleotide sequence
- `qual` (str, optional): Quality scores

**Methods**:

##### `subsequence(start, end)`

Extract subsequence.

```python
def subsequence(start: int, end: int) -> Read
```

**Returns**: New Read object with subsequence

##### `reverse_complement_sequence()`

Calculate reverse complement.

```python
def reverse_complement_sequence() -> str
```

**Returns**: Reverse complement of sequence

##### `reverse_read()`

Create reverse-complemented read.

```python
def reverse_read() -> Read
```

**Returns**: New Read object with RC sequence

##### `get_next_from_file(f)`

Parse next FASTQ record from file.

```python
def get_next_from_file(f: file) -> Read | bool
```

**Returns**: Self if successful, False if EOF

##### `__str__()`

Format as FASTQ.

```python
def __str__() -> str
```

**Returns**: Four-line FASTQ format string

---

### tiptoft.InputTypes

Input validation utilities.

#### Class: `InputTypes`

Static methods for command-line validation.

##### `is_fastq_file_valid(filename)`

Validate FASTQ file.

```python
@staticmethod
def is_fastq_file_valid(filename: str) -> str
```

**Raises**: argparse.ArgumentTypeError if invalid

##### `is_plasmid_file_valid(filename)`

Validate plasmid database file.

```python
@staticmethod
def is_plasmid_file_valid(filename: str) -> str
```

**Raises**: argparse.ArgumentTypeError if invalid

##### `is_kmer_valid(value_str)`

Validate k-mer size.

```python
@staticmethod
def is_kmer_valid(value_str: str) -> int
```

**Returns**: int k-mer size (7-31)

**Raises**: argparse.ArgumentTypeError if invalid

---

### tiptoft.RefGenesGetter

Database downloader.

#### Class: `RefGenesGetter`

```python
class RefGenesGetter(verbose=False)
```

**Description**: Downloads PlasmidFinder database.

**Parameters**:
- `verbose` (bool): Keep temporary files

**Methods**:

##### `run(outprefix)`

Download and process database.

```python
def run(outprefix: str) -> None
```

**Parameters**:
- `outprefix` (str): Prefix for output files

**Creates**:
- `{outprefix}.fa`: FASTA sequences
- `{outprefix}.tsv`: Metadata

---

### tiptoft.TipToftDatabaseDownloader

Database downloader driver.

#### Class: `TipToftDatabaseDownloader`

```python
class TipToftDatabaseDownloader(options)
```

**Description**: Coordinates database download.

**Parameters**:
- `options`: Parsed command-line arguments

**Methods**:

##### `run()`

Execute download.

```python
def run() -> None
```

---

## Helper Functions

### homopolymer_compression

```python
def homopolymer_compression_of_sequence(sequence: str) -> str
```

Compress homopolymer runs to single bases.

**Parameters**:
- `sequence` (str): DNA sequence

**Returns**: Compressed sequence

**Example**:
```python
from homopolymer_compression import homopolymer_compression_of_sequence

compressed = homopolymer_compression_of_sequence("AAATTTGGG")
print(compressed)  # "ATG"
```

---

## Data Structures

### K-mer Dictionary Format

```python
{
    'ATCG': [0, 5, 10],  # K-mer with positions
    'TCGA': [1, 6],
    'CGAT': [2, 7]
}
```

### Gene Results Format

```python
{
    'gene_id': {
        'kmers_with_coverage': 95,
        'kmers_without_coverage': 5,
        'percentage': 95
    }
}
```

### Output Format

Tab-delimited text:
```
GENE    COMPLETENESS    %COVERAGE    ACCESSION    DATABASE    PRODUCT
rep7.1  Full            100          AB037671     plasmidfinder    rep7.1_repC...
```

---

## Usage Examples

### Basic Analysis

```python
#!/usr/bin/env python3
from tiptoft.TipToft import TipToft
from argparse import Namespace

# Create options
options = Namespace(
    input_fastq='reads.fastq.gz',
    plasmid_data=None,  # Use bundled
    kmer=13,
    verbose=False,
    min_fasta_hits=8,
    print_interval=None,
    output_file='results.txt',
    filtered_reads_file=None,
    max_gap=3,
    min_block_size=50,
    margin=10,
    min_kmers_for_onex_pass=5,
    min_perc_coverage=85,
    no_hc_compression=False,
    max_kmer_count=10,
    no_gene_filter=False
)

# Run analysis
tiptoft = TipToft(options)
tiptoft.run()
```

### Custom K-mer Extraction

```python
from tiptoft.Kmers import Kmers

# Extract k-mers from sequence
sequence = "ATCGATCGATCG"
kmers = Kmers(sequence, k=4, homopolyer_compression=False)

# Get all k-mers with positions
all_kmers = kmers.get_all_kmers_filtered(max_kmer_count=10)
print(all_kmers)
# {'ATCG': [0, 4, 8], 'TCGA': [1, 5], 'CGAT': [2, 6], 'GATC': [3, 7]}

# Get 1x coverage sample
onex = kmers.get_one_x_coverage_of_kmers()
print(onex)
# ['ATCG', 'ATCG', 'ATCG']
```

### Gene Analysis

```python
from tiptoft.Gene import Gene

# Create gene with coverage data
gene = Gene(
    name="rep7.1_repC_AB037671",
    kmers_with_coverage=95,
    kmers_without_coverage=5
)

# Check coverage
print(f"Coverage: {gene.percentage_coverage()}%")
print(f"Complete: {gene.is_full_coverage()}")
print(f"Status: {gene.completeness()}")

# Output formatting
print(gene)  # Tab-delimited output line
```

---

## Error Handling

### Common Exceptions

- `FileNotFoundError`: Input file not found
- `argparse.ArgumentTypeError`: Invalid command-line argument
- `IOError`: File read/write error
- `MemoryError`: Insufficient memory
- `KeyboardInterrupt`: User cancelled

### Example Error Handling

```python
try:
    tiptoft = TipToft(options)
    tiptoft.run()
except FileNotFoundError as e:
    print(f"Error: Input file not found - {e}")
except MemoryError:
    print("Error: Insufficient memory. Try smaller input or more RAM.")
except KeyboardInterrupt:
    print("\nAnalysis cancelled by user")
```

---

## Performance Considerations

### Memory Usage

- **Database loading**: ~50-100 MB depending on database size
- **Per read processing**: ~1-10 KB per read
- **Total**: Linear with input size

### Time Complexity

- **Database indexing**: O(D × K) where D = database size, K = k-mer size
- **Read processing**: O(R × L × K) where R = read count, L = read length
- **Overall**: Linear with input size

### Optimization Tips

1. Use larger k-mer for faster analysis (if error rate permits)
2. Increase min_fasta_hits to filter more reads early
3. Process multiple samples in parallel
4. Use SSD storage for large FASTQ files

---

## Version Information

This API reference is for TipToft version 1.0+. Check your version:

```bash
tiptoft --version
```

For the latest documentation, visit:
https://github.com/andrewjpage/tiptoft
