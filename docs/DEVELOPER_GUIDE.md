# TipToft Developer Guide

## Architecture Overview

TipToft is a Python-based tool for plasmid detection in long-read sequencing data using k-mer matching. The architecture follows a modular design with clear separation of concerns.

### Core Components

```
┌─────────────────────────────────────────────────────────────┐
│                     Command Line Scripts                     │
│  tiptoft                    tiptoft_database_downloader     │
└──────────────────┬───────────────────────┬──────────────────┘
                   │                       │
                   v                       v
         ┌─────────────────┐    ┌──────────────────────┐
         │    TipToft      │    │ TipToftDatabase      │
         │    (Driver)     │    │ Downloader           │
         └────────┬────────┘    └──────────┬───────────┘
                  │                        │
                  v                        v
         ┌────────────────┐      ┌──────────────────┐
         │  Fasta/Fastq   │      │  RefGenesGetter  │
         │  (Processors)  │      │  (Downloader)    │
         └────────┬───────┘      └──────────────────┘
                  │
         ┌────────┴──────────┬──────────┬─────────┐
         v                   v          v         v
    ┌────────┐         ┌────────┐  ┌──────┐  ┌──────┐
    │ Kmers  │         │ Blocks │  │ Gene │  │ Read │
    └────────┘         └────────┘  └──────┘  └──────┘
```

## Module Documentation

### TipToft.py
**Purpose**: Main workflow coordinator

**Responsibilities**:
- Parse command-line options
- Initialize Fasta database loader
- Initialize Fastq processor
- Coordinate execution

**Key Methods**:
- `__init__(options)`: Configure from command-line args
- `run()`: Execute the plasmid detection workflow

### Fasta.py
**Purpose**: Load and index plasmid database

**Responsibilities**:
- Parse FASTA file of plasmid sequences
- Extract all k-mers from each sequence
- Build k-mer → gene lookup tables
- Provide k-mer access interface

**Key Methods**:
- `sequence_kmers()`: Extract k-mers from all sequences
- `all_kmers_in_file()`: Build master k-mer dictionary
- `all_kmers_to_seq_in_file()`: Create k-mer → sequence mapping

**Data Structures**:
```python
sequences_to_kmers = {
    'gene1': {'ATCG': 0, 'TCGA': 0, ...},
    'gene2': {'CGAT': 0, 'GATC': 0, ...}
}

kmers_to_genes = {
    'ATCG': ['gene1', 'gene3', ...],
    'TCGA': ['gene1', 'gene2', ...]
}
```

### Fastq.py
**Purpose**: Process reads and match against database

**Responsibilities**:
- Parse FASTQ input (file or stdin)
- Extract k-mers from each read
- Match read k-mers against database
- Identify and report plasmid sequences
- Optionally save matching reads

**Algorithm**:
1. **Initial Screening**: Fast 1x k-mer pass to filter reads
2. **Block Finding**: Identify contiguous matching regions
3. **Detailed Matching**: Align largest blocks against database
4. **Coverage Calculation**: Compute k-mer coverage for each gene
5. **Filtering**: Apply coverage threshold and gene filters
6. **Reporting**: Output results in tab-delimited format

**Key Methods**:
- `read_filter_and_map()`: Main processing loop
- `does_read_contain_quick_pass_kmers()`: Fast pre-filter
- `find_largest_region_and_analyse()`: Detailed analysis
- `output_gene_results()`: Format and print results

### Kmers.py
**Purpose**: K-mer extraction and management

**Responsibilities**:
- Extract k-mers from DNA sequences
- Apply optional homopolymer compression
- Filter repetitive k-mers
- Provide various k-mer access methods

**Key Features**:
- Sliding window k-mer extraction
- Homopolymer compression for error tolerance
- Position tracking for each k-mer
- Frequency counting

**Key Methods**:
- `get_all_kmers_counter()`: K-mers with zero counts
- `get_all_kmers_freq()`: K-mers with frequencies
- `get_all_kmers_filtered()`: K-mers with positions, filtered
- `get_one_x_coverage_of_kmers()`: Non-overlapping k-mers

### Blocks.py
**Purpose**: Identify and merge k-mer match regions

**Responsibilities**:
- Find contiguous blocks of k-mer matches
- Merge nearby blocks (error tolerance)
- Adjust block boundaries with margins
- Select largest block for analysis

**Algorithm**:
```python
# Sequence:  AAAAA-----TTTTTT--GGGGG
# K-mers:    11111-----111111--11111
# Blocks:    [----]     [-----]  [---]
# Merged:    [------------------]  [---]
# Largest:   [------------------]
```

**Key Methods**:
- `find_all_blocks()`: Identify contiguous regions
- `merge_blocks()`: Combine nearby blocks
- `find_largest_block()`: Select best candidate
- `adjust_block_start/end()`: Add margins and bounds check

### Gene.py
**Purpose**: Represent detected genes with coverage

**Responsibilities**:
- Store gene identification and coverage data
- Calculate coverage percentage
- Determine completeness (Full/Partial)
- Parse gene names for output formatting

**Key Methods**:
- `percentage_coverage()`: Calculate % k-mers detected
- `is_full_coverage()`: Check if all k-mers present
- `short_name()`: Extract human-readable name
- `accession()`: Extract GenBank accession

### Read.py
**Purpose**: Represent FASTQ reads

**Responsibilities**:
- Store read ID, sequence, and quality
- Extract subsequences
- Calculate reverse complement
- Parse FASTQ format

**Key Methods**:
- `subsequence()`: Extract region of read
- `reverse_complement_sequence()`: DNA RC
- `reverse_read()`: Create RC read object
- `get_next_from_file()`: Parse FASTQ record

### InputTypes.py
**Purpose**: Command-line argument validation

**Responsibilities**:
- Validate file paths exist
- Validate k-mer size in valid range
- Provide clear error messages

### RefGenesGetter.py
**Purpose**: Download PlasmidFinder database

**Responsibilities**:
- Download FASTA files from BitBucket
- Combine multiple files
- Normalize sequence IDs
- Create metadata TSV

### TipToftDatabaseDownloader.py
**Purpose**: Database downloader driver

**Responsibilities**:
- Parse command-line options
- Initialize RefGenesGetter
- Execute download workflow

## Development Setup

### Prerequisites

```bash
# System dependencies
sudo apt-get install python3 python3-pip python3-dev gcc

# Python dependencies
pip3 install cython biopython pyfastaq nose pytest pytest-cov
```

### Building from Source

```bash
# Clone repository
git clone https://github.com/andrewjpage/tiptoft.git
cd tiptoft

# Install in development mode
pip3 install -e .
```

### Running Tests

```bash
# Run all tests
python3 -m pytest tiptoft/tests/

# Run with coverage
python3 -m pytest tiptoft/tests/ --cov=tiptoft --cov-report=html

# Run specific test file
python3 -m pytest tiptoft/tests/Kmers_test.py -v

# Run specific test
python3 -m pytest tiptoft/tests/Kmers_test.py::TestKmers::test_four_kmers -v
```

### Code Style

- Follow PEP 8 guidelines
- Use docstrings for classes and methods
- Comment complex algorithms
- Keep functions focused and small

### Testing Guidelines

- Write tests for all new functionality
- Maintain >80% code coverage
- Use descriptive test names
- Include edge cases
- Mock external dependencies

Example test structure:
```python
import unittest
from tiptoft.Module import Class

class TestClass(unittest.TestCase):
    def setUp(self):
        # Setup test fixtures
        pass
        
    def test_specific_behavior(self):
        # Test one specific behavior
        obj = Class(...)
        result = obj.method()
        self.assertEqual(result, expected)
        
    def tearDown(self):
        # Clean up
        pass
```

## Algorithm Details

### K-mer Matching Strategy

1. **Database Indexing** (startup):
   - Extract all k-mers from plasmid sequences
   - Build k-mer → gene mapping
   - Filter overly repetitive k-mers

2. **Read Processing** (per read):
   - Extract read k-mers with homopolymer compression
   - Quick 1x pass: Check if minimum k-mers present
   - If pass: Extract all k-mers (overlapping)
   - Match against database k-mers

3. **Block Finding**:
   - Create array: position → k-mer match count
   - Find contiguous regions of matches
   - Merge blocks within max_gap distance
   - Select largest block for analysis

4. **Gene Coverage**:
   - Extract largest block from read
   - Match all k-mers in block against database
   - For each gene: count matched vs. total k-mers
   - Calculate coverage percentage
   - Mark as Full if 100%, else Partial

5. **Filtering and Output**:
   - Filter genes below min_perc_coverage
   - Optionally filter duplicate groups
   - Sort and output results

### Homopolymer Compression

Long reads often have errors in homopolymer lengths:
- `AAAA` might be sequenced as `AAA` or `AAAAA`

Solution: Compress homopolymers to single base:
- `AAAA` → `A`
- `TTTTGGGCC` → `TGC`

Apply to both database and reads consistently.

### Error Tolerance

Multiple mechanisms handle sequencing errors:

1. **K-mer size**: Smaller k-mers less affected by errors
2. **Homopolymer compression**: Handles length errors
3. **Block merging**: Bridges error-rich regions
4. **Coverage calculation**: Partial matches still reported

## Performance Optimization

### Speed Improvements

1. **Two-pass strategy**: Fast initial filter, detailed analysis only on candidates
2. **K-mer sets**: O(1) lookup for k-mer matching
3. **Early termination**: Skip reads with too few k-mer matches
4. **Filtering repetitive k-mers**: Reduces false positives

### Memory Optimization

1. **Streaming**: Process reads one at a time
2. **Efficient data structures**: Dictionaries for k-mer storage
3. **Optional C extension**: Homopolymer compression in C

### Typical Performance

- **Memory**: ~80 MB for bundled database
- **Speed**: ~1 minute for 800 MB FASTQ
- **Scaling**: Linear with input size

## Contributing

### Workflow

1. Fork the repository
2. Create a feature branch
3. Make changes with tests
4. Ensure tests pass and coverage maintained
5. Submit pull request with description

### Pull Request Checklist

- [ ] Tests added for new functionality
- [ ] All tests pass
- [ ] Code coverage ≥80%
- [ ] Code follows style guidelines
- [ ] Documentation updated
- [ ] Commit messages are clear

### Release Process

1. Update VERSION file
2. Update CHANGELOG
3. Tag release: `git tag -a v1.2.3 -m "Release 1.2.3"`
4. Push tag: `git push origin v1.2.3`
5. GitHub Actions builds and publishes to PyPI

## Debugging

### Enable Verbose Mode

```bash
tiptoft -v sample.fastq.gz
```

Shows:
- Detailed logging
- Performance profiling
- K-mer statistics

### Common Issues

**Import errors**: Check dependencies installed
```bash
pip3 install biopython pyfastaq cython
```

**C extension errors**: Use Python fallback
```bash
# Already handled automatically in current version
```

**Memory errors**: Process in chunks or use smaller database

### Profiling

Built-in profiling with `-v` flag shows:
- Function call counts
- Cumulative time
- Per-call time

For deeper profiling:
```bash
python3 -m cProfile -s cumulative scripts/tiptoft sample.fastq.gz
```

## Future Development

Potential enhancements:
- Multi-threading for parallel processing
- GPU acceleration for k-mer matching
- Additional database support
- Variant calling for known plasmids
- Assembly integration
- Quality score incorporation

## Additional Resources

- **Repository**: https://github.com/andrewjpage/tiptoft
- **Issues**: https://github.com/andrewjpage/tiptoft/issues
- **PyPI**: https://pypi.org/project/tiptoft/
- **Docker**: https://hub.docker.com/r/andrewjpage/tiptoft
