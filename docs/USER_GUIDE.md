# TipToft User Guide

## Overview

TipToft is a bioinformatics tool that predicts which plasmids are present in uncorrected long-read sequencing data from PacBio or Oxford Nanopore technologies. It uses k-mer matching to identify plasmid incompatibility groups and replicon types even when plasmids are missing from genome assemblies.

## Why Use TipToft?

Long read assemblies can miss plasmids when they are:
- Very small or very large
- Present at high or low copy numbers relative to the chromosome
- Circular and difficult to assemble

TipToft analyzes raw reads to predict expected plasmids, helping identify assembly issues.

## Quick Start

### Installation

```bash
# Install dependencies
pip3 install cython biopython pyfastaq

# Install TipToft
pip3 install tiptoft

# Or install latest development version
pip3 install git+git://github.com/andrewjpage/tiptoft.git
```

### Basic Usage

```bash
# Analyze reads with default settings
tiptoft reads.fastq.gz

# Use custom k-mer size (for PacBio)
tiptoft -k 15 reads.fastq.gz

# Save results to file
tiptoft -o results.txt reads.fastq.gz

# Save matching reads for further analysis
tiptoft -f plasmid_reads.fastq reads.fastq.gz
```

## Detailed Usage

### Input Requirements

**FASTQ File**: Your long reads in FASTQ format (optionally gzipped)
- Quality scores are ignored (only sequence is used)
- Can read from stdin using `-` as filename
- Supports both PacBio and Oxford Nanopore reads

**Plasmid Database**: FASTA file with plasmid typing sequences
- Default: Uses bundled PlasmidFinder database
- Update: Use `tiptoft_database_downloader` for latest version

### K-mer Selection

The k-mer size (`-k`) is the most critical parameter:

| Read Type | Recommended k-mer | Rationale |
|-----------|-------------------|-----------|
| Nanopore | 13 (default) | ~1 error per 10 bases |
| PacBio | 15 | Lower error rate |
| High quality | 17-31 | Longer accurate stretches |

**Rule of thumb**: Set k to the average error-free stretch length in your reads.

### Output Format

Tab-delimited output with columns:

1. **GENE**: Short name of the replicon/incompatibility group
2. **COMPLETENESS**: Full (all k-mers found) or Partial (some missing)
3. **%COVERAGE**: Percentage of k-mers detected (0-100)
4. **ACCESSION**: GenBank/EMBL accession number
5. **DATABASE**: Source database (plasmidfinder)
6. **PRODUCT**: Full description from database

Example:
```
GENE         COMPLETENESS  %COVERAGE  ACCESSION    DATABASE       PRODUCT
rep7.1       Full          100        AB037671     plasmidfinder  rep7.1_repC(Cassette)_AB037671
rep20.11     Full          100        AP003367     plasmidfinder  rep20.11_repA(VRSAp)_AP003367
```

### Common Workflows

#### 1. Basic Plasmid Detection

```bash
tiptoft sample.fastq.gz > plasmids.txt
```

#### 2. Extract Plasmid Reads for Assembly

```bash
# Identify and save plasmid-containing reads
tiptoft -f plasmid_reads.fastq sample.fastq.gz

# Assemble just the plasmid reads
# (use your preferred assembler)
flye --nano-raw plasmid_reads.fastq --out-dir plasmid_assembly
```

#### 3. Real-time Streaming Analysis

For long sequencing runs, get intermediate results:

```bash
# Print results every 1000 reads
tiptoft -p 1000 sample.fastq.gz

# Stop when you have enough information
# Results separated by ****
```

#### 4. Comprehensive Analysis with Latest Database

```bash
# Download latest PlasmidFinder database
tiptoft_database_downloader my_db

# Use custom database
tiptoft -d my_db.fa sample.fastq.gz
```

## Advanced Parameters

### Sensitivity Tuning

**`--min_fasta_hits` (-m)**: Minimum k-mer matches (default: 8)
- Lower = more sensitive, slower
- Higher = less sensitive, faster
- Adjust based on read length and coverage

**`--min_perc_coverage` (-c)**: Minimum coverage to report (default: 85%)
- Lower = report partial matches
- Higher = only report strong hits
- Useful for filtering low-confidence hits

**`--min_kmers_for_onex_pass`**: Initial screening threshold (default: 5)
- First pass uses 1x k-mer coverage for speed
- Reads passing this go to detailed analysis
- Lower values = more thorough but slower

### Algorithm Parameters

**`--max_gap`**: Gap tolerance for block merging (default: 3)
- Measured in k-mer multiples
- Allows bridging error-rich regions
- Higher = more tolerant of errors

**`--min_block_size`**: Minimum block size (default: 50 bases)
- Filters small spurious matches
- Increase if getting too many false positives

**`--margin`**: Flanking region around blocks (default: 10 bases)
- Added context for mapping
- Helps catch k-mers near block boundaries

**`--max_kmer_count`**: Max k-mer repetition (default: 10)
- Filters highly repetitive k-mers
- Prevents false matches to repetitive regions

### Special Options

**`--no_hc_compression`**: Disable homopolymer compression
- Default: Enabled (recommended for long reads)
- Helps handle homopolymer length errors
- Disable for high-quality reads

**`--no_gene_filter`**: Disable gene filtering
- Default: Filters lower coverage genes from same group
- Disable to see all matches

**`--verbose` (-v)**: Debug mode
- Enables detailed logging
- Shows performance profiling
- Useful for troubleshooting

## Interpreting Results

### Coverage Interpretation

- **100% Full**: Entire typing sequence detected → High confidence
- **95-99% Full**: Nearly complete → Likely present
- **85-94% Partial**: Partial match → May be variant or related
- **<85%**: Not reported by default → Consider false positive

### Multiple Hits

Finding multiple plasmids is normal:
- **Same group, different alleles**: Common variants
- **Different groups**: Multiple plasmids present
- **Related groups**: May be misclassification

### No Hits Found

Possible reasons:
1. No plasmids in sample
2. K-mer size too large for error rate
3. Novel plasmid not in database
4. Very low plasmid coverage

Try:
- Reduce k-mer size
- Update database
- Check read quality

## Troubleshooting

### Error: "Cannot access input file"

- Check file path is correct
- Verify file permissions
- For stdin, use `-` as filename

### No Results or Very Few Results

**Try**:
- Reduce k-mer size: `-k 11` or `-k 9`
- Lower coverage threshold: `-c 70`
- Check if database is correct: `-d` option

### Too Many False Positives

**Try**:
- Increase k-mer size: `-k 15` or `-k 17`
- Raise coverage threshold: `-c 95`
- Increase min hits: `-m 15`

### Performance Issues

**Speed up**:
- Increase `--min_kmers_for_onex_pass`
- Increase `--min_fasta_hits`
- Increase `--min_block_size`

**More thorough** (slower):
- Decrease thresholds above
- Use smaller k-mer size

## Citations

When using TipToft, please cite:

**TipToft**:
```
AJ Page, T Seemann (2019). TipToft: detecting plasmids contained in 
uncorrected long read sequencing data. Journal of Open Source Software, 
4(35), 1021, https://doi.org/10.21105/joss.01021
```

**PlasmidFinder** (database):
```
Carattoli et al, In Silico Detection and Typing of Plasmids using 
PlasmidFinder and Plasmid Multilocus Sequence Typing, 
Antimicrob Agents Chemother. 2014;58(7):3895–3903.
```

## Getting Help

- **Issues**: https://github.com/andrewjpage/tiptoft/issues
- **Documentation**: https://github.com/andrewjpage/tiptoft
- **Email**: andrew.page@quadram.ac.uk

## Example Data

Test TipToft with example data:

```bash
# Download example data
wget https://github.com/andrewjpage/tiptoft/raw/master/example_data/ERS654932_plasmids.fastq.gz

# Run TipToft
tiptoft ERS654932_plasmids.fastq.gz

# Compare to expected output
wget https://raw.githubusercontent.com/andrewjpage/tiptoft/master/example_data/expected_output
```
