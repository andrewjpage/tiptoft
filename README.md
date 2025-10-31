# TipToft

> **Plasmid incompatibility group prediction from uncorrected long reads**

TipToft identifies plasmids in raw, uncorrected long-read sequencing data from PacBio or Oxford Nanopore platforms. Using k-mer matching against the PlasmidFinder database, it predicts which plasmids are present even when they are missing from genome assemblies.

[![Build Status](https://travis-ci.org/andrewjpage/tiptoft.svg?branch=master)](https://travis-ci.org/andrewjpage/tiptoft)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/andrewjpage/tiptoft/blob/master/LICENSE)
[![codecov](https://codecov.io/gh/andrewjpage/tiptoft/branch/master/graph/badge.svg)](https://codecov.io/gh/andrewjpage/tiptoft)
[![Docker Build Status](https://img.shields.io/docker/build/andrewjpage/tiptoft.svg)](https://hub.docker.com/r/andrewjpage/tiptoft)
[![Docker Pulls](https://img.shields.io/docker/pulls/andrewjpage/tiptoft.svg)](https://hub.docker.com/r/andrewjpage/tiptoft)

## Why TipToft?

Long-read assemblies frequently miss plasmids due to:
- **Size extremes**: Very small or very large plasmids
- **Copy number issues**: Too high or too low relative to chromosome
- **Circularity**: Difficult to resolve circular structures

TipToft analyzes **raw reads** to predict expected plasmids, helping you:
- Identify missing plasmids in assemblies
- Validate assembly completeness
- Detect antimicrobial resistance plasmids
- Guide targeted plasmid assembly

## Key Features

- ✅ **Direct read analysis** - No assembly required
- ✅ **Error tolerant** - Handles high error rates in long reads
- ✅ **Fast** - Processes 800 MB in under 1 minute
- ✅ **Memory efficient** - ~80 MB RAM usage
- ✅ **Comprehensive** - Detects incompatibility groups and replicons
- ✅ **Well tested** - 92% code coverage with comprehensive test suite  

## Quick Start

```bash
# Install TipToft
pip3 install cython
pip3 install tiptoft

# Analyze your reads
tiptoft my_reads.fastq.gz

# Save results and matching reads
tiptoft -o results.txt -f plasmid_reads.fastq my_reads.fastq.gz
```

## Table of Contents

- [Paper & Citation](#paper--citation)
- [Installation](#installation)
  - [PyPI](#pypi)
  - [Debian/Ubuntu](#debianubuntu-trustyxenial)
  - [Docker](#docker)
  - [Homebrew](#homebrew)
  - [Bioconda](#bioconda)
- [Usage](#usage)
- [Documentation](#documentation)
- [Output](#output)
- [Example](#example-usage)
- [Resource Usage](#resource-usage)
- [Contributing](#contribute-to-the-software)
- [License](#license)
- [Contact](#feedbackissues)

## Paper & Citation
[![DOI](http://joss.theoj.org/papers/10.21105/joss.01021/status.svg)](https://doi.org/10.21105/joss.01021)

AJ Page, T Seemann (2019). TipToft: detecting plasmids contained in uncorrected long read sequencing data. Journal of Open Source Software, 4(35), 1021, https://doi.org/10.21105/joss.01021

Please remember to cite the plasmidFinder paper as their database makes this software work:

Carattoli *et al*, *In Silico Detection and Typing of Plasmids using PlasmidFinder and Plasmid Multilocus Sequence Typing*, **Antimicrob Agents Chemother.** 2014;58(7):3895–3903. [view](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4068535/)


## Installation

### PyPI

The easiest installation method using pip:
```bash
# Install dependencies
pip3 install cython

# Install TipToft
pip3 install tiptoft
```

**Development version:**
```bash
pip3 install git+git://github.com/andrewjpage/tiptoft.git
```

**Requirements:**
- Python 3.6+
- Cython
- BioPython ≥ 1.68
- pyfastaq ≥ 3.12.0
- C compiler (gcc, clang, etc.)

### Debian/Ubuntu (Trusty/Xenial)
```bash
sudo apt-get update -qq
sudo apt-get install -y git python3 python3-setuptools python3-biopython python3-pip
pip3 install cython
pip3 install tiptoft
```

### Docker
```bash
# Pull Docker image
docker pull andrewjpage/tiptoft

# Run on your data
docker run --rm -it -v /path/to/data:/data andrewjpage/tiptoft tiptoft /data/reads.fastq.gz
```

Example with bundled test data:
```bash
docker run --rm -it -v /path/to/example_data:/example_data andrewjpage/tiptoft tiptoft /example_data/ERS654932_plasmids.fastq.gz
```

### Homebrew
For macOS users:

```bash
# Install Homebrew (if needed)
# See https://brew.sh/

# Install Python 3
brew install python

# Install TipToft
pip3 install cython
pip3 install tiptoft
```

### Bioconda
```bash
# Install Bioconda (if needed)
# See http://bioconda.github.io/

# Install TipToft
conda install tiptoft
```

### Windows (WSL)

Windows users can use Ubuntu on Windows 10 (WSL):
```bash
# In Ubuntu on Windows terminal
sudo apt-get update
sudo apt-get install python3 python3-pip
pip3 install cython
pip3 install tiptoft
```

**Note:** Native Windows is not officially supported. Use WSL or Docker for best results.

## Documentation

Comprehensive documentation is available in the `docs/` directory:

- **[User Guide](docs/USER_GUIDE.md)** - Installation, usage, parameters, troubleshooting
- **[Developer Guide](docs/DEVELOPER_GUIDE.md)** - Architecture, API, contributing
- **[Work Instructions](docs/WORK_INSTRUCTION.md)** - SOPs for production use
- **[API Reference](docs/API_REFERENCE.md)** - Complete API documentation

## Usage
## tiptoft_database_downloader script
First of all you need plasmid database from PlasmidFinder. There is a snapshot bundled with this repository for your convenience, or alternatively you can use the downloader script to get the latest data. You will need internet access for this step. Please remember to cite the PlasmidFinder paper.

```
usage: tiptoft_database_downloader [options] output_prefix

Download PlasmidFinder database

positional arguments:
  output_prefix  Output prefix

optional arguments:
  -h, --help     show this help message and exit
  --verbose, -v  Turn on debugging (default: False)
  --version      show program's version number and exit
```

Just run:
```
tiptoft_database_downloader 
```
You will now have a file called 'plasmid_files.fa' which can be used with the main script.

## tiptoft script
This is the main script of the application. The mandatory inputs are a FASTQ file of long reads, which can be optionally gzipped.
```
usage: tiptoft [options] input.fastq

plasmid incompatibility group prediction from uncorrected long reads

positional arguments:
  input_fastq           Input FASTQ file (optionally gzipped)

optional arguments:
  -h, --help            show this help message and exit

Optional input arguments:
  --plasmid_data PLASMID_DATA, -d PLASMID_DATA
                        FASTA file containing plasmid data from downloader
                        script, defaults to bundled database (default: None)
  --kmer KMER, -k KMER  k-mer size (default: 13)

Optional output arguments:
  --filtered_reads_file FILTERED_READS_FILE, -f FILTERED_READS_FILE
                        Filename to save matching reads to (default: None)
  --output_file OUTPUT_FILE, -o OUTPUT_FILE
                        Output file [STDOUT] (default: None)
  --print_interval PRINT_INTERVAL, -p PRINT_INTERVAL
                        Print results every this number of reads (default:
                        None)
  --verbose, -v         Turn on debugging [False]
  --version             show program's version number and exit

Optional advanced input arguments:
  --max_gap MAX_GAP     Maximum gap for blocks to be contigous, measured in
                        multiples of the k-mer size (default: 3)
  --margin MARGIN       Flanking region around a block to use for mapping
                        (default: 10)
  --min_block_size MIN_BLOCK_SIZE
                        Minimum block size in bases (default: 130)
  --min_fasta_hits MIN_FASTA_HITS, -m MIN_FASTA_HITS
                        Minimum No. of kmers matching a read (default: 10)
  --min_perc_coverage MIN_PERC_COVERAGE, -c MIN_PERC_COVERAGE
                        Minimum percentage coverage of typing sequence to
                        report (default: 85)
  --min_kmers_for_onex_pass MIN_KMERS_FOR_ONEX_PASS
                        Minimum No. of kmers matching a read in 1st pass
                        (default: 10)
```

### Required argument

__input_fastq__: This is a single FASTQ file. It can be optionally gzipped. Alternatively input can be read from stdin by using the dash character (-) as the input file name. The file must contain long reads, such as those from PacBio or Oxford Nanopore. The quality scores are ignored.

### Optional input arguments

__plasmid_data__: This is a FASTA file containing all of the plasmid typing sequences. This is generated by the tiptoft_database_downloader script. It comes from the PlasmidFinder website, so please be sure to cite their paper (citation gets printed every time you run the script).

__kmer__:  The most important parameter. 13 works well for Nanopore, 15 works well for PacBio, but you may need to play around with it for your data. Long reads have a high error rate, so if you set this too high, nothing will match (because it will contain errors). If you set it too low, everything will match, which isnt much use to you. Thinking about your data, on average how long of a stretch of bases can you get in your read without errors? This is what you should set your kmer to. For example, if you have an average of 1 error every 10 bases, then the ideal kmer would be 9.

### Optional output arguments

__filtered_reads_file__: Save the reads which contain the rep/inc sequences to a new FASTQ file. This is useful if you want to undertake a further assembly just on the plasmids.This file should not already exist. 

__output_file OUTPUT_FILE__: By default the results are printed to STDOUT. If you provide an output filename (which must not exist already), it will print the results to the file.

__print_interval__: By default the whole file is processed and the final results are printed out. However you can get intermediate results printed after every X number of reads, which is useful if you are doing real time streaming of data into the application and can halt when you have enough information. They are separated by "****". 

__verbose__: Enable debugging mode where lots of extra output is printed to STDOUT.

__version__: Print the version number and exit.


### Optional advanced input arguments

__max_gap__: Maximum gap for blocks to be contigous, measured in multiples of the k-mer size. This allows for short regions of elevated errors in the reads to be spanned.

__margin__:  Expand the analysis to look at a few bases on either side of where the sequence is predicted to be on the read. This allows for k-mers to overlap the ends.

__min_block_size__:  This is the minimum sub read size of a read to consider for indepth analysis after matching k-mers have been identified in the read. This speeds up the analysis quite a bit, but there is the risk that some reads may be missed, particularly if they have partial rep/inc sequences.

__min_fasta_hits__: This is the minimum number of matching kmers in a read, for the read to be considered for analysis. It is a hard minimum threshold to speed up analysis.

__min_perc_coverage__: Only report rep/inc sequences above this percentage coverage. Coverage in this instance is kmer coverage of the underlying sequence (rather than depth of coverage).

__min_kmers_for_onex_pass__: The number of k-mers that must be present in the read for the initial onex pass of the database to be considered for further analysis. This speeds up the analysis quite a bit, but there is the risk that some reads may be missed, particularly if they have partial rep/inc sequences.

# Output
The output is tab delmited and printed to STDOUT by default. You can optionally print it to a file using the '-o' parameter. If you would like to see intermediate results, you can tell it to print every X reads with the '-p' parameter, separated by '****'.   An example of the output is:

```
GENE	COMPLETENESS	%COVERAGE	ACCESSION	DATABASE	PRODUCT
rep7.1	Full	100	AB037671	plasmidfinder	rep7.1_repC(Cassette)_AB037671
rep7.5	Partial	99	AF378372	plasmidfinder	rep7.5_CDS1(pKC5b)_AF378372
rep7.6	Partial	94	SAU38656	plasmidfinder	rep7.6_ORF(pKH1)_SAU38656
rep7.9	Full	100	NC007791	plasmidfinder	rep7.9_CDS3(pUSA02)_NC007791
rep7.10	Partial	91	NC_010284.1	plasmidfinder	rep7.10_repC(pKH17)_NC_010284.1
rep7.12	Partial	93	GQ900417.1	plasmidfinder	rep7.12_rep(SAP060B)_GQ900417.1
rep7.17	Full	100	AM990993.1	plasmidfinder	rep7.17_repC(pS0385-1)_AM990993.1
rep20.11	Full	100	AP003367	plasmidfinder	rep20.11_repA(VRSAp)_AP003367
repUS14.	Full	100	AP003367	plasmidfinder	repUS14._repA(VRSAp)_AP003367
```

__GENE__: The first column is the first part of the product name. 

__COMPLETENESS__: If all of the k-mers in the gene are found in the reads, the completeness is noted as 'Full', otherwise if there are some k-mers missing, it is noted as 'Partial'. 

__%COVERAGE__: The percentage coverage is the number of underlying k-mers in the gene where at least 1 matching k-mer has been found in the reads. 100 indicates that every k-mer in the gene is covered. Low coverage results are not shown (controlled by the --min_perc_coverage parameter).

__ACCESSION__: This is the accession number from where the typing sequence originates. You can look this up at NCBI or EBI.

__DATABASE__: This is where the data has come from, which is currently always plasmidfinder.

__PRODUCT__: This is the full product of the gene as found in the database.

# Example usage
A real [test file](https://github.com/andrewjpage/tiptoft/raw/master/example_data/ERS654932_plasmids.fastq.gz) is bundled in the repository. Download it then run:

```
tiptoft ERS654932_plasmids.fastq.gz
```

The [expected output](https://raw.githubusercontent.com/andrewjpage/tiptoft/master/example_data/expected_output) is in the repository. This uses a bundled database; for the latest database, run the `tiptoft_database_downloader` script.

## Resource Usage

**Performance characteristics:**
- **Memory**: ~80 MB RAM for bundled database
- **Speed**: ~1 minute for 800 MB uncompressed FASTQ
- **Scaling**: Linear with input file size
- **Disk**: Minimal temporary space needed

**Tested on:**
- Oxford Nanopore MinION data (Salmonella)
- PacBio Sequel data (various bacteria)
- Read lengths: 1-50 kb
- Coverage: 5x - 100x

## Testing and Quality

TipToft has comprehensive test coverage:

```bash
# Run tests
python3 -m pytest tiptoft/tests/

# Check coverage
python3 -m pytest tiptoft/tests/ --cov=tiptoft --cov-report=html
```

**Current metrics:**
- **Test coverage**: 92% (54 tests)
- **Modules tested**: All core functionality
- **CI/CD**: Travis CI with automated testing

## License
TipToft is free software, licensed under [GPLv3](https://github.com/andrewjpage/tiptoft/blob/master/GPL-LICENSE).

## Feedback/Issues

Having problems? We're here to help:

1. **Check documentation**: Review the [User Guide](docs/USER_GUIDE.md)
2. **Search issues**: https://github.com/andrewjpage/tiptoft/issues
3. **Report bugs**: Open a new issue with:
   - TipToft version (`tiptoft --version`)
   - Command used
   - Error message
   - Input file characteristics
4. **Email support**: andrew.page@quadram.ac.uk

## Contribute to the Software

We welcome contributions! Here's how to get involved:

### Reporting Issues
- Use the [issue tracker](https://github.com/andrewjpage/tiptoft/issues)
- Include reproduction steps
- Provide example data if possible

### Contributing Code

We use [GitHub Flow](https://guides.github.com/introduction/flow/) for development:

1. **Fork** the repository
2. **Create** a feature branch (`git checkout -b feature/amazing-feature`)
3. **Make** your changes with:
   - Tests for new functionality
   - Updated documentation
   - Code comments
4. **Ensure** all tests pass and coverage ≥80%
5. **Commit** with clear messages
6. **Push** to your fork
7. **Submit** a Pull Request with description

### Development Setup

```bash
# Clone repository
git clone https://github.com/andrewjpage/tiptoft.git
cd tiptoft

# Install in development mode
pip3 install -e .

# Install test dependencies
pip3 install pytest pytest-cov

# Run tests
python3 -m pytest tiptoft/tests/
```

See [DEVELOPER_GUIDE.md](docs/DEVELOPER_GUIDE.md) for detailed development information.

### Contribution Guidelines

- Follow PEP 8 style guidelines
- Write comprehensive docstrings
- Add tests for new features
- Maintain test coverage ≥80%
- Update documentation
- Keep commits focused and atomic

## Acknowledgments

- **PlasmidFinder team** for the plasmid typing database
- **Martin Hunt** for ARIBA code contributions
- **Contributors** to TipToft development

## Related Tools

- **PlasmidFinder**: https://cge.cbs.dtu.dk/services/PlasmidFinder/
- **ARIBA**: https://github.com/sanger-pathogens/ariba
- **MOB-suite**: https://github.com/phac-nml/mob-suite
- **Unicycler**: https://github.com/rrwick/Unicycler

---

**Developed by:** Andrew J. Page and Torsten Seemann  
**Institution:** Quadram Institute Bioscience  
**Contact:** andrew.page@quadram.ac.uk  
**Repository:** https://github.com/andrewjpage/tiptoft

