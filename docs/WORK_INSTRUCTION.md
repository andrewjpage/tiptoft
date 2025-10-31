# TipToft Work Instruction

## Standard Operating Procedure for Plasmid Detection

This document provides step-by-step instructions for routine use of TipToft in a production environment.

## Scope

This procedure applies to:
- Bacterial genome sequencing projects
- Quality control of assemblies
- Plasmid characterization studies
- Antimicrobial resistance profiling

## Prerequisites

### Required Software
- TipToft v1.0.0 or later
- Python 3.6 or later
- 8 GB RAM minimum
- 20 GB disk space for databases

### Required Skills
- Basic command-line proficiency
- Understanding of FASTQ file format
- Knowledge of plasmid biology

## Procedure

### 1. Environment Setup

#### 1.1 Verify Installation

```bash
# Check TipToft is installed
tiptoft --version

# Expected output: tiptoft, version x.y.z
```

#### 1.2 Create Working Directory

```bash
# Create project directory
mkdir -p /path/to/project/tiptoft_analysis
cd /path/to/project/tiptoft_analysis

# Create subdirectories
mkdir -p input output logs
```

#### 1.3 Database Verification

```bash
# Check bundled database exists
tiptoft --help | grep "defaults to bundled database"

# For latest database (optional):
tiptoft_database_downloader output/plasmidfinder_db
```

**Note**: Update database monthly for latest plasmid types.

### 2. Input Data Preparation

#### 2.1 Verify Input Files

```bash
# Check FASTQ file exists and is readable
ls -lh input/sample.fastq.gz

# Verify file is not corrupted (optional)
gunzip -t input/sample.fastq.gz
```

**Quality Criteria**:
- File size >100 KB (indicates data present)
- Compression intact (no errors from gunzip -t)
- Read count >1000 recommended

#### 2.2 Input File Naming Convention

Use consistent naming:
```
[SampleID]_[SeqType]_[Date].fastq.gz

Examples:
- Sample001_Nanopore_20240101.fastq.gz
- Sample002_PacBio_20240115.fastq.gz
```

### 3. Analysis Execution

#### 3.1 Standard Analysis (Nanopore)

```bash
# Run with default parameters for Nanopore
tiptoft \
    -o output/Sample001_results.txt \
    -f output/Sample001_plasmid_reads.fastq \
    --verbose \
    input/Sample001_Nanopore_20240101.fastq.gz \
    2>&1 | tee logs/Sample001.log
```

**Parameters explained**:
- `-o`: Output file for results
- `-f`: Save matching reads
- `--verbose`: Detailed logging
- `2>&1 | tee`: Save log to file

#### 3.2 PacBio Analysis

```bash
# Adjust k-mer for PacBio's lower error rate
tiptoft \
    -k 15 \
    -o output/Sample002_results.txt \
    -f output/Sample002_plasmid_reads.fastq \
    input/Sample002_PacBio_20240115.fastq.gz \
    2>&1 | tee logs/Sample002.log
```

#### 3.3 High-Throughput Batch Processing

```bash
# Process multiple samples
for sample in input/*.fastq.gz; do
    basename=$(basename "$sample" .fastq.gz)
    echo "Processing $basename..."
    
    tiptoft \
        -o "output/${basename}_results.txt" \
        -f "output/${basename}_plasmid_reads.fastq" \
        "$sample" \
        2>&1 | tee "logs/${basename}.log"
done
```

### 4. Quality Control

#### 4.1 Verify Output Files

```bash
# Check all expected outputs exist
ls -lh output/Sample001_results.txt
ls -lh output/Sample001_plasmid_reads.fastq
ls -lh logs/Sample001.log

# Verify results file is not empty
wc -l output/Sample001_results.txt
```

**Expected**:
- Results file: >0 lines (or 0 if no plasmids)
- Log file: Contains "Starting analysis"
- Plasmid reads: Present if plasmids detected

#### 4.2 Review Log for Errors

```bash
# Check for errors in log
grep -i "error" logs/Sample001.log
grep -i "warning" logs/Sample001.log

# Check completion
grep -i "finished" logs/Sample001.log
```

**Action if errors found**:
1. Review error message
2. Check input file integrity
3. Verify disk space available
4. Consult troubleshooting section

### 5. Results Interpretation

#### 5.1 Review Results File

```bash
# Display results
column -t output/Sample001_results.txt
```

#### 5.2 Classification Guidelines

| Coverage | Completeness | Confidence | Action |
|----------|--------------|------------|--------|
| 100% | Full | High | Confirm plasmid present |
| 95-99% | Full/Partial | High | Likely present, possible variant |
| 85-94% | Partial | Medium | Verify with additional analysis |
| <85% | Partial | Low | Not reported by default |

#### 5.3 Document Results

Create summary report:

```bash
cat > output/Sample001_summary.txt << EOF
Sample: Sample001
Date: $(date +%Y-%m-%d)
Analysis: TipToft v$(tiptoft --version | cut -d' ' -f3)

Plasmids Detected:
$(cat output/Sample001_results.txt)

Notes:
- [Add any observations]
- [Add any quality concerns]

Analyst: [Name]
EOF
```

### 6. Follow-up Actions

#### 6.1 If Plasmids Detected

**Option A: Assembly Verification**
```bash
# Check if plasmids are in assembly
# (requires assembly FASTA and BLAST+)
makeblastdb -in assembly.fasta -dbtype nucl
blastn -query output/Sample001_plasmid_reads.fastq \
       -db assembly.fasta -outfmt 6
```

**Option B: Targeted Assembly**
```bash
# Assemble just the plasmid reads
# Example with Flye:
flye --nano-raw output/Sample001_plasmid_reads.fastq \
     --out-dir output/Sample001_plasmid_assembly \
     --threads 8
```

#### 6.2 If No Plasmids Detected

**Troubleshooting steps**:

1. Verify read quality:
```bash
# Check read statistics (requires seqtk or similar)
seqtk fqchk input/Sample001_Nanopore_20240101.fastq.gz
```

2. Try lower k-mer:
```bash
tiptoft -k 11 -c 70 -o output/Sample001_sensitive.txt \
    input/Sample001_Nanopore_20240101.fastq.gz
```

3. Check database:
```bash
# Verify database file exists and is not empty
ls -lh /path/to/plasmid_data.fa
```

### 7. Data Management

#### 7.1 Archive Results

```bash
# Create archive
tar -czf Sample001_tiptoft_$(date +%Y%m%d).tar.gz \
    output/Sample001_* \
    logs/Sample001.log

# Move to archive location
mv Sample001_tiptoft_*.tar.gz /path/to/archive/
```

#### 7.2 Clean Up Temporary Files

```bash
# Remove large intermediate files if not needed
# ONLY if results are verified and archived
rm output/Sample001_plasmid_reads.fastq
```

### 8. Reporting

#### 8.1 Generate Summary Report

```bash
# Create CSV summary for all samples
echo "Sample,Date,Plasmids_Detected,Full_Coverage,Partial_Coverage" > summary.csv

for result in output/*_results.txt; do
    sample=$(basename "$result" _results.txt)
    date=$(date +%Y-%m-%d)
    total=$(wc -l < "$result")
    full=$(grep -c "Full" "$result" || echo 0)
    partial=$(grep -c "Partial" "$result" || echo 0)
    
    echo "$sample,$date,$total,$full,$partial" >> summary.csv
done
```

#### 8.2 Quality Metrics

Document these metrics:
- Total samples processed
- Samples with plasmids detected
- Average coverage percentage
- Analysis time per sample
- Any failed analyses

### 9. Troubleshooting

#### Common Issues and Solutions

| Issue | Cause | Solution |
|-------|-------|----------|
| No output file created | Input file not found | Verify input path |
| Empty results | No plasmids present | Try -k 11 -c 70 |
| "Cannot access input file" | Wrong path/permissions | Check file exists and readable |
| Killed/Out of memory | Insufficient RAM | Use smaller chunks or more RAM |
| Very slow processing | Large file | Normal; use -p 1000 for progress |

#### 9.1 Performance Issues

If analysis is too slow:
```bash
# Increase thresholds for speed
tiptoft \
    -k 15 \
    -m 12 \
    --min_kmers_for_onex_pass 8 \
    --min_block_size 100 \
    input/sample.fastq.gz
```

#### 9.2 Sensitivity Issues

If missing expected plasmids:
```bash
# Increase sensitivity (slower)
tiptoft \
    -k 11 \
    -c 70 \
    -m 5 \
    --min_kmers_for_onex_pass 3 \
    input/sample.fastq.gz
```

### 10. Validation and QC

#### 10.1 Positive Control

Run known positive sample monthly:
```bash
# Download example data
wget https://github.com/andrewjpage/tiptoft/raw/master/example_data/ERS654932_plasmids.fastq.gz

# Run analysis
tiptoft ERS654932_plasmids.fastq.gz > control_results.txt

# Compare to expected
wget https://raw.githubusercontent.com/andrewjpage/tiptoft/master/example_data/expected_output
diff control_results.txt expected_output
```

**Expected**: Should match expected_output file

#### 10.2 Negative Control

Process a sample known to have no plasmids:
- Should produce empty or minimal results
- Verifies no contamination

### 11. Documentation Requirements

For each analysis, document:
1. Sample ID and source
2. Date of analysis
3. TipToft version used
4. Parameters used (if non-default)
5. Results summary
6. Any anomalies or issues
7. Follow-up actions taken
8. Analyst name and date

### 12. Safety and Compliance

- Follow institutional biosafety guidelines
- Maintain data privacy per regulations
- Secure computational resources
- Back up results regularly
- Document analysis parameters for reproducibility

## References

- TipToft User Guide: docs/USER_GUIDE.md
- TipToft Paper: https://doi.org/10.21105/joss.01021
- PlasmidFinder: https://cge.cbs.dtu.dk/services/PlasmidFinder/

## Revision History

| Version | Date | Changes | Author |
|---------|------|---------|--------|
| 1.0 | 2024 | Initial SOP | TipToft Team |

## Contact

For technical support:
- GitHub Issues: https://github.com/andrewjpage/tiptoft/issues
- Email: andrew.page@quadram.ac.uk
