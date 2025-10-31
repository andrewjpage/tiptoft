# TipToft Documentation

This directory contains comprehensive documentation for TipToft, a tool for plasmid detection in uncorrected long-read sequencing data.

## Documentation Structure

### [USER_GUIDE.md](USER_GUIDE.md)
Complete user documentation covering:
- Installation instructions for multiple platforms
- Basic and advanced usage examples
- Parameter descriptions and tuning
- Workflow examples
- Troubleshooting guide
- Output interpretation
- Citations

**Audience**: End users, bioinformaticians, lab technicians

### [DEVELOPER_GUIDE.md](DEVELOPER_GUIDE.md)
Technical documentation for developers:
- Architecture overview and design
- Module descriptions and interactions
- Algorithm details
- Development setup
- Testing guidelines
- Performance optimization
- Contributing guidelines
- Debugging techniques

**Audience**: Software developers, contributors, maintainers

### [WORK_INSTRUCTION.md](WORK_INSTRUCTION.md)
Standard operating procedure for routine use:
- Step-by-step procedures
- Quality control checkpoints
- Batch processing workflows
- Data management practices
- Validation procedures
- Safety and compliance
- Reporting templates

**Audience**: Lab operators, QC staff, production environments

### [API_REFERENCE.md](API_REFERENCE.md)
Complete API documentation:
- Class and method signatures
- Parameter descriptions
- Return types
- Usage examples
- Data structures
- Error handling

**Audience**: Python developers, API users, automation developers

## Quick Links

- **Main Repository**: https://github.com/andrewjpage/tiptoft
- **Issues**: https://github.com/andrewjpage/tiptoft/issues
- **PyPI**: https://pypi.org/project/tiptoft/
- **Docker Hub**: https://hub.docker.com/r/andrewjpage/tiptoft
- **Paper**: https://doi.org/10.21105/joss.01021

## Quick Start

For most users, start with:
1. [USER_GUIDE.md](USER_GUIDE.md) - Installation and basic usage
2. Run example: `tiptoft example_data/ERS654932_plasmids.fastq.gz`

For developers:
1. [DEVELOPER_GUIDE.md](DEVELOPER_GUIDE.md) - Architecture and setup
2. [API_REFERENCE.md](API_REFERENCE.md) - Detailed API docs

For production environments:
1. [WORK_INSTRUCTION.md](WORK_INSTRUCTION.md) - SOPs and workflows

## Getting Help

- **Documentation Questions**: Review appropriate guide above
- **Usage Questions**: See USER_GUIDE.md troubleshooting section
- **Bug Reports**: https://github.com/andrewjpage/tiptoft/issues
- **Feature Requests**: https://github.com/andrewjpage/tiptoft/issues
- **Email Support**: andrew.page@quadram.ac.uk

## Contributing to Documentation

Documentation improvements are welcome! To contribute:

1. Fork the repository
2. Edit documentation files (Markdown format)
3. Test any code examples
4. Submit a pull request

Documentation should be:
- Clear and concise
- Include examples where helpful
- Use consistent formatting
- Be technically accurate

## License

TipToft and its documentation are licensed under GNU GPL v3.

## Citations

If you use TipToft, please cite:

**TipToft**:
```
AJ Page, T Seemann (2019). TipToft: detecting plasmids contained in 
uncorrected long read sequencing data. Journal of Open Source Software, 
4(35), 1021, https://doi.org/10.21105/joss.01021
```

**PlasmidFinder**:
```
Carattoli et al, In Silico Detection and Typing of Plasmids using 
PlasmidFinder and Plasmid Multilocus Sequence Typing, 
Antimicrob Agents Chemother. 2014;58(7):3895–3903.
```
