"""
Database downloader for PlasmidFinder sequences.

This module downloads plasmid typing sequences from the PlasmidFinder
BitBucket repository and prepares them for use with TipToft.

Derived from ARIBA:
github.com/sanger-pathogens/ariba/blob/master/ariba/ref_genes_getter.py
bfd6cc9a409828ecc940b554164fabc8e6cc6b9a
Author: Martin Hunt
Copyright (c) 2015 - 2018 by Genome Research Ltd.
GNU GPL version 3
"""

import os
import shutil
import subprocess

import pyfastaq


class Error(Exception):
    """Custom exception for RefGenesGetter errors."""
    pass


class RefGenesGetter:
    """
    Downloads and processes plasmid typing database from PlasmidFinder.
    
    This class handles:
    1. Downloading FASTA files from PlasmidFinder BitBucket repository
    2. Combining multiple files into a single database
    3. Normalizing sequence IDs and creating metadata
    
    Attributes:
        ref_db (str): Name of the reference database ("plasmidfinder")
        verbose (bool): Whether to keep temporary files for debugging
        genetic_code (int): Genetic code for translation (11 = bacterial)
        max_download_attempts (int): Max retries for failed downloads
        sleep_time (int): Seconds to wait between retries
    """
    
    def __init__(self, verbose=False):
        """
        Initialize the database downloader.
        
        Args:
            verbose (bool): If True, keeps temporary download directory
                           for debugging purposes
        """
        self.ref_db = 'plasmidfinder'
        self.verbose = verbose
        self.genetic_code = 11  # Bacterial/Archaeal genetic code
        self.max_download_attempts = 3
        self.sleep_time = 2

    def _validate_fasta_file(self, filepath):
        """
        Validate that a downloaded file contains FASTA data.

        BitBucket returns HTTP 200 with an HTML body for missing files,
        so a successful download may still contain unusable content.

        Args:
            filepath (str): Path to the downloaded file

        Raises:
            Error: If the file cannot be read, is empty, or starts with '<'
                   (indicating an HTML error page rather than FASTA data)
        """
        try:
            with open(filepath, 'r') as fh:
                first_line = fh.readline()
        except OSError as e:
            raise Error('Could not read downloaded file {}: {}'.format(filepath, e))

        if not first_line:
            raise Error(
                'Downloaded file is empty: {}. '
                'Check the URL and try again.'.format(filepath)
            )

        if first_line.strip().startswith('<'):
            raise Error(
                'Downloaded file appears to be an HTML error page, not a FASTA file: {}. '
                'First line: {}. '
                'Check that the filename still exists on the PlasmidFinder '
                'repository.'.format(filepath, first_line.strip())
            )

    def _get_from_plasmidfinder(self, outprefix):
        """
        Download and process PlasmidFinder database files.
        
        Downloads multiple FASTA files from the PlasmidFinder BitBucket
        repository, combines them, normalizes sequence IDs, and creates
        a TSV metadata file.
        
        Args:
            outprefix (str): Prefix for output files (will create .fa and .tsv)
            
        Raises:
            Error: If unable to create or access temporary directory
            
        Note:
            Creates two output files:
            - <outprefix>.fa: Combined FASTA sequences
            - <outprefix>.tsv: Metadata with original IDs
        """
        outprefix = os.path.abspath(outprefix)
        final_fasta = outprefix + '.fa'
        final_tsv = outprefix + '.tsv'
        tmpdir = outprefix + '.tmp.download'
        current_dir = os.getcwd()

        # Create temporary download directory
        try:
            os.mkdir(tmpdir)
            os.chdir(tmpdir)
        except Error:
            raise Error('Error mkdir/chdir ' + tmpdir)

        # PlasmidFinder database files to download
        # Note: 'enterobacteriales.fsa' replaced 'enterobacteriaceae.fsa' (Issue #28)
        files_to_download = [
            'enterobacteriales.fsa',
            'Inc18.fsa',
            'NT_Rep.fsa',
            'Rep1.fsa',
            'Rep2.fsa',
            'Rep3.fsa',
            'RepA_N.fsa',
            'RepL.fsa',
            'Rep_trans.fsa',
        ]

        base_url = ('https://bitbucket.org/genomicepidemiology'
                     '/plasmidfinder_db/raw/master/')

        for filename in files_to_download:
            cmd = 'curl -o ' + filename + ' ' + base_url + filename
            print('Downloading data with:', cmd, sep='\n')
            subprocess.check_call(cmd, shell=True)
            self._validate_fasta_file(os.path.join(tmpdir, filename))

        # Combine downloaded files into final FASTA and TSV
        print('Combining downloaded fasta files...')
        fout_fa = pyfastaq.utils.open_file_write(final_fasta)
        fout_tsv = pyfastaq.utils.open_file_write(final_tsv)
        name_count = {}

        for filename in os.listdir(tmpdir):
            if filename.endswith('.fsa'):
                print('   ', filename)
                file_reader = pyfastaq.sequences.file_reader(
                    os.path.join(tmpdir, filename))
                for seq in file_reader:
                    original_id = seq.id
                    seq.id = seq.id.replace('_', '.', 1)  # Normalize ID
                    seq.seq = seq.seq.upper()

                    # De-duplicate sequence IDs by appending occurrence count
                    if seq.id in name_count:
                        name_count[seq.id] += 1
                        seq.id = seq.id + '.' + str(name_count[seq.id])
                    else:
                        name_count[seq.id] = 1

                    print(seq, file=fout_fa)
                    print(seq.id, '0', '0', '.', '.',
                          'Original name was ' + original_id,
                          sep='\t', file=fout_tsv)

        pyfastaq.utils.close(fout_fa)
        pyfastaq.utils.close(fout_tsv)
        print('\nFinished combining files\n')
        
        # Clean up temporary directory unless verbose mode
        os.chdir(current_dir)
        if not self.verbose:
            shutil.rmtree(tmpdir)
            
        # Print citation information
        print('Finished. Final files are:', final_fasta,
              final_tsv, sep='\n\t', end='\n\n')
        print('If you use this downloaded data, please cite:')
        print('"PlasmidFinder and pMLST: in silico detection and typing')
        print(' of plasmids", Carattoli et al 2014, PMID: 24777092\n')

    def run(self, outprefix):
        """
        Execute the database download.
        
        Args:
            outprefix (str): Prefix for output files
            
        Note:
            Currently only supports downloading from PlasmidFinder.
            The method name is constructed dynamically to allow for
            future support of additional databases.
        """
        getattr(self, '_get_from_' + self.ref_db)(outprefix)
