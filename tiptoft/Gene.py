"""
Gene representation and coverage analysis.

This module contains the Gene class which represents a plasmid typing sequence
(gene/allele) and tracks its k-mer coverage in the sequencing reads.
"""
import re


class Gene:
    """
    Represents a plasmid typing sequence with coverage information.
    
    A Gene object tracks both the k-mers that were found in the reads
    (kmers_with_coverage) and those that were not found (kmers_without_coverage).
    This allows calculation of completeness and percentage coverage.
    
    Attributes:
        name (str): Full name of the gene from the database
        kmers_with_coverage (int): Number of k-mers found in reads
        kmers_without_coverage (int): Number of k-mers not found in reads
    """
    
    def __init__(self, name, kmers_with_coverage, kmers_without_coverage):
        """
        Initialize a Gene with coverage information.
        
        Args:
            name (str): Full gene name from the plasmid database
            kmers_with_coverage (int): Count of k-mers detected in reads
            kmers_without_coverage (int): Count of k-mers not detected in reads
        """
        self.name = name
        self.kmers_with_coverage = kmers_with_coverage
        self.kmers_without_coverage = kmers_without_coverage

    def __str__(self):
        """
        Generate tab-delimited string representation for output.
        
        Returns:
            str: Tab-separated values: short_name, completeness, coverage%, 
                 accession, database, full_name
        """
        return "\t".join((
            self.short_name(),
            self.completeness(),
            str(self.percentage_coverage()),
            self.accession(),
            'plasmidfinder',
            self.name))

    def completeness(self):
        """
        Determine if gene has full or partial k-mer coverage.
        
        Returns:
            str: "Full" if all k-mers were found, "Partial" otherwise
        """
        if self.is_full_coverage():
            return "Full"
        else:
            return "Partial"

    def is_full_coverage(self):
        """
        Check if all k-mers in the gene were detected.
        
        Returns:
            bool: True if no k-mers are missing, False otherwise
        """
        if self.kmers_without_coverage == 0:
            return True
        else:
            return False

    def percentage_coverage(self):
        """
        Calculate the percentage of k-mers that were detected.
        
        Coverage is calculated as (kmers_with_coverage / total_kmers) * 100.
        This represents what proportion of the gene sequence was found in the reads.
        
        Returns:
            int: Percentage coverage (0-100)
        """
        total_kmers = self.kmers_with_coverage + self.kmers_without_coverage
        if total_kmers > 0:
            return int((self.kmers_with_coverage*100)/total_kmers)
        else:
            return 0

    def prefix_short_name(self):
        """
        Extract the prefix of the short name (text before first dot).
        
        Returns:
            str: Prefix of short name, or empty string if not found
            
        Example:
            For "rep7.1_repC_AB037671", returns "rep7"
        """
        regex = r"^([^\.]+)\."

        m = re.search(regex, self.short_name())
        if m and m.group:
            return str(m.group(1))
        else:
            return ''

    def short_name(self):
        """
        Extract the human-readable short name from the full gene name.
        
        The short name is the first field before the underscore in the gene name.
        
        Returns:
            str: Short name, or empty string if pattern doesn't match
            
        Example:
            For "rep7.1_repC_AB037671", returns "rep7.1"
        """
        regex = r"^([^_]+)_"

        m = re.search(regex, self.name)
        if m and m.group:
            return str(m.group(1))
        else:
            return ''

    def accession(self):
        """
        Extract the accession number from the gene name.
        
        The accession is embedded in the database sequence name and refers
        to the original NCBI/EBI accession number.
        
        Returns:
            str: Accession number, or empty string if pattern doesn't match
            
        Example:
            For "rep7.1_repC_AB037671", returns "AB037671"
        """
        regex = r"^([^_]+)_([^_]*)_(.+)$"

        m = re.search(regex, self.name)
        if m and m.group:
            return str(m.group(3))
        else:
            return ''
