"""
File:   read_genome.py
Description:    Reads a Fasta file using Bio package, and creates a dictionary of the genome in memory
Author: Óscar González-Velasco
Date:   1 - August - 2022
"""

from Bio import SeqIO
import gzip
# from nucleotide_to_int import nucleotid_to_int

def read_genome(file_path: str):
    """
    Index and access fastas like a database using the Python dictionary data type (like a hash in Perl). This is very useful for moderately large files where you only need to access certain elements of the file, and makes for a nice quick ’n dirty database. 
    You can use the function Bio.SeqIO.to_dict() to make a SeqRecord dictionary (in memory). By default this will use each record’s identifier (i.e. the .id attribute) as the key.
    
    :param file_path: string, path and name of the fasta file to be loaded
    :returns: a dictionary with the fasta entries
    :raises keyError: raises an exception    
    """
    in_handle = gzip.open(file_path, "rt")
    genome_dict = SeqIO.to_dict(SeqIO.parse(in_handle, "fasta"))
    in_handle.close()
    #list(genome_dict.keys())
    return(genome_dict)
    #nucleotid_to_int(seq_record.seq)
