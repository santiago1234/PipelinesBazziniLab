#!/usr/bin/python3

""""
count_k_mers.py
count six-mers of seqs in fasta file

usage:
    python count_k_mers.py <in_fasta_file>
"""


from collections import defaultdict
from Bio import SeqIO
import sys

# function to count the k-mers --------------------------------------------

def generate_codones():
    """ generates the codons"""
    nuc = ['A', 'C', 'G', 'T']
    codons = list()
    for x in nuc:
        for y in nuc:
            for z in nuc:
                codons.append(x + y + z)
    return codons


def generate_dicodones():
    for x in generate_codones():
        for y in generate_codones():
            yield x + y


def iterate_over_kmer(dna, k):
    """
    returns a generator that iterates over the kmers
    """

    for i in range(0, len(dna) - k + 1):
        yield dna[i: i + k]

def count_k_mers(dna, k):
    """
    counts the k-mers in dna
    Returns:
        a counter dict with counts
    """
    counts = defaultdict(int)
    for kmer in iterate_over_kmer(dna, k):
        counts[kmer.upper()] += 1
    return [counts[x] for x in generate_dicodones()]


# code --------------------------------------------------------------------

fasta_file = sys.argv[1]
out_counts = open('counts.tab', 'w') 
header = ['seq_id'] + [x for x in generate_dicodones()]
out_counts.write('\t'.join(x for x in header))
out_counts.write('\n')

for record in SeqIO.parse(fasta_file, 'fasta'):
    print('counting %s ...' % record.id)
    row = [record.id]
    counts = count_k_mers(str(record.seq), 6)
    row = row + counts 
    out_counts.write('\t'.join(str(x) for x in row) + '\n')

out_counts.close()
