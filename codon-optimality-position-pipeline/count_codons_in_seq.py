#!/usr/bin/python
from collections import defaultdict


def codons():
    """
    Returns a list of the 64 codons
    """
    codons = ['TTT', 'TTC', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG',
        'ATT', 'ATC', 'ATA', 'ATG', 'GTT', 'GTC', 'GTA', 'GTG',
        'TCT', 'TCC', 'TCA', 'TCG', 'CCT', 'CCC', 'CCA', 'CCG',
        'ACT', 'ACC', 'ACA', 'ACG', 'GCT', 'GCC', 'GCA', 'GCG',
        'TAT', 'TAC', 'TAA', 'TAG', 'CAT', 'CAC', 'CAA', 'CAG',
        'AAT', 'AAC', 'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 'GAG',
        'TGT', 'TGC', 'TGA', 'TGG', 'CGT', 'CGC', 'CGA', 'CGG',
        'AGT', 'AGC', 'AGA', 'AGG', 'GGT', 'GGC', 'GGA', 'GGG']
    return codons


def iterate_over_tripletes(dna_seq):
    """
    Iterates over the codons of a DNA sequence
    args:
        dna_seq: a string of uppercase nucleotides
    returns:
        an iterator over the triples in the first frame
        of the dna_seq
    """
    return (dna_seq[i: i + 3] for i in range(0, len(dna_seq) -3 + 1, 3))


def count_codons(dna_seq):
    """
    count the codons in first frame on dna_seq
    Args:
        dna_seq, str dna sequence to count codons on first frame
    Returns:
        diccionary mapping codons to frequencies
    """
    codon_frequency = defaultdict(int)
    def add_count(triplet):
        codon_frequency[triplet] += 1
    [add_count(x) for x in iterate_over_tripletes(dna_seq)]
    return codon_frequency

def print_counts(dna_seq, position):
    counts = count_codons(dna_seq)
    for codon in codons():
        out = '\t'.join([str(position), codon, str(counts[codon])])
        print(out)

