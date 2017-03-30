#!/usr/bin/python3

## helper function to count dicodons and write counts to files
from collections import defaultdict

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


def iter_dicodons(dna):
    """
    iterates over the dicodones of a dna string
    """
    for i in range(0, len(dna) - 5, 3):
        yield dna[i: i + 6]


def count_dicodons(dna):
    """
    counts the dicodon composition of a dna string
    """
    counts_dict = defaultdict(int)

    for dicodon in iter_dicodons(dna):
        counts_dict[dicodon] += 1
    return counts_dict

def write_counts_to_file(seq_id, dna, file_to_writte):
    """
    write dicodon counts to file
    """
    counts = count_dicodons(dna)
    counts = [counts[x] for x in generate_dicodones()]
    seq_len = len(dna)
    with open(file_to_writte, "a") as myfile:
        myfile.write(seq_id + '\t')
        myfile.write('\t' + str(seq_len) + '\t')
        myfile.write('\t'.join([str(x) for x in counts]))
        myfile.write('\n')






