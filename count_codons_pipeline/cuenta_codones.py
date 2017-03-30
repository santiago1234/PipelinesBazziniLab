#!/usr/bin/python

###
# Count the codont frequency in a DNA string
# given a DNA string counts the codon frequency on the first 
# frame
## #


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
    Counts the codon frquency of dna_seq
    args:
        dna_seq: a str of uppercase nucleotides
    returns:
        a dicctionary of the codon frequency
    """
    tripletes = iterate_over_tripletes(dna_seq)
    codon_frequency = defaultdict(int)
    for triplete in tripletes:
        codon_frequency[triplete] += 1

    transcript_length = int(len(dna_seq) / 3)
    for k, v in codon_frequency.items():
        codon_frequency[k] = v / transcript_length
    return codon_frequency


def write_codon_counts(dna_seq, file_name, id_transcrit):
    """
    Writtes codon density / frquency to file_name
    """
    codon_frequency = count_codons(dna_seq)

    def print_freq():
        s = [codon_frequency[c] for c in codons()]
        return '\t'.join(str(_) for _ in s)

    with open(file_name, 'a') as file:
        file.write(id_transcrit + '\t')
        file.write(print_freq())
        file.write('\t' + str(len(dna_seq)))
        file.write('\n')


