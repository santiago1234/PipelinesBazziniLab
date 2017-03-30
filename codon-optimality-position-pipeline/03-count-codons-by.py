from multiprocessing import Pool
from count_codons_in_seq import count_codons, codons

out_file = "tmp_counts"

def read_sequence_from_file():
    """Reads the sequences in seqs_tmp"""
    return (line for line in open('seqs_tmp'))

def procces_line(line):
    """
    procces a line of seqs_tmp
    """
    return line.split()

def make_out_file():
    f = open(out_file, 'w')
    f.close()

def write_count_table(line):
    """
    counts the codons and writtes the result to a file
    """
    position, seq = procces_line(line)
    codon_counts = count_codons(seq)
    with open(out_file, 'a') as out:
        for codon in codons():
            to_write = '\t'.join([str(position), codon, str(codon_counts[codon])])
            out.write(to_write)
            out.write('\n')


make_out_file()
p = Pool(90)
p.map(write_count_table, read_sequence_from_file())
print("counting codons completed")

