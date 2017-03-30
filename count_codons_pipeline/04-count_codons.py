from multiprocessing import Pool 
from cuenta_codones import write_codon_counts, codons
import sys

# Script usr parameters
try:
    scrp_name, input_file, output_file, taskings = sys.argv
except ValueError:
    print('Invalid arguments')
    print('usage: python script.py <input_file> <output_file> <taskings>')
    sys.exit()

##### Functions        #########
def read_input_file(file_name):
    """ Reads the input file with no premature stop codons,
    Args:
        file_name: str of path to input file
    Returns:
        Iterator of input file
    """
    try:
        return (x for x in open(file_name))
    except FileNotFoundError:
        print("unable to open file: %s" % file_name)

def proccess_line(file_line):
    """
    Process a line of the input file
    Args:
        file_line: a str representing a single line of the input file
    Returns:
        a list with the following information
        [id, seq, length]
    """
    file_line = file_line.split()
    return [file_line[0], file_line[1], file_line[3]]


def make_out_file(out_filename):
    """
    Creates the output file for the codons and adds
    the columns header
    """
    fo = open(out_filename, 'w')
    print("Writting codon frequncy/density table to file: %s" % fo.name)
    header = ["transcript_id"] + codons() + ["transcript_len"]
    header = '\t'.join(_ for _ in header)
    fo.write(header)
    fo.write('\n')
    fo.close()


def run_codon_count(file_line):
    """
    main function for counting the codons
    in each transcit and writing to output file
    """
    id_transc, seq_transc = proccess_line(file_line)[:2]
    write_codon_counts(seq_transc, output_file, id_transc)

##### End of functions #########
make_out_file(output_file)
p = Pool(int(taskings))
p.map(run_codon_count, read_input_file(input_file))
