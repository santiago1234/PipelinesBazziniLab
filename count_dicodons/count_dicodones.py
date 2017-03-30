 #!/usr/bin/python3

from multiprocessing import Pool
from helper import write_counts_to_file, generate_dicodones
import sys

# functions ---------------------------------------------------------------

def get_sequences(file_name):
    """
    procces the input file to get the sequences
    for counting the codons
    """
    i = 0 # counter to skip firs line (file header)
    for line in open(file_name, 'r'):
        if i > 0:
            data = line.split()
            yield [data[0], data[1]]
        i += 1


def make_out_file(out_file_name):
    """
    creates the out filename and adds the header
    """
    dicodons = [x for x in generate_dicodones()]
    header = ["id", "len"] + dicodons
    handler = open(out_file_name, "w")
    handler.write('\t'.join(x for x in header))
    handler.write('\n')
    handler.close()


# parse script arguments --------------------------------------------------

if len(sys.argv) != 3:
    print('too few arguments 2 needed')
    print('count_dicodons.py -i input_library -o output_counts')
    sys.exit(1)

input_library = sys.argv[1]
output_counts = sys.argv[2]

# run counts --------------------------------------------------------------

print('runing dicodon counts on file: ', input_library, ' ...')
make_out_file(output_counts)

def run_dicodon_counts(seq_line):
    write_counts_to_file(seq_line[0], seq_line[1], output_counts)


p = Pool(13)
p.map(run_dicodon_counts, get_sequences(input_library))

print('job finished!!!')

