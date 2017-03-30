from multiprocessing import Pool 
import sys

# Script usr parameters
try:
    scrp_name, input_file, number_of_codons = sys.argv
except ValueError:
    print('Invalid arguments')
    print('usage: python script.py <input_file> <number_of_codons>')
    sys.exit()

out_file = "tmp_seqs"
## begining of functions
def read_transcripts(file_name):
    """
    Reads the transcripts, all transcritps should be equal length
    Args:
        file_name, str path to input files
    Returns:
        iterators over the lines of file_name
    """
    try:
        itr = (x for x in open(file_name))
        next(itr) # skeep header
        return itr
    except FileNotFoundError:
        print("unable to open file: %s" % file_name)

def get_Sequence_from_line(line):
    """
    given a sequence retrives the line
    Args:
        line, str output from previous function
    Reutrns:
        sequence, str
    """
    return line.split()[1]

def get_codons_by_pos(sequence):
    """
    gets the codons in 1st frame of sequence
    """
    return (sequence[i: i + 3] for i in range(0, len(sequence) - 2, 3))

def write_by_pos(sequence):
    """ write codons in sequence by position to file """
    sequence = get_Sequence_from_line(sequence)
    with open(out_file, 'a') as my_file:
        my_file.write('\t'.join(x for x in get_codons_by_pos(sequence)))
        my_file.write('\n')

def make_out_file():
    oo = open(out_file, 'w')
    oo.close()
## end of functions

print("sequences by codons positions running on file: %s" % input_file)
make_out_file()
p = Pool(35)
p.map(write_by_pos, read_transcripts(input_file))

pos = open('positions', 'w')
number_of_codons = int(number_of_codons) * 2 # begining and end positions

for i in range(number_of_codons):
    pos.write(str(i + 1))
    pos.write('\n')
pos.close()
print("script completed")
