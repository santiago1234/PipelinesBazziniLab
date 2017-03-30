from multiprocessing import Pool
from find_premature_stop import find_premature_stop
import sys

# Script user parameters
try:
    script_name, input_file, output_file, taskings = sys.argv
except ValueError:
    print('Invalid arguments')
    print('usage: python detect_premature_stop.py <input_file> <output_file> <taskings>')


###### FUNCTIONS ########

def read_transcrits(file_name):
    """
    Generates an iterator over file_name
    Args:
        file_name: str, path to file name
    Returns:
        an iterator over the lines of file_name
    """
    try:
        return (_ for _ in open(file_name))
    except FileNotFoundError:
        print("unable to open file: %s" % file_name)


def read_sequence(file_line):
    """
    reads a line of input file
    Args:
        file_line: a str representing a line of input file
    Returns:
        tuple: (id, sequence)
    """
    id_transcrit = file_line.split()[0]
    seq = file_line.split()[1]
    return [id_transcrit, seq]


def make_out_file(output_file):
    """
    Creates the output file for the script
    """
    output_f = open(output_file, 'w')
    header = ["transcript_id", "transcrit_sequence", "stop_position",
    "transcrit_leng"]
    header = '\t'.join(_ for _ in header)
    output_f.write(header)
    output_f.write('\n')
    output_f.close()


def detect_premature_stop(file_line):
    """"
    main function to detect the position
    of the premature stop codon and write to
    output file
    """
    id_transcrit, seq = read_sequence(file_line)
    stop_position = find_premature_stop(seq)
    data = (id_transcrit, seq, stop_position, len(seq))
    with open(output_file, 'a') as file_out:
        file_out.write('\t'.join(str(_) for _ in data))
        file_out.write('\n')


###### END FUNCTIONS ####


def main():
    make_out_file(output_file)
    p = Pool(int(taskings))
    p.map(detect_premature_stop, read_transcrits(input_file))


main()
print("script finished >>>")

