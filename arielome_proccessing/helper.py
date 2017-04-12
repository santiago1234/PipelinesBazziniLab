"""
hlper.py

helper functions and variables for: arielome_pipe.py

"""

import os

# functions ---------------------------------------------------------------

def check_file(file_name):
    """
    check whether file exist, if not raises and exception
    """
    if not os.path.isfile(file_name):
        raise NameError('file: %s not found' % file_name)


def barcode_split_cmd(read, out_prefix, r1 = True):
    """
    returns the cmd of the barcode split task
    """
    barcodes_file = "data/barcodes/barcodes_r1.txt" if r1 else "data/barcodes/barcodes_r2.txt"
    command_line = [
                   'fastx_barcode_splitter.pl', read,
                   "--bcfile", barcodes_file,
                   "--prefix", out_prefix,
                   "--suffix", ".fq",
                   "--bol"
                   ]
    return ' '.join(_ for _ in command_line)


