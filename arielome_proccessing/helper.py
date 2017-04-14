"""
helper.py

helper functions and variables for: arielome_pipe.py

"""

import os
import subprocess

# variables ---------------------------------------------------------------


out_fastq = {
            'fish': ('fish_01.fq', 'fish_02.fq'),
            'fly': ('fly_01.fq', 'fly_02.fq'),
            'minigene': ('minigene_01.fq', 'minigene_02.fq'),
            'planaria': ('planaria_01.fq', 'planaria_02.fq'),
            'pombe': ('pombe_01.fq', 'pombe_02.fq')
            }



# functions ---------------------------------------------------------------
# #### --------------------------------------------------------------------

# Initializer -------------------------------------------------------------

def check_file(file_name):
    """
    check whether file exist, if not raises and exception
    """
    if not os.path.isfile(file_name):
        raise NameError('file: %s not found' % file_name)


# BarcodeSpliting ---------------------------------------------------------

def barcode_split_cmd(read, out_prefix):
    """
    returns the cmd of the barcode split task
    """
    barcodes_file = "data/barcodes_spitting.txt"
    command_line = [
                   'cat', read, '|',
                   'fastx_barcode_splitter.pl',
                   "--bcfile", barcodes_file,
                   "--prefix", out_prefix,
                   "--suffix", "_01.fq",
                   "--bol"
                   ]
    return ' '.join(_ for _ in command_line)


def get_r2_from_subset(R2, out_prefix_dir):
    """
    given a fastq file r1 that is a subset from R1
    extracts the corresponding mates r2 from R2
    """

    def get_r2_command(r1_input, r2_output):
        command_line = [
                      "bioawk -cfastx '{print $name}'",
                      out_prefix_dir + r1_input,
                      "|",
                      "seqtk subseq",
                      R2, "-",
                      ">", out_prefix_dir + r2_output
                      ]
        return ' '.join(_ for _ in command_line)

    for specie in out_fastq:
        # check if file exist
        if not os.path.isfile(out_prefix_dir + out_fastq[specie][0]):
            raise NameError('file: %s not found' % out_fastq[specie][0])
        yield get_r2_command(out_fastq[specie][0], out_fastq[specie][1])

