"""
helper.py

helper functions and variables for: arielome_pipe.py

"""

import os
import subprocess

# variables ---------------------------------------------------------------


fastq_files = {
            'fish': ('fish_01.fq', 'fish_02.fq'),
            'fly': ('fly_01.fq', 'fly_02.fq'),
            'minigene': ('minigene_01.fq', 'minigene_02.fq'),
            'planaria': ('planaria_01.fq', 'planaria_02.fq'),
            'pombe': ('pombe_01.fq', 'pombe_02.fq')
            }


adapters = {
           "fish": ("TCTAACGGCGAAATGGC", "TTAGTCACCTA"),
           "fly": ("CGCTGCATTGAGGTGCCATGGC", "ACTAGTTAGTCA"),
           "pombe": ("GTGGTTACCGACAATGGC", "CGGTCAGCTACTTA"),
           "planaria": ("AAACCGTTTCATCATGGC", "CGGTCAGCTACTTA"),
           "minigene": ("CCTGCAGGCACCATG", "CTA")
           }


trimmed_fastq = {
                "fish": ("fish_01_trimmed.fq", "fish_02_trimmed.fq"),
                "fly": ("fly_01_trimmed.fq", "fly_02_trimmed.fq"),
                "minigene" :("minigene_01_trimmed.fq", "minigene_02_trimmed.fq"),
                "planaria": ("planaria_01_trimmed.fq", "planaria_02_trimmed.fq"),
                "pombe": ("pombe_01_trimmed.fq", "pombe_02_trimmed.fq")
                }


references = {
             "fish": "data/transcriptome_references/fish",
             "fly": "data/transcriptome_references/fly",
             "planaria": "data/transcriptome_references/planaria",
             "pombe": "data/transcriptome_references/pombe"
             }

bam_files = {
            "fish": ("fish.bam", "fish_filter.bam"),
            "fly": ("fly.bam", "fly_filter.bam"),
            "planaria": ("planaria.bam", "planaria_filter.bam"),
            "pombe": ("pombe.bam", "pombe_filter.bam")
            }

fasta_bed = {
            "fish": ("fish.bed", "fish.fasta"),
            "fly": ("fly.bed", "fly.fasta"),
            "planaria": ("planaria.bed", "planaria.fasta"),
            "pombe": ("pombe.bed", "pombe.fasta")
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
    if os.stat(file_name).st_size == 0:
        raise NameError('file: %s is empty' % file_name)


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

    for specie in fastq_files:
        # check if file exist
        if not os.path.isfile(out_prefix_dir + fastq_files[specie][0]):
            raise NameError('file: %s not found' % fastq_files[specie][0])
        yield get_r2_command(fastq_files[specie][0], fastq_files[specie][1])



# RemoveAdapters ----------------------------------------------------------

def remove_adapters(out_prefix_dir):
    """Removes adapters from fastq

    Removes the adpaters for each specie in the 5' position
    using the cutadapt software

    Args:
        out_prefix_dir: a direcoty indicating where the files are located.

    Returns:
        a list where each element is the command to run in the pipeline

     Raises:
        NameError: An error ocurrer when the input fastq files are not present
    """

    def remove_adapter_cmd(specie):
        check_file(out_prefix_dir + fastq_files[specie][0])
        remove_from_r1 = ['cutadapt -g',
                         adapters[specie][0],
                         out_prefix_dir + fastq_files[specie][0],
                         '-o', out_prefix_dir + trimmed_fastq[specie][0]
                         ]

        check_file(out_prefix_dir + fastq_files[specie][0])
        remove_from_r2 = [
                         'cutadapt -g',
                         adapters[specie][1],
                         out_prefix_dir + fastq_files[specie][1],
                         '-o', out_prefix_dir + trimmed_fastq[specie][1]
                         ]
        return (
               ' '.join(_ for _ in remove_from_r1),
               ' '.join(_ for _ in remove_from_r2)
               )

    commands_to_run = []
    for specie in trimmed_fastq:
        commands_to_run.append(remove_adapter_cmd(specie)[0])
        commands_to_run.append(remove_adapter_cmd(specie)[1])
 
    return commands_to_run


# collapse minigenes ------------------------------------------------------

def collapse_minigenes(out_prefix_dir):
    """Collapses the minigenes

    Collapses the trimmer minigenes with FASTA/Q Collapser

    Args:
        out_prefix_dir: a direcoty indicating where the files are located.

    Returns:
        the command to run the collapser

    Raises:
        NameError: An error ocurrer when the input fastq files are not present
    """

    [check_file(out_prefix_dir + trimmed_fastq["minigene"][i]) for i in (0, 1)]
    collapse_r1 = [
                  "fastx_collapser -v",
                  "-i", out_prefix_dir + trimmed_fastq["minigene"][0],
                  "-o", out_prefix_dir + "minigene_r1_collased.fq"
                  ]
    collapse_r2 = [
                  "fastx_collapser -v",
                  "-i", out_prefix_dir + trimmed_fastq["minigene"][1],
                  "-o", out_prefix_dir + "minigene_r2_collased.fq"
                  ]
    return (
           ' '.join(_ for _ in collapse_r1),
           ' '.join(_ for _ in collapse_r2)
           )

# QuantifyMinigenes -------------------------------------------------------

def quantify_minigenes(out_prefix_dir):
    """ Quantify the minigenes

    Computes the number of counts for each minigene using the fastqx barcode-sppliter,
    the file data/minigenes_barcodes.txt cotains the barcodes of the first 15 bp of the
    mingenes.

    Args:
        out_prefix_dir: a direcoty indicating where the file  minigene_01_trimmed.fq is located

    Returns:
        str, the command line to compute the counts

    Raises:
        NameError: An error ocurrer when the input file is not present
    """

    minigene_file = out_prefix_dir + "minigene_01_trimmed.fq"
    check_file(minigene_file)
    cuantify_minigenes = [
                         "cat", minigene_file,
                         "|",
                         "fastx_barcode_splitter.pl",
                         "--bcfile data/minigenes_barcodes.txt",
                         "--prefix", out_prefix_dir + "MiniGenes",
                         "--bol",
                         ">", out_prefix_dir + "minigene_quantifycation.txt"
                         ]

    # we only need the quantifications, so lets remove the files generated by fastx_barcode_splitter.pl

    rm_files = ['rm', out_prefix_dir + "MiniGenes*"]
    return (
           ' '.join(_ for _ in cuantify_minigenes),
           ' '.join(_ for _ in rm_files)
           )


# mapping -----------------------------------------------------------------

def mapping_bowtie(out_prefix_dir):
    """Maps reads to transcriptome

    Maps each of the species reads to its corresponding transcriptome

    Args:
        out_prefix_dir: a direcoty indicating where the trimmed fastq files are found

    Returns:
        iterator with the commands to run

    Raises:
        NameError: An error ocurrer when the input fastq files are not present
    """

    def bowtie_run(specie):
        r1 = out_prefix_dir + trimmed_fastq[specie][0]
        r2 = out_prefix_dir + trimmed_fastq[specie][1]
        reference = references[specie]
        check_file(r1)
        check_file(r2)
        out_bam = out_prefix_dir + specie + ".bam"
        bowtie_params = [
                        "bowtie",
                        "-n 2", # maximum two mistmaches
                        "--seedlen 10",
                        "-I 200",
                        "-X 800",
                        "--threads 10",
                        reference,
                        "-1", r1,
                        "-2", r2,
                        "-S", out_bam
                        ]
        return ' '.join(_ for _ in bowtie_params)

    for spe in ["fish", "fly", "pombe", "planaria"]:
        yield bowtie_run(spe)


def filter_mapped_reads(out_prefix_dir):
    """filter mapped reads from bam

    Filters the mapped reads with samtools flags

    Args:
        out_prefix_dir: a direcoty indicating where the bam files are found

    Returns:
        iterator with the commands to run
    """

    def filter_bam(specie):
        input_bam = out_prefix_dir + bam_files[specie][0]
        out_bam = out_prefix_dir + bam_files[specie][1]
        check_file(input_bam)
        filter_mapped = [
                        'samtools view',
                        '-f 0x02',
                        '-Sb',
                        input_bam,
                        '|',
                        'samtools sort',
                        '-m 6G -@ 4',
                        '-o', out_bam
                        ]
        return ' '.join(_ for _ in filter_mapped)

    for specie in bam_files:
        yield filter_bam(specie)



# ExtractSeqs -------------------------------------------------------------

# some error here have to change the variables
def to_bed(out_prefix_dir):
    """converts to bed formato

    filtered bam to bed
    """
    def bam_to_bed(specie):
        filter_bam = out_prefix_dir + bam_files[specie][1]
        check_file(filter_bam)
        bed_file = out_prefix_dir + fasta_bed[specie][0]

        cmds = [
               'bamToBed -bedpe',
               '-mate1',
               '-bed12',
               '-i', filter_bam,
               '|',
               'cut',
               '-f 1,2,6,7',
               '>', bed_file
               ]

        return ' '.join(_ for _ in cmds)

    for specie in bam_files:
       yield bam_to_bed(specie)


def extract_seqs(out_prefix_dir):
    """
    extract the mapped reads form the transcriptome
    """
    def get_seq(specie):
        transcriptome = references[specie] + ".fa"
        check_file(transcriptome)
        bed_file = out_prefix_dir + fasta_bed[specie][0]
        check_file(bed_file)
        out_fasta = out_prefix_dir + fasta_bed[specie][1]

        cmd = [
              'bedtools getfasta',
              '-fi', transcriptome,
              '-bed', bed_file,
              '-name',
              '>', out_fasta
              ]
        return ' '.join(_ for _ in cmd)

    for specie in bam_files:
        yield get_seq(specie)


