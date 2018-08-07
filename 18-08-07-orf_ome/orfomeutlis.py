from os import path


def sample_name(fq_file):
    """
    gets the sample name from fastq file
    """
    return path.splitext(path.basename(fq_file))[0]


def run_cutadapt(fqfile, adapter5, adapter3, outdir):
    """
    generates comand to run cut adapt
    """
    cmd = f"""cutadapt
    -j 15
    --discard-untrimmed
    -a {adapter5}...{adapter3}
    {fqfile} > {path.join(outdir, 'trimed.fq')} 2>
    {path.join(outdir, 'report1.txt')}
    """.replace('\n', ' ')
    return cmd
