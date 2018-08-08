import luigi
from os import path
import subprocess


class Params(luigi.Config):
    """
    input parameters for the pipeline
    Inputs:
        inputFq: path to the fastq file to process
        adapter5: 5' adapter
        adapter3: 3' adapter
        idprefix: prefix string for output file names
    """
    inputFq = luigi.Parameter()  # input fastq file
    adapter5 = luigi.Parameter()  # 5' adapter
    adapter3 = luigi.Parameter()  # 3' adapter
    idprefix = luigi.Parameter(default='p-')  # prefix for naming output files


class TrimAdapters(luigi.Task):
    """
    trims the 3' and 5' adapter
    The 5’ adapter is removed if it occurs. If a 3’ adapter occurs, it is
    removed only when also a 5’ adapter is present.
    """

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget(Params().idprefix + '_trimed.fq')

    def run(self):
        """
        run cutadapt to trim the adapters
        """
        cutadapt_cmd = f"""
        cutadapt -j 15 --discard-untrimmed
        -a {Params().adapter5}...{Params().adapter3}
        {Params().inputFq} > {Params().idprefix + '_trimed.fq'} 2>
        {Params().idprefix + '.report1'}
        """.replace('\n', ' ')
        subprocess.call(cutadapt_cmd, shell=True)


class RunQC(luigi.Task):
    """"
    runs quality analysis in input fastq file
    """

    def requires(self):
        return TrimAdapters()

    def output(self):
        return luigi.LocalTarget(Params().idprefix + '_fastqc.html')

    def run(self):
        qc_cmd = f"fastqc {Params().inputFq}"
        subprocess.call(qc_cmd, shell=True)


class MakeIndex(luigi.Task):
    """
    builts a salmon index for the barcodes
    NOTE: I set k = 5, since the barcodes are short.
    """

    def requires(self):
        return list()

    def output(self):
        return luigi.LocalTarget('data/barcodes_index')

    def run(self):
        built_index = 'salmon index -t data/barcodes.fasta -i data/barcodes_index --type quasi -k 5'
        subprocess.call(built_index, shell=True)


class QuantifyReads(luigi.Task):
    """
    quantify the trimmed reads with salmon
    """

    def requires(self):
        return [MakeIndex(), TrimAdapters(), RunQC()]

    def output(self):
        return luigi.LocalTarget(Params().idprefix + '-quants')

    def run(self):
        """
        runs salmon
        """
        cmd = f"""
        salmon quant -i data/barcodes_index
        -l A -r {Params().idprefix + '_trimed.fq'}
        -o {Params().idprefix + '_quants'}
        """.replace('\n', ' ')
        subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    luigi.run()
