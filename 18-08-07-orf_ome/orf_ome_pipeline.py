import luigi
from os import path
import subprocess
import orfomeutlis as utils


class Params(luigi.Config):
    """
    input parameters for the pipeline
    """
    inputFq = luigi.Parameter()  # input fastq file
    adapter5 = luigi.Parameter()  # 5' adapter
    adapter3 = luigi.Parameter()  # 3' adapter


class MakeOutputDir(luigi.Task):
    """
    makes output dir to save results
    """

    def requires(self):
        return list()

    def output(self):
        return luigi.LocalTarget(utils.sample_name(Params().inputFq))

    def run(self):
        subprocess.call(f'mkdir {utils.sample_name(Params().inputFq)}',
                        shell=True)


class TrimAdapters(luigi.Task):
    """
    trims the 3' and 5' adapter
    """

    def requires(self):
        return MakeOutputDir()

    def output(self):
        return luigi.LocalTarget(
            path.join(utils.sample_name(Params().inputFq),
                      'trimed.fq')
        )

    def run(self):
        """
        run cutadapt to trim the adapters
        """
        cutadapt_comand = utils.run_cutadapt(
            Params().inputFq,
            Params().adapter5,
            Params().adapter3,
            utils.sample_name(Params().inputFq)
        )
        subprocess.call(cutadapt_comand, shell=True)


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
        return [MakeIndex(), TrimAdapters()]

    def output(self):
        return luigi.LocalTarget(path.join(
            utils.sample_name(Params().inputFq),
            'quants'
        ))

    def run(self):
        """
        runs salmon
        """
        cmd = f"""
        salmon quant -i data/barcodes_index
        -l A -r {path.join(utils.sample_name(Params().inputFq),'trimed.fq')}
        -o {path.join(utils.sample_name(Params().inputFq), 'quants')}
        """.replace('\n', ' ')
        subprocess.call(cmd, shell=True)



if __name__ == '__main__':
    luigi.run()
