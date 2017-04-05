import subprocess
import luigi
import os


# initialize parameters ---------------------------------------------------

OUT_PREFIX = "out/"
adapter_r1 = "ATGGTGAGCAAGGGCGAGGAGCTGTC"
adapter_r2 = "CTAGCTACCTA"
reference_transcriptome = "data/transcriptome/zfish"


# run command line --------------------------------------------------------

def run_cmd(cmd):
    p = subprocess.Popen(cmd, shell = False, universal_newlines = True,
                        stdout = subprocess.PIPE)
    ret_code = p.wait()
    output =  p.communicate()[0]
    return output

def run_cutadapt(adapter, fastq_file):
    """
    Runs the cuatdapt shell command
    """
    return run_cmd(["cutadapt", "-g", adapter, fastq_file])


##.......................... PIPELINE ...................................## 


class FastQC(luigi.Task):
    """
    runs fastqc quality check
    """
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget(OUT_PREFIX + 'fastqc.log')

    def run(self):
        shell_run = ' '.join(x for x in ['fastqc', self.r1, self.r2, '-o', OUT_PREFIX])
        subprocess.call(shell_run, shell = True)
        with self.output().open('w') as fastqc:
            fastqc.write('job completed!!')


class RemoveAdapters(luigi.Task):
    """
    remove 5' adapters from r1 and r2
    """

    r1 = luigi.Parameter()
    r2 = luigi.Parameter()

    def requires(self):
        yield FastQC(r1 = self.r1, r2 = self.r2)

    def output(self):
        return [luigi.LocalTarget(OUT_PREFIX + "r1.fq"),
                luigi.LocalTarget(OUT_PREFIX + "r2.fq")]

    def run(self):
        r1_trimmed = run_cutadapt(adapter_r1, self.r1)
        with self.output()[0].open('w') as cutadapt_out:
            cutadapt_out.write(r1_trimmed)

        r2_trimmed = run_cutadapt(adapter_r2, self.r2)
        with self.output()[1].open('w') as cutadapt_out:
            cutadapt_out.write(r2_trimmed)


class BuiltIndex(luigi.Task):
    """
    builts the transcriptome index
    """

    def requires(self):
        return []

    def output(self):
        return [luigi.LocalTarget(OUT_PREFIX + 'index.log')]

    def run(self):
        shell_run = 'bowtie2-build --threads 4 data/transcriptome/Zebrafish.fasta' + " " + reference_transcriptome
        print("runing on shell: %s" % shell_run)
        subprocess.call(shell_run, shell = True)
        with self.output()[0].open('w') as index:
            index.write('transcriptome index completed!!!')
        print("INDEX BUILT")


class RunAll(luigi.task.WrapperTask):
    """
    run the compleate pipeline
    """
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()

    def requires(self):
        yield RemoveAdapters(r1 = self.r1, r2 = self.r2)
        yield BuiltIndex()


# run pipeline ------------------------------------------------------------

if __name__ == '__main__':
    luigi.run()
