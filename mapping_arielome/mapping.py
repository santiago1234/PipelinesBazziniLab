import subprocess
import luigi


# initialize parameters ---------------------------------------------------

OUT_PREFIX = "out/"
FILE_PREFIX = "test_" # a pefix to differentiate files
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
        return luigi.LocalTarget(OUT_PREFIX + FILE_PREFIX + 'fastqc.log')

    def run(self):
        shell_run = ' '.join(x for x in ['fastqc', self.r1, self.r2, '-o', OUT_PREFIX])
        subprocess.call(shell_run, shell = True)
        with self.output().open('w') as fastqc:
            fastqc.write('fastqc was run with the following parameters:\n %s' % shell_run)


class RemoveAdapters(luigi.Task):
    """
    remove 5' adapters from r1 and r2
    """

    r1 = luigi.Parameter()
    r2 = luigi.Parameter()

    def requires(self):
        yield FastQC(r1 = self.r1, r2 = self.r2)

    def output(self):
        return luigi.LocalTarget(OUT_PREFIX + FILE_PREFIX +  'cutadapt.log')

    def run(self):
        shell_run_1 = ' '.join(['cutadapt', '-g', adapter_r1, self.r1, '-o', OUT_PREFIX +  FILE_PREFIX + 'r1_trimmed.fq'])
        shell_run_2 = ' '.join(['cutadapt', '-g', adapter_r2, self.r2, '-o', OUT_PREFIX + FILE_PREFIX +  'r2_trimmed.fq'])
        subprocess.call(shell_run_1, shell = True)
        subprocess.call(shell_run_2, shell = True)
        with self.output().open('w') as cutadapt:
            cutadapt.write('cutadapt was run with the following params: \n %s \n %s' % (shell_run_1, shell_run_2))


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


class Bowtie2Mapping(luigi.Task):
    """"
    maps the trimmed reads to the transcriptome
    """
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()
    def requires(self):
        return [BuiltIndex(), RemoveAdapters(r1 = self.r1, r2 = self.r2)]

    def output(self):
        return luigi.LocalTarget(OUT_PREFIX + FILE_PREFIX +  'bowtie2.log')

    def run(self):
        bowtie_params = ['bowtie2', '-k 1', '--local', '--no-mixed',
                        '--no-discordant', '--no-overlap', '--no-unal',
                        '--threads 4',
                        '-I 200', '-X 800', '-x data/transcriptome/zfish',
                        '-1', OUT_PREFIX + FILE_PREFIX +  'r1_trimmed.fq',
                        '-2', OUT_PREFIX + FILE_PREFIX +  'r2_trimmed.fq',
                        '-S', OUT_PREFIX + FILE_PREFIX +  'aligment.sam',
                        '--rg-id ' + FILE_PREFIX]

        # map with bowtie2
        print('mapping ...')
        map_bowtie = ' '.join(_ for _ in bowtie_params)
        subprocess.call(map_bowtie, shell = True)

        # convert also to bam format
        print('converting sam to bam ...')
        bam_out_file = OUT_PREFIX + FILE_PREFIX + "aligment.bam"
        shell_sam_to_bam = "samtools view -b " + OUT_PREFIX + FILE_PREFIX + 'aligment.sam' + " > " + bam_out_file
        subprocess.call(shell_sam_to_bam, shell = True)

        # sort bam
        print('sorting bam file ...')
        sort_bam = 'samtools sort -m 5G -@ 4 ' + bam_out_file + ' > ' + bam_out_file
        subprocess.call(sort_bam, shell = True)

        # index
        print('indexing bam file ...')
        index_bam = 'samtools index ' + bam_out_file 
        subprocess.call(index_bam, shell = True)

        with self.output().open('w') as bowtie2:
            bowtie2.write('bowtie2 was run with the following params: \n %s\n %s\n %s\n %s' % (map_bowtie, shell_sam_to_bam, sort_bam, index_bam))


class RunAll(luigi.task.WrapperTask):
    """
    run the compleate pipeline
    """
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()

    def requires(self):
        yield Bowtie2Mapping(r1 = self.r1, r2 = self.r2)

# run pipeline ------------------------------------------------------------

if __name__ == '__main__':
    luigi.run()
