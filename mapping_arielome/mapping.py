import subprocess
import luigi


# initialize parameters ---------------------------------------------------

OUT_PREFIX = "out/"
FILE_PREFIX = "test_" # a pefix to differentiate files
adapter_r1 = "ATGGTGAGCAAGGGCGAGGAGCTGTC"
adapter_r2 = "CTAGCTACCTA"
reference_transcriptome = "data/transcriptome/zfish"


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
        shell_run = 'bowtie-build  data/transcriptome/Zebrafish.fasta' + " " + reference_transcriptome
        print("runing on shell: %s" % shell_run)
        subprocess.call(shell_run, shell = True)
        with self.output()[0].open('w') as index:
            index.write('transcriptome index completed!!!')
        print("INDEX BUILT")


class BowtieMapping(luigi.Task):
    """"
    maps the trimmed reads to the transcriptome
    """
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()
    def requires(self):
        return [BuiltIndex(), RemoveAdapters(r1 = self.r1, r2 = self.r2)]

    def output(self):
        return luigi.LocalTarget(OUT_PREFIX + FILE_PREFIX +  'bowtie.log')

    def run(self):

        out_bam_file = OUT_PREFIX + FILE_PREFIX +  'aligment.sam'
        bowtie_params = ['bowtie',
                        '-n 2', # allow a maximum of 2 mistmaches
                        '--seedlen 10', # seed length
                        '-I 200',
                        '-X 800',
                        '--threads 20',
                        reference_transcriptome,
                        '-1', OUT_PREFIX + FILE_PREFIX +  'r1_trimmed.fq',
                        '-2', OUT_PREFIX + FILE_PREFIX +  'r2_trimmed.fq',
                        '-S', out_bam_file
                        ]

        # map with bowtie
 
        print('mapping ...')
        map_bowtie = ' '.join(_ for _ in bowtie_params)
        subprocess.call(map_bowtie, shell = True)

        # mapping statistics
        print('computing mapping stats ...')
        samtools_stats = ['samtools',
                         'flagstat',
                         out_bam_file,
                         '>',
                         OUT_PREFIX + FILE_PREFIX + "alg.stats"]
        samtools_stats = ' '.join(_ for _ in samtools_stats)
        subprocess.call(samtools_stats, shell = True)

        # filter mapped reads, sort
        print('filtering and sorting ...')
        filter_sort = ['samtools view',
                      '-f 0x02',
                      '-Sb', out_bam_file,
                      '|',
                      'samtools sort',
                      '-m 6G -@ 4',
                      '-o', OUT_PREFIX + FILE_PREFIX + 'alg_sorted.bam'
                      ]

        filter_sort = ' '.join(_ for _ in filter_sort)
        subprocess.call(filter_sort, shell = True)

        with self.output().open('w') as bowtie2:
            bowtie2.write('bowtie2 was run with the following params: \n %s\n %s\n %s' % (
                         map_bowtie, samtools_stats, filter_sort
                         ))




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
