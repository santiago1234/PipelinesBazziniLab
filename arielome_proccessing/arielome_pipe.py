""""
arielome_pipe.py

Pipeline

"""

import os
import subprocess
import time
import luigi
import helper
# global parameters configuration -----------------------------------------

class GlobalConfig(luigi.Config):
    """
    define the global parameters for the pipeline
    """
    outdir = luigi.Parameter(default = './')
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()



# define pipe tasks -------------------------------------------------------

class Initializer(luigi.Task):
    """
    initialize the output dir and an info file
    """

    def requires(self):
        return[]

    def output(self):
        return luigi.LocalTarget(GlobalConfig().outdir + 'info.log')

    def run(self):
        with self.output().open('w') as info_file:
            info_file.write('pipeline info\nsample: %s\nproccessed on: %s\nautor: Santiago Medina'% (
                            GlobalConfig().outdir, time.strftime("%c")))


class BarcodeSpliting(luigi.Task):
    """
    split the fastqc libraries by barcode
    """

    def requires(self):
        return Initializer()

    def output(self):
        return luigi.LocalTarget(GlobalConfig().outdir + 'barcodesplit.log')

    def run(self):
        r1 = GlobalConfig().r1
        r2 = GlobalConfig().r2
        # check file exist
        helper.check_file(r1)
        helper.check_file(r2)

        barcode_r1 = helper.barcode_split_cmd(r1, GlobalConfig().outdir)
        print("splitting library by barcode ...")
        subprocess.call(barcode_r1, shell = True)

        # retrive read pair r2
        commands_to_run = open(GlobalConfig().outdir + "RunGetMaete.exe", 'w')
        for command_to_retrive in helper.get_r2_from_subset(r2, GlobalConfig().outdir):
            commands_to_run.write(command_to_retrive)
            commands_to_run.write("\n")
            print("runing: %s ..."  % command_to_retrive)
        commands_to_run.close()
        get_mate = GlobalConfig().outdir + "RunGetMaete.exe"
        get_mate = "cat" + " " + get_mate + " | " + "parallel --verbose"
        subprocess.call(get_mate, shell = True)

        with self.output().open('w') as barsplit:
            barsplit.write(barcode_r1 + '\n')
            barsplit.write(get_mate)


class RemoveAdapters(luigi.Task):
    """
    remove adapters task
    """

    def requires(self):
        return BarcodeSpliting()

    def output(self):
        return luigi.LocalTarget(GlobalConfig().outdir + 'barcode_trimming.log')

    def run(self):
        cutadpat_cmd = helper.remove_adapters(GlobalConfig().outdir)
        commands_to_run = open(GlobalConfig().outdir + "RemoveAdapters.exe", 'w')
        for cmd in cutadpat_cmd:
            print("runing: %s ..." % cmd)
            commands_to_run.write(cmd)
            commands_to_run.write("\n")
        commands_to_run.close()

        run_cutadapt = "cat" + " " + GlobalConfig().outdir + "RemoveAdapters.exe" + " | " + "parallel --verbose"
        subprocess.call(run_cutadapt, shell = True)


        with self.output().open('w') as cutadapt:
            cutadapt.write(run_cutadapt)


class CollapseMiniGenes(luigi.Task):
    """
    collapse the minigenes
    """

    def requires(self):
        return RemoveAdapters()

    def output(self):
        return luigi.LocalTarget(GlobalConfig().outdir + 'collapse_minigenes.log')

    def run(self):
        print('collapsing minigenes ...')
        collapse = helper.collapse_minigenes(GlobalConfig().outdir)
        subprocess.call(collapse[0], shell = True) #collapse r1
        subprocess.call(collapse[1], shell = True) #collapse r2

        with self.output().open('w') as fastx_collapse:
            fastx_collapse.write(collapse[0] + '\n' + collapse[1])


class QuantifyMinigenes(luigi.Task):
    """
    quantify the minigenes with fake minigene barcod,
    only reads R1 will be used
    """

    def requires(self):
        return RemoveAdapters()

    def output(self):
        return luigi.LocalTarget(GlobalConfig().outdir + "quantify_minigenes.log")

    def run(self):
        print('computing minigene counts ...')

        quantify = helper.quantify_minigenes(GlobalConfig().outdir)
        subprocess.call(quantify[0], shell = True)
        subprocess.call(quantify[1], shell = True)

        with self.output().open('w') as counts:
            counts.write(quantify[0])
            counts.write('\n')
            counts.write(quantify[1])
        counts.close()

        print("task completed!!!")


class Mapping(luigi.Task):
    """
    map the trimmed fastq to their corresponfing transcriptome with bowtie
    """

    def requires(self):
        return RemoveAdapters()

    def output(self):
        return luigi.LocalTarget(GlobalConfig().outdir + "bowtie_mapping.log")

    def run(self):
        print("mapping reads to transcriptome ...")
        bowtie_exe = open(GlobalConfig().outdir + 'bowtie.exe', 'w')
        bowtie_run = helper.mapping_bowtie(GlobalConfig().outdir)
        for cmd in bowtie_run:
            print("runing bowtie ...")
            subprocess.call(cmd, shell = True)
            bowtie_exe.write(cmd + '\n')

        with self.output().open("w") as bowtie:
            bowtie.write("bowtie was run succesfully")


class FilterAndSort(luigi.Task):
    """
    filters only reads mapped in proper pair &
    sorts the bam file
    """

    def requires(self):
        return Mapping()

    def output(self):
        return luigi.LocalTarget(GlobalConfig().outdir + "filter_sort.log")

    def run(self):
        print('filtering reads mapped in proper pair ...')
        filter_exe = open(GlobalConfig().outdir + 'filter.exe', 'w')
        for cmd in helper.filter_mapped_reads(GlobalConfig().outdir):
            print(cmd)
            filter_exe.write(cmd + '\n')
        filter_exe.close()
        filter_cmd = 'cat ' + GlobalConfig().outdir + 'filter.exe' + ' | ' + 'parallel --verbose'
        print('runing: %s ...' % filter_cmd)
        subprocess.call(filter_cmd, shell = True)
        with self.output().open('w') as filter_sort:
            filter_sort.write('filer and run')


class ExtractSeqs(luigi.Task):
    """
    extract the sequences of the mapped reads form the transcriptome
    """

    def requires(self):
        return [FilterAndSort(), QuantifyMinigenes()]

    def output(self):
        return luigi.LocalTarget(GlobalConfig().outdir + "extarct_seqs.log")

    def run(self):

        # convert bam to bed
        print("runing bam to bed ...")
        to_bed_exe = open(GlobalConfig().outdir + 'to_bed.exe', 'w')
        for cmd in helper.to_bed(GlobalConfig().outdir):
            print(cmd)
            to_bed_exe.write(cmd + '\n')
        to_bed_exe.close()

        to_bed_cmd = 'cat ' + GlobalConfig().outdir + 'to_bed.exe' + ' | ' + 'parallel --verbose'
        subprocess.call(to_bed_cmd, shell = True)

        # extract seqs
        print('extracting fasta sequences')
        extract = open(GlobalConfig().outdir + 'extract_seqs.exe', 'w')

        for cmd in helper.extract_seqs(GlobalConfig().outdir):
            extract.write(cmd + '\n')
        extract.close()
        extract_cmd = 'cat ' + GlobalConfig().outdir + 'extract_seqs.exe' + ' | ' + 'parallel --verbose'
        subprocess.call(extract_cmd, shell = True)

        # put sequences in correct fram and format in tabular
        print("putting sequences in frame ...")
        helper.format_in_frame(GlobalConfig().outdir)

        with self.output().open('w') as getseqs:
            getseqs.write('task completed!!')

if __name__ == '__main__':
    luigi.run()
