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
        get_mate = "cat" + " " + get_mate + " | " + "parallel"
        subprocess.call(get_mate, shell = True)

        with self.output().open('w') as barsplit:
            barsplit.write(barcode_r1 + '\n')
            barsplit.write(get_mate)

class RemoveAdapters(luigi.Task):
    pass


if __name__ == '__main__':
    luigi.run()
