""""
arielome_pipe.py

Pipeline

"""

import os
import subprocess
import luigi
from helper import *

# global parameters configuration -----------------------------------------

class GlobalConfig(luigi.Config):
    """
    define the global parameters for the pipeline
    """
    outdir = luigi.Parameter(default = './')
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()
    prefix = luigi.Parameter(default = '')


# define pipe tasks -------------------------------------------------------

class BarcodeSpliting(luigi.Task):
    """
    split the fastqc libraries by barcode
    """

    def requires(self):
        return []

    def output(self):
        return luigi.LocalTarget(GlobalConfig().outdir + GlobalConfig().prefix + 'barcodesplit.log')

    def run(self):
        r1 = GlobalConfig().r1
        r2 = GlobalConfig().r2
        # check file exist
        check_file(r1)
        check_file(r2)

        barcode_r1 = barcode_split_cmd(r1, GlobalConfig().outdir + GlobalConfig().prefix, True)
        print("splitting library by barcode ...")
        #subprocess.call(barcode_r1, shell = True)

        with self.output().open('w') as barsplit:
            barsplit.write(barcode_r1)


if __name__ == '__main__':
    luigi.run()
