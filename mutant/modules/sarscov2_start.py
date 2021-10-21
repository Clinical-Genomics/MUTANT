""" This script renames the output files for the ARTIC pipeline and Pangolin analysis,
    and creates a deliverables file for Clinical Genomics Infrastructure"""

import os
import sys
import click
import json
import subprocess
from mutant import version, log
from mutant.modules.generic_parser import get_json


class RunSC2:
    def __init__(
        self, input_folder, caseID, prefix, profiles, timestamp, WD, config_artic=""
    ):

        self.fastq = input_folder
        self.timestamp = timestamp
        self.WD = WD
        self.case = caseID
        self.prefix = prefix
        self.config_artic = config_artic
        self.profiles = profiles

    def get_results_dir(self, config, outdir):

        """Return result output directory"""

        if outdir != "":
            resdir = outdir
        elif config != "":
            general_config = get_json(config)
            resdir = os.path.abspath(
                os.path.join(
                    general_config["SARS-CoV-2"]["folders"]["results"],
                    "{}_{}".format(self.case, self.timestamp),
                )
            )
        else:
            resdir = os.path.abspath("results")
        return resdir

    def run_case(self, resdir, nanopore):

        """Run SARS-CoV-2 analysis"""

        resultsline = "--outdir {}".format(resdir)
        workline = "-work-dir {}".format(os.path.join(resdir, "work"))
        nflog = os.path.join(resdir, "nextflow.log")
        confline = ""
        if self.config_artic != "":
            confline = "-C {0}".format(self.config_artic)

        if nanopore:
            cmd = "nextflow {0} -log {1} run {2} {3}/externals/gms-artic/main.nf -profile singularity --nanopolish --prefix 211021_nanopore --basecalled_fastq /home/hiseq.clinical/HO_data_processing/nanopore/210811_47CoV_SABasecalled/CS5/20210811_1157_MC-111732_0_FAQ57606_c89872a3/fastq_pass -fast5_pass /home/hiseq.clinical/HO_data_processing/nanopore/210811_47CoV_SABasecalled/CS5/20210811_1157_MC-111732_0_FAQ57606_c89872a3/fast5_pass --sequencing_summary /home/hiseq.clinical/HO_data_processing/nanopore/210811_47CoV_SABasecalled/CS5/20210811_1157_MC-111732_0_FAQ57606_c89872a3/sequencing_summary_FAQ57606_71c83ae0.txt".format(
                confline,
                nflog,
                workline,
                self.WD,
            )
        else:
            cmd = "nextflow {0} -log {1} run {2} {3}/externals/gms-artic/main.nf -profile {4} --illumina --prefix {5} --directory {6} {7}".format(
                confline,
                nflog,
                workline,
                self.WD,
                self.profiles,
                self.prefix,
                self.fastq,
                resultsline,
            )
        log.debug("Command ran: {}".format(cmd))
        proc = subprocess.Popen(cmd.split())
        out, err = proc.communicate()
        log.info(out)
        log.info(err)
