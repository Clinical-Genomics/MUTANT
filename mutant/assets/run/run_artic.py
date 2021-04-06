""" This script renames the output files for the ARTIC pipeline and Pangolin analysis,
    and creates a deliverables file for Clinical Genomics Infrastructure"""

import os
import sys
import click
import json
import subprocess
from mutant import version, log
from mutant.assets.utils.parse import get_json

class RunSC2:

    def __init__(self, input_folder, config_artic, caseID, prefix, profiles, timestamp, WD):

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
            general_config = self.get_json_data(config)
            resdir = os.path.join(general_config["SARS-CoV-2"]["folders"]["results"], "{}_{}".format(
                self.case, self.timestamp))
        else:
            resdir = "results"
        return resdir

    def run_case(self, resdir):

        """Run SARS-CoV-2 analysis"""

        resultsline = "--outdir {}".format(resdir)
        workline = "-work-dir {}".format(os.path.join(resdir, "work"))
        nflog = os.path.join(resdir, "nextflow.log")
        confline = ""
        if self.config_artic != "":
            confline = "-C {0}".format(self.config_artic)

        cmd = 'nextflow {0} -log {1} run {2} {3}/externals/gms-artic/main.nf -profile {4} --illumina --prefix {5} ' \
              '--directory {6} {7}'.format(confline, nflog, workline, self.WD, self.profiles, self.prefix, self.fastq, resultsline)
        log.debug("Command ran: {}".format(cmd))
        proc = subprocess.Popen(cmd.split())
        out, err = proc.communicate()
        log.info(out)
        log.info(err)

    def get_json_data(self, config):
        if os.path.exists(config):
            """Get sample information as json object"""
            try:
                with open(config) as json_file:
                    data = json.load(json_file)
            except Exception as e:
                click.echo("Unable to read provided json file: {}. Exiting..".format(config))
                click.echo(e)
                sys.exit(-1)
        else:
            click.echo("Could not find supplied config: {}. Exiting..".format(config))
            sys.exit(-1)
        return data


