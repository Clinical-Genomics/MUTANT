""" This script renames the output files for the ARTIC pipeline and Pangolin analysis,
    and creates a deliverables file for Clinical Genomics Infrastructure"""

import os
import subprocess
from mutant import log
from mutant.modules.generic_parser import get_json


class RunSC2:
    def __init__(self, input_folder, caseID, prefix, profiles, timestamp, WD, config_artic=""):
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
                    f"{self.case}_{self.timestamp}",
                )
            )
        else:
            resdir = os.path.abspath("results")
        return resdir

    def run_case(self, resdir, nanopore):
        """Run SARS-CoV-2 analysis"""

        resultsline = f"--outdir {resdir}"
        workline = f"-work-dir {os.path.join(resdir, 'work')}"
        nflog = os.path.join(resdir, "nextflow.log")
        confline = ""
        if self.config_artic != "":
            confline = f"-C {self.config_artic}"

        if nanopore:
            cmd = (
                f"nextflow {confline} -log {nflog} run {workline} "
                f"{self.WD}/externals/gms-artic/main.nf -profile {self.profiles} --medaka "
                f"--prefix {self.prefix} --basecalled_fastq {self.fastq} {resultsline}"
            )
        else:
            cmd = (
                f"nextflow {confline} -log {nflog} run {workline} "
                f"{self.WD}/externals/gms-artic/main.nf -profile {self.profiles} --illumina "
                f"--prefix {self.prefix} --directory {self.fastq} {resultsline}"
            )

        log.debug(f"Command ran: {cmd}")
        proc = subprocess.Popen(cmd.split())
        out, err = proc.communicate()
        log.info(out)
        log.info(err)
