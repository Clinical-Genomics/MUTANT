import glob
import os
from datetime import date
from pathlib import Path

import yaml

from mutant.modules.generic_parser import get_sarscov2_config, get_results_paths


class GenericReporter:
    def __init__(
        self,
        caseinfo: dict,
        casefile: str,
        indir: str,
        nanopore: bool,
        config_artic: str,
    ):
        self.casefile = casefile
        self.caseinfo = get_sarscov2_config(self.casefile)
        self.indir = indir
        self.config_artic = config_artic
        self.nanopore = nanopore
        self.ticket = caseinfo[0]["Customer_ID_project"]
        self.case = caseinfo[0]["case_ID"]
        self.filepaths = get_results_paths(self.indir, self.case, self.ticket, nanopore)

    def get_finished_slurm_ids(self) -> list:
        """Get slurm IDs"""

        trace_file_path = Path(self.indir, "pipeline_info", "execution_trace.txt")
        slurm_id_list = []
        with open(trace_file_path, "r") as trace_file_contents:
            for line in trace_file_contents:
                slurm_id = line.split()[2]
                try:
                    slurm_id_list.append(int(slurm_id))
                except Exception:
                    continue
        return slurm_id_list

    def create_trailblazer_config(self) -> None:
        """Create Trailblazer config file"""

        trailblazer_config_path = Path(self.indir, "trailblazer_config.yaml")
        finished_slurm_ids = self.get_finished_slurm_ids()
        if not finished_slurm_ids:
            return
        with open(trailblazer_config_path, "w") as trailblazer_config_file:
            yaml.dump(data={"jobs": finished_slurm_ids}, stream=trailblazer_config_file)

    def create_concat_consensus(self, target_files: str) -> None:
        """Concatenate consensus files"""

        concat_consensus = "{0}/{1}.consensus.fa".format(self.indir, self.ticket)
        concat = open(concat_consensus, "w+")
        for item in glob.glob(target_files):
            single = open(item, "r")
            concat.write(single.read())
            concat.write("\n")
        concat.close()

    def create_deliveryfile(self, fastq_dir: str) -> None:
        """Create deliverables file"""

        deliv = {"files": []}
        delivfile = "{}/{}_deliverables.yaml".format(self.indir, self.case)

        ## Per Case
        # Instrument properties
        deliv["files"].append(
            {
                "format": "txt",
                "id": self.case,
                "path": "{}/instrument.properties".format(self.indir),
                "path_index": "~",
                "step": "report",
                "tag": "instrument-properties",
            }
        )
        # KS Report
        deliv["files"].append(
            {
                "format": "csv",
                "id": self.case,
                "path": "{}/sars-cov-2_{}_results.csv".format(self.indir, self.ticket),
                "path_index": "~",
                "step": "report",
                "tag": "ks-results",
            }
        )
        # KS Aux report
        deliv["files"].append(
            {
                "format": "csv",
                "id": self.case,
                "path": "{}/sars-cov-2_{}_variants.csv".format(self.indir, self.ticket),
                "path_index": "~",
                "step": "report",
                "tag": "ks-aux-results",
            }
        )
        # Pangolin typing
        deliv["files"].append(
            {
                "format": "csv",
                "id": self.case,
                "path": "{}/{}.pangolin.csv".format(self.indir, self.ticket),
                "path_index": "~",
                "step": "report",
                "tag": "pangolin-typing",
            }
        )
        # Pangolin typing for FOHM (only qcpass files)
        deliv["files"].append(
            {
                "format": "csv",
                "id": self.case,
                "path": "{}/{}_{}_pangolin_classification_format4.txt".format(
                    self.indir, self.ticket, str(date.today())
                ),
                "path_index": "~",
                "step": "report",
                "tag": "pangolin-typing-fohm",
            }
        )
        # Consensus file
        deliv["files"].append(
            {
                "format": "fasta",
                "id": self.case,
                "path": "{}/{}.consensus.fa".format(self.indir, self.ticket),
                "path_index": "~",
                "step": "analysis",
                "tag": "consensus",
            }
        )
        if not self.nanopore:
            # Multiqc report
            deliv["files"].append(
                {
                    "format": "html",
                    "id": self.case,
                    "path": "{}/{}_multiqc.html".format(self.indir, self.ticket),
                    "path_index": "~",
                    "step": "report",
                    "tag": "multiqc-html",
                }
            )
            # MultiQC json
            deliv["files"].append(
                {
                    "format": "json",
                    "id": self.case,
                    "path": "{}/{}_multiqc.json".format(self.indir, self.ticket),
                    "path_index": "~",
                    "step": "report",
                    "tag": "multiqc-json",
                }
            )
            # Artic yaml (Vogue) data
            deliv["files"].append(
                {
                    "format": "yaml",
                    "id": self.case,
                    "path": self.filepaths[self.case]["vogue-metrics"],
                    "path_index": "~",
                    "step": "result_aggregation",
                    "tag": "metrics",
                }
            )
        # Provided CG CASE info from StatusDB
        deliv["files"].append(
            {
                "format": "json",
                "id": self.case,
                "path": os.path.abspath(self.casefile),
                "path_index": "~",
                "step": "runinfo",
                "tag": "sampleinfo",
            }
        )
        # Input settings dump
        deliv["files"].append(
            {
                "format": "txt",
                "id": self.case,
                "path": os.path.abspath(self.config_artic),
                "path_index": "~",
                "step": "runinfo",
                "tag": "runtime-settings",
            }
        )
        # Software versions
        deliv["files"].append(
            {
                "format": "csv",
                "id": self.case,
                "path": self.filepaths[self.case]["versions-file"],
                "path_index": "~",
                "step": "runinfo",
                "tag": "software-versions",
            }
        )
        # Execution log
        deliv["files"].append(
            {
                "format": "txt",
                "id": self.case,
                "path": "{}/nextflow.log".format(self.indir),
                "path_index": "~",
                "step": "runinfo",
                "tag": "logfile",
            }
        )
        # FoHM delivery file
        deliv["files"].append(
            {
                "format": "csv",
                "id": self.case,
                "path": os.path.join(
                    self.indir, "{}_komplettering.csv".format(self.ticket)
                ),
                "path_index": "~",
                "step": "report",
                "tag": "SARS-CoV-2-info",
            }
        )

        # Per sample
        for record in self.caseinfo:
            sampleID = record["CG_ID_sample"]
            sample = record["Customer_ID_sample"]
            region = record["region_code"]
            lab = record["lab_code"]
            base_sample = "{0}_{1}_{2}".format(region, lab, sample)
            if self.nanopore:
                # Concat reads forwards
                deliv["files"].append(
                    {
                        "format": "fastq",
                        "id": sampleID,
                        "path": "{0}/{1}/{2}_1.fastq.gz".format(
                            fastq_dir, sample, base_sample
                        ),
                        "path_index": "~",
                        "step": "concatination",
                        "tag": "forward-reads",
                    }
                )
                # Variants (vcf)
                deliv["files"].append(
                    {
                        "format": "vcf",
                        "id": sampleID,
                        "path": "{}/articNcovNanopore_Genotyping_typeVariants/vcf/{}.vcf".format(
                            self.indir, base_sample
                        ),
                        "path_index": "~",
                        "step": "genotyping",
                        "tag": "vcf-covid",
                    }
                )
                # Single-file fasta
                deliv["files"].append(
                    {
                        "format": "fasta",
                        "id": sampleID,
                        "path": "{}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka/{}.consensus.fasta".format(
                            self.indir, base_sample
                        ),
                        "path_index": "~",
                        "step": "consensus",
                        "tag": "consensus-sample",
                    }
                )
            else:
                if not record["sequencing_qc_pass"]:
                    continue
                # Concat reads forwards
                deliv["files"].append(
                    {
                        "format": "fastq",
                        "id": sampleID,
                        "path": "{0}/{1}_1.fastq.gz".format(fastq_dir, base_sample),
                        "path_index": "~",
                        "step": "concatination",
                        "tag": "forward-reads",
                    }
                )
                # Concat reads reverse
                deliv["files"].append(
                    {
                        "format": "fastq",
                        "id": sampleID,
                        "path": "{0}/{1}_2.fastq.gz".format(fastq_dir, base_sample),
                        "path_index": "~",
                        "step": "concatination",
                        "tag": "reverse-reads",
                    }
                )
                # Variants (vcf)
                deliv["files"].append(
                    {
                        "format": "vcf",
                        "id": sampleID,
                        "path": "{}/ncovIllumina_Genotyping_typeVariants/vcf/{}.vcf".format(
                            self.indir, base_sample
                        ),
                        "path_index": "~",
                        "step": "genotyping",
                        "tag": "vcf-covid",
                    }
                )
                # Single-file fasta
                deliv["files"].append(
                    {
                        "format": "fasta",
                        "id": sampleID,
                        "path": "{}/ncovIllumina_sequenceAnalysis_makeConsensus/{}.consensus.fasta".format(
                            self.indir, base_sample
                        ),
                        "path_index": "~",
                        "step": "consensus",
                        "tag": "consensus-sample",
                    }
                )
        with open(delivfile, "w") as out:
            yaml.dump(deliv, out)
