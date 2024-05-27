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

        concat_consensus = f"{self.indir}/{self.ticket}.consensus.fa"
        concat = open(concat_consensus, "w+")
        for item in glob.glob(target_files):
            single = open(item, "r")
            concat.write(single.read())
            concat.write("\n")
        concat.close()

    def create_delivery_metrics(self, fastq_dir: str) -> None:
        """Create deliverables file"""

        deliverables = {"files": []}
        deliverables_file = f"{self.indir}/{self.case}_deliverables.yaml"

        # Instrument properties
        deliverables["files"].append(
            {
                "format": "txt",
                "id": self.case,
                "path": f"{self.indir}/instrument.properties",
                "path_index": "~",
                "step": "report",
                "tag": "instrument-properties",
            }
        )
        # KS Report
        deliverables["files"].append(
            {
                "format": "csv",
                "id": self.case,
                "path": f"{self.indir}/sars-cov-2_{self.ticket}_results.csv",
                "path_index": "~",
                "step": "report",
                "tag": "ks-results",
            }
        )
        # KS Aux report
        deliverables["files"].append(
            {
                "format": "csv",
                "id": self.case,
                "path": f"{self.indir}/sars-cov-2_{self.ticket}_variants.csv",
                "path_index": "~",
                "step": "report",
                "tag": "ks-aux-results",
            }
        )
        # Pangolin typing
        deliverables["files"].append(
            {
                "format": "csv",
                "id": self.case,
                "path": f"{self.indir}/{self.ticket}.pangolin.csv",
                "path_index": "~",
                "step": "report",
                "tag": "pangolin-typing",
            }
        )
        # Pangolin typing for FOHM (only qcpass files)
        deliverables["files"].append(
            {
                "format": "csv",
                "id": self.case,
                "path": f"{self.indir}/{self.ticket}_{date.today()}_pangolin_classification_format4.txt",
                "path_index": "~",
                "step": "report",
                "tag": "pangolin-typing-fohm",
            }
        )
        # Consensus file
        deliverables["files"].append(
            {
                "format": "fasta",
                "id": self.case,
                "path": f"{self.indir}/{self.ticket}.consensus.fa",
                "path_index": "~",
                "step": "analysis",
                "tag": "consensus",
            }
        )
        # Multiqc report
        deliverables["files"].append(
            {
                "format": "html",
                "id": self.case,
                "path": f"{self.indir}/{self.ticket}_multiqc.html",
                "path_index": "~",
                "step": "report",
                "tag": "multiqc-html",
            }
        )
        # MultiQC json
        deliverables["files"].append(
            {
                "format": "json",
                "id": self.case,
                "path": f"{self.indir}/{self.ticket}_multiqc.json",
                "path_index": "~",
                "step": "report",
                "tag": "multiqc-json",
            }
        )
        deliverables["files"].append(
            {
                "format": "csv",
                "id": self.case,
                "path": f"{self.indir}/nextclade_summary.csv",
                "path_index": "~",
                "step": "report",
                "tag": "nextclade-summary",
            }
        )
        if not self.nanopore:
            # Artic yaml (Vogue) data
            deliverables["files"].append(
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
        deliverables["files"].append(
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
        deliverables["files"].append(
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
        deliverables["files"].append(
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
        deliverables["files"].append(
            {
                "format": "txt",
                "id": self.case,
                "path": f"{self.indir}/nextflow.log",
                "path_index": "~",
                "step": "runinfo",
                "tag": "logfile",
            }
        )
        # FoHM delivery file
        deliverables["files"].append(
            {
                "format": "csv",
                "id": self.case,
                "path": os.path.join(
                    self.indir, f"{self.ticket}_komplettering.csv"
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
            base_sample = f"{region}_{lab}_{sample}"
            if self.nanopore:
                # Concat reads forwards
                deliverables["files"].append(
                    {
                        "format": "fastq",
                        "id": sampleID,
                        "path": f"{fastq_dir}/{sample}/{base_sample}_1.fastq.gz",
                        "path_index": "~",
                        "step": "concatination",
                        "tag": "forward-reads",
                    }
                )
                # Variants (vcf)
                deliverables["files"].append(
                    {
                        "format": "vcf",
                        "id": sampleID,
                        "path": f"{self.indir}/articNcovNanopore_Genotyping_typeVariants/vcf/{base_sample}.vcf",
                        "path_index": "~",
                        "step": "genotyping",
                        "tag": "vcf-covid",
                    }
                )
                # Single-file fasta
                deliverables["files"].append(
                    {
                        "format": "fasta",
                        "id": sampleID,
                        "path": f"{self.indir}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka/{base_sample}.consensus.fasta",
                        "path_index": "~",
                        "step": "consensus",
                        "tag": "consensus-sample",
                    }
                )
            else:
                if not record["sequencing_qc_pass"]:
                    continue
                # Concat reads forwards
                deliverables["files"].append(
                    {
                        "format": "fastq",
                        "id": sampleID,
                        "path": f"{fastq_dir}/{base_sample}_1.fastq.gz",
                        "path_index": "~",
                        "step": "concatination",
                        "tag": "forward-reads",
                    }
                )
                # Concat reads reverse
                deliverables["files"].append(
                    {
                        "format": "fastq",
                        "id": sampleID,
                        "path": f"{fastq_dir}/{base_sample}_2.fastq.gz",
                        "path_index": "~",
                        "step": "concatination",
                        "tag": "reverse-reads",
                    }
                )
                # Variants (vcf)
                deliverables["files"].append(
                    {
                        "format": "vcf",
                        "id": sampleID,
                        "path": f"{self.indir}/ncovIllumina_Genotyping_typeVariants/vcf/{base_sample}.vcf",
                        "path_index": "~",
                        "step": "genotyping",
                        "tag": "vcf-covid",
                    }
                )
                # Single-file fasta
                deliverables["files"].append(
                    {
                        "format": "fasta",
                        "id": sampleID,
                        "path": (
                            f"{self.indir}/ncovIllumina_sequenceAnalysis_makeConsensus/{base_sample}.consensus.fasta"
                        ),
                        "path_index": "~",
                        "step": "consensus",
                        "tag": "consensus-sample",
                    }
                )
        with open(deliverables_file, "w") as out:
            yaml.dump(deliverables, out)
