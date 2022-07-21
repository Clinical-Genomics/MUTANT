""" This class creates reports. Specifically it acts on the sarscov2 pipeline,
    and creates report files for Clinical Genomics Infrastructure

    By: Isak Sylvin & Tanja Normark
"""
import os
import sys
import csv
import glob
import re
import yaml
import json
from datetime import date
from pathlib import Path

from mutant.constants.artic import ILLUMINA_FILES_CASE, NANOPORE_FILES_CASE
from mutant.modules.generic_parser import (
    get_sarscov2_config,
    read_filelines,
    get_artic_results,
)
from mutant.modules.artic_illumina.parser import get_vogue_multiqc_data


def get_results_paths(indir, case, ticket, nanopore) -> dict:
    """Get paths for all reports and output files from GMS-Artic and MUTANT"""

    # Case paths
    path_dict = dict()
    case_paths = dict()
    if nanopore:
        casefiles = NANOPORE_FILES_CASE
    else:
        casefiles = ILLUMINA_FILES_CASE
    for stepname, filepath in casefiles.items():
        fullpath = filepath.format(resdir=indir, case=case, ticket=ticket)
        if "*" in fullpath:
            fullpath = glob.glob(fullpath)[0]
        case_paths[stepname] = fullpath
    path_dict[case] = case_paths
    return path_dict


class ReportSC2:
    def __init__(self, caseinfo, indir, config_artic, fastq_dir, timestamp, nanopore):
        self.casefile = caseinfo
        caseinfo = get_sarscov2_config(caseinfo)
        self.caseinfo = caseinfo
        self.case = caseinfo[0]["case_ID"]
        self.ticket = caseinfo[0]["Customer_ID_project"]
        self.project = caseinfo[0]["Customer_ID_project"]
        self.indir = indir
        self.config_artic = config_artic
        self.time = timestamp
        today = date.today().strftime("%Y%m%d")
        self.today = today
        self.fastq_dir = fastq_dir
        self.filepaths = get_results_paths(self.indir, self.case, self.ticket, nanopore)
        self.articdata = dict()
        self.nanopore = nanopore

    def create_all_files(self, nanopore: bool):
        self.create_trailblazer_config()
        self.load_lookup_dict()
        self.create_concat_pangolin(nanopore=nanopore)
        self.create_concat_pangolin_fohm(nanopore=nanopore)
        # This works off concat pango and needs to occur after
        self.create_concat_consensus(nanopore=nanopore)
        self.create_deliveryfile(nanopore=nanopore)
        if not nanopore:
            self.create_vogue_metrics_file()
        self.create_fohm_csv()
        self.create_jsonfile()
        self.create_instrument_properties(nanopore=nanopore)

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

    def load_lookup_dict(self):
        """Loads articdata with data from various sources. Atm, artic output and the case config input file"""

        self.load_artic_results()
        self.load_case_config()

    def load_case_config(self):
        """Appends additional data to articdata dictionary"""

        casekeys = self.caseinfo[0].keys()
        packing = dict(zip(casekeys, "-" * len(casekeys)))
        # Updates existing samples with defaults for case-config
        for k, v in self.articdata.items():
            self.articdata[k].update(packing)
        # Updates existing samples with provided case-config info
        for entry in self.caseinfo:
            k = entry["Customer_ID_sample"]
            if k in self.articdata.keys():
                self.articdata[k].update(entry)
            else:
                self.articdata[k] = entry

    def load_artic_results(self, nanopore: bool):
        """Parse artic output directory for analysis results. Returns dictionary data object"""
        artic_data = get_artic_results(self.indir, nanopore)
        self.articdata.update(artic_data)

    def create_vogue_metrics_file(self):
        """Create file with metrics for vogue"""

        out_yaml = {"metrics": []}
        # Get multiqc-data
        out_yaml["metrics"].extend(
            get_vogue_multiqc_data(
                self.filepaths[self.case]["multiqc-json"], self.caseinfo
            )
        )
        # Add info about failed/passed QC
        for record in self.caseinfo:
            if record["sequencing_qc_pass"]:
                sample_yaml = {
                    "header": "~",
                    "id": record["CG_ID_sample"],
                    "input": os.path.basename(
                        self.filepaths[self.case]["results-file"]
                    ),
                    "name": "passed-qc",
                    "step": "QC-threshold",
                    "value": self.articdata[record["Customer_ID_sample"]]["qc"],
                }
            # Sample failed QC and has not been analyzed
            else:
                sample_yaml = {
                    "header": "~",
                    "id": record["CG_ID_sample"],
                    "input": os.path.basename(
                        self.filepaths[self.case]["results-file"]
                    ),
                    "name": "passed-qc",
                    "step": "QC-threshold",
                    "value": "FALSE",
                }
            out_yaml["metrics"].append(sample_yaml)
        # Create output file
        with open(self.filepaths[self.case]["vogue-metrics"], "w") as out:
            yaml.dump(out_yaml, out)

    def create_concat_pangolin(self, nanopore: bool):
        """Concatenate pangolin results"""
        if nanopore:
            indir = (
                "{0}/articNcovNanopore_sequenceAnalysisMedaka_pangolinTyping".format(
                    self.indir
                )
            )
        else:
            indir = "{0}/ncovIllumina_sequenceAnalysis_pangolinTyping".format(
                self.indir
            )
        concatfile = "{0}/{1}.pangolin.csv".format(self.indir, self.ticket)
        pangolins = glob.glob("{0}/*.pangolin.csv".format(indir))
        # Copy header
        header = read_filelines(pangolins[0])[0]
        with open(concatfile, "w") as concat:
            concat.write(header)
            # Parse sample pangolin data
            for pango in pangolins:
                data = read_filelines(pango)[1:]
                for line in data:
                    # Use sample name at taxon field
                    taxon_regex = "(\w+)_(\w+)_(\w+)_(?P<name>\w+).(\S+)"
                    sample, subs = re.subn(taxon_regex, r"\g<name>", line.split(",")[0])
                    if subs == 0:
                        print(
                            "Unable to rename taxon - using original: {}".format(pango)
                        )
                    else:
                        line = sample + line[line.find(",") :]
                    concat.write(line)

    def create_concat_pangolin_fohm(self, nanopore: bool):
        """Concatenate pangolin results and format for fohm"""
        if nanopore:
            indir = (
                "{0}/articNcovNanopore_sequenceAnalysisMedaka_pangolinTyping".format(
                    self.indir
                )
            )
        else:
            indir = "{0}/ncovIllumina_sequenceAnalysis_pangolinTyping".format(
                self.indir
            )
        concatfile = "{0}/{1}_{2}_pangolin_classification_format4.txt".format(
            self.indir, self.ticket, str(date.today())
        )
        pangolins = glob.glob("{0}/*.pangolin.csv".format(indir))
        # Copy header
        header = read_filelines(pangolins[0])[0]
        with open(concatfile, "w") as concat:
            concat.write(header)
            # Parse sample pangolin data
            for pango in pangolins:
                with open(pango, "r") as pangolinfile:
                    pango_data = pangolinfile.readlines()[1]
                    csv_items = pango_data.split(",")
                    for sample, data in self.articdata.items():
                        if sample in csv_items[0]:
                            if data.get("qc") != "TRUE":
                                continue
                            csv_items[0] = sample
                            concat.write(",".join(csv_items))

    def create_concat_consensus(self, nanopore: bool):
        """Concatenate consensus files"""

        if nanopore:
            indir = (
                "{0}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka".format(
                    self.indir
                )
            )
        else:
            indir = "{0}/ncovIllumina_sequenceAnalysis_makeConsensus".format(self.indir)
        concat = open("{0}/{1}.consensus.fa".format(self.indir, self.ticket), "w+")
        if nanopore:
            consensus_postfix = "{0}/*.consensus.fasta".format(indir)
        else:
            consensus_postfix = "{0}/*.consensus.fa".format(indir)
        for item in glob.glob(consensus_postfix):
            single = open(item, "r")
            concat.write(single.read())
            concat.write("\n")
        concat.close()

    def create_fohm_csv(self):
        """Creates a summary file for FoHMs additional info"""

        sumfile = os.path.join(
            self.indir,
            "{}_komplettering.csv".format(self.ticket),
        )
        with open(sumfile, "w") as out:
            summary = csv.writer(out)
            # Add header
            summary.writerow(["provnummer", "urvalskriterium", "GISAID_accession"])
            # Write sample information
            for sample, data in self.articdata.items():
                if data.get("qc") != "TRUE":
                    continue
                summary.writerow(
                    [
                        sample,
                        data.get("selection_criteria", "Information saknas"),
                    ]
                )

    def create_instrument_properties(self, nanopore: bool):
        """Creates a properties file with instrument information"""

        propfile = os.path.join(self.indir, "instrument.properties")

        if nanopore:
            plfrm = "nanopore"
        else:
            plfrm = "illumina"
        lanes = "1"
        ms = "N/A"
        ml = "N/A"

        for sample, data in self.articdata.items():
            if ml == "N/A":
                ml = data["method_libprep"]
            elif ml != data["method_libprep"]:
                ml = "INCONSISTENT"
            if ms == "N/A":
                ms = data["method_sequencing"]
            elif ms != data["method_sequencing"]:
                ms = "INCONSISTENT"

        with open(propfile, "w") as prop:
            prop.write("instrument={}\n".format(ms))
            prop.write("plattform={}\n".format(plfrm))
            prop.write("biblioteksmetod={}\n".format(ml))
            prop.write("lanes={}\n".format(lanes))

    def create_jsonfile(self):
        """Output all result data in a json format for easy parsing"""

        if self.articdata == dict():
            print("No artic results loaded. Quitting create_jsonfile")
            sys.exit(-1)

        with open(
            "{}/{}_artic.json".format(self.indir, self.ticket, self.today), "w"
        ) as outfile:
            json.dump(self.articdata, outfile)

    def create_deliveryfile(self, nanopore: bool):
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
        if not nanopore:
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
        if not nanopore:
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
        if not nanopore:
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
            if not record["sequencing_qc_pass"]:
                continue
            # Concat reads forwards
            deliv["files"].append(
                {
                    "format": "fastq",
                    "id": sampleID,
                    "path": "{0}/{1}_1.fastq.gz".format(self.fastq_dir, base_sample),
                    "path_index": "~",
                    "step": "concatination",
                    "tag": "forward-reads",
                }
            )
            if not nanopore:
                # Concat reads reverse
                deliv["files"].append(
                    {
                        "format": "fastq",
                        "id": sampleID,
                        "path": "{0}/{1}_2.fastq.gz".format(
                            self.fastq_dir, base_sample
                        ),
                        "path_index": "~",
                        "step": "concatination",
                        "tag": "reverse-reads",
                    }
                )
            if nanopore:
                vcf_path = (
                    "{}/articNcovNanopore_Genotyping_typeVariants/vcf/{}.vcf".format(
                        self.indir, base_sample
                    )
                )
            else:
                vcf_path = "{}/ncovIllumina_Genotyping_typeVariants/vcf/{}.vcf".format(
                    self.indir, base_sample
                )
            # Variants (vcf)
            deliv["files"].append(
                {
                    "format": "vcf",
                    "id": sampleID,
                    "path": vcf_path,
                    "path_index": "~",
                    "step": "genotyping",
                    "tag": "vcf-covid",
                }
            )
            if nanopore:
                fasta_path = "{}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka/{}.consensus.fasta".format(
                    self.indir, base_sample
                )
            else:
                fasta_path = "{}/ncovIllumina_sequenceAnalysis_makeConsensus/{}.consensus.fasta".format(
                    self.indir, base_sample
                )
            # Single-file fasta
            deliv["files"].append(
                {
                    "format": "fasta",
                    "id": sampleID,
                    "path": fasta_path,
                    "path_index": "~",
                    "step": "consensus",
                    "tag": "consensus-sample",
                }
            )
        with open(delivfile, "w") as out:
            yaml.dump(deliv, out)
