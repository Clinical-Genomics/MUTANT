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
from mutant.modules.generic_parser import (
    get_sarscov2_config,
    read_filelines,
)
from mutant.modules.artic_illumina.parser import (
    get_vogue_multiqc_data,
    get_artic_results,
    get_results_paths,
)
from mutant.modules.generic_reporter import GenericReporter


class ReportSC2:
    def __init__(self, caseinfo, indir, config_artic, fastq_dir, timestamp):
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
        self.filepaths = get_results_paths(self.indir, self.case, self.ticket)
        self.articdata = dict()
        self.consensus_path = "{0}/ncovIllumina_sequenceAnalysis_makeConsensus".format(
            self.indir
        )
        self.consensus_target_files = "{0}/*.consensus.fa".format(self.consensus_path)

    def create_all_files(self):
        generic_reporter = GenericReporter(
            caseinfo=self.caseinfo,
            indir=self.indir,
            nanopore=False,
        )
        generic_reporter.create_trailblazer_config()
        self.load_lookup_dict()
        self.create_concat_pangolin()
        self.create_concat_pangolin_fohm()
        # This works off concat pango and needs to occur after
        generic_reporter.create_concat_consensus(
            target_files=self.consensus_target_files
        )
        self.create_deliveryfile()
        self.create_vogue_metrics_file()
        self.create_fohm_csv()
        self.create_sarscov2_resultfile()
        self.create_sarscov2_variantfile()
        self.create_jsonfile()
        self.create_instrument_properties()

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

    def load_artic_results(self):
        """Parse artic output directory for analysis results. Returns dictionary data object"""
        artic_data = get_artic_results(self.indir)
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

    def create_concat_pangolin(self):
        """Concatenate pangolin results"""

        indir = "{0}/ncovIllumina_sequenceAnalysis_pangolinTyping".format(self.indir)
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

    def create_concat_pangolin_fohm(self):
        """Concatenate pangolin results and format for fohm"""

        indir = "{0}/ncovIllumina_sequenceAnalysis_pangolinTyping".format(self.indir)
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

    def create_instrument_properties(self):
        """Creates a properties file with instrument information"""

        propfile = os.path.join(self.indir, "instrument.properties")

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

    def create_sarscov2_resultfile(self):
        """Write summary csv report of Artic and Pangolin results"""

        ticket = self.ticket
        if self.articdata == dict():
            print("No artic results loaded. Quitting sarscov2_resultfile")
            sys.exit(-1)

        with open(self.filepaths[self.case]["results-file"], mode="w") as out:
            summary = csv.writer(out)
            summary.writerow(
                [
                    "Sample",
                    "Selection",
                    "Region Code",
                    "Ticket",
                    "%N_bases",
                    "%10X_coverage",
                    "QC_pass",
                    "Lineage",
                    "Pangolin_data_version",
                    "VOC",
                    "Mutations",
                ]
            )
            for sample, data in self.articdata.items():
                n_bases = tenx_bases = "0.0"
                qc_status = "FALSE"
                lineage = "None"
                verzion = "1970-01-01"
                vocs = vocs_aa = "-"
                selection = "-"

                if "selection_criteria" in data:
                    selection = data["selection_criteria"]
                if "region_code" in data:
                    region = data["region_code"]
                if "pct_n_bases" in data:
                    n_bases = data["pct_n_bases"]
                if "pct_10X_bases" in data:
                    tenx_bases = data["pct_10X_bases"]
                if "qc" in data:
                    if data["lineage"] == "None":
                        qc_status = "FALSE"
                    else:
                        qc_status = data["qc"]
                if "lineage" in data:
                    lineage = data["lineage"]
                if (
                    "pangolin_data_version" in data
                    and data["pangolin_data_version"] != ""
                ):
                    verzion = data["pangolin_data_version"]
                if "VOC" in data:
                    vocs = data["VOC"]
                if "VOC_aa" in data:
                    vocs_aa = data["VOC_aa"]

                row = [
                    sample,
                    selection,
                    region,
                    ticket,
                    n_bases,
                    tenx_bases,
                    qc_status,
                    lineage,
                    verzion,
                    vocs,
                    vocs_aa,
                ]

                summary.writerow(row)

    def create_sarscov2_variantfile(self):
        """Write variant csv report of identified variants
        I am literally just variant_summary.csv but with sample names"""

        indir = self.indir
        ticket = self.ticket
        today = self.today
        varRep = glob.glob(os.path.join(indir, "*variant_summary.csv"))[0]
        varout = os.path.join(indir, "sars-cov-2_{}_variants.csv".format(ticket, today))
        if os.stat(varRep).st_size != 0:
            with open(varRep) as f, open(varout, mode="w") as out:
                variants = f.readlines()
                varsummary = csv.writer(out)
                varsummary.writerow(variants[0].strip().split(","))
                for line in variants[1:]:
                    line = line.strip().split(",")
                    varsummary.writerow([line[0].split("_")[-1]] + line[1:])
        else:
            try:
                open(varout, "a").close()
            except Exception as e:
                print("Failed creating file {}\n{}".format(varout, e))

    def create_jsonfile(self):
        """Output all result data in a json format for easy parsing"""

        if self.articdata == dict():
            print("No artic results loaded. Quitting create_jsonfile")
            sys.exit(-1)

        with open(
            "{}/{}_artic.json".format(self.indir, self.ticket, self.today), "w"
        ) as outfile:
            json.dump(self.articdata, outfile)

    def create_deliveryfile(self):
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
            # Concat reads reverse
            deliv["files"].append(
                {
                    "format": "fastq",
                    "id": sampleID,
                    "path": "{0}/{1}_2.fastq.gz".format(self.fastq_dir, base_sample),
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
