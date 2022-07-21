""" This class creates reports. Specifically it acts on the sarscov2 pipeline,
    and creates report files for Clinical Genomics Infrastructure

    By: Isak Sylvin & Tanja Normark
"""
import os
import sys
import csv
import glob
import json
from datetime import date
from mutant.modules.generic_parser import get_sarscov2_config
from mutant.modules.generic_reporter import get_results_paths


class IlluminaReporter:
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

    def create_reports(self):
        self.create_sarscov2_resultfile()
        self.create_sarscov2_variantfile()

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
