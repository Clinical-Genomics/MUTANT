""" This class creates reports. Specifically it acts on the sarscov2 pipeline,
    and creates report files for Clinical Genomics Infrastructure

    By: Isak Sylvin & Tanja Normark
"""

import csv
import glob
import json
import pandas
import re
import os
import sys
import yaml 

from datetime import date
from pathlib import Path

from mutant import WD
from mutant.modules.generic_parser import get_sarscov2_config, get_json, append_dict


class ReportSC2:
    def __init__(self, caseinfo, indir, config_artic, fastq_dir, timestamp):
        self.casefile = caseinfo
        caseinfo = get_sarscov2_config(caseinfo)

        regionlab_list = []
        for record in caseinfo:
            regionlab = "{}_{}".format(record["region_code"], record["lab_code"])
            if regionlab not in regionlab_list:
                regionlab_list.append(regionlab)
        self.caseinfo = caseinfo
        self.case = caseinfo[0]["case_ID"]
        self.ticket = caseinfo[0]["Customer_ID_project"]
        self.project = caseinfo[0]["Customer_ID_project"]
        self.regionlabs = regionlab_list
        self.indir = indir
        self.config_artic = config_artic
        self.time = timestamp
        today = date.today().strftime("%Y%m%d")
        self.today = today
        self.fastq_dir = fastq_dir
        self.articdata = dict()

    def create_all_files(self):
        self.create_trailblazer_config()
        self.create_concat_pangolin()
        self.create_concat_consensus()
        self.create_deliveryfile()
        self.create_fohm_csv()
        self.load_lookup_dict()
        self.create_sarscov2_resultfile()
        self.create_sarscov2_variantfile()
        self.create_jsonfile()

    def get_finished_slurm_ids(self) -> list:
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
        trailblazer_config_path = Path(self.indir, "trailblazer_config.yaml")
        finished_slurm_ids = self.get_finished_slurm_ids()
        if not finished_slurm_ids:
            return
        with open(trailblazer_config_path, "w") as trailblazer_config_file:
            yaml.dump(data={"jobs": finished_slurm_ids}, stream=trailblazer_config_file)

    def create_concat_pangolin(self):

        indir = "{0}/ncovIllumina_sequenceAnalysis_pangolinTyping".format(self.indir)

        concat = open("{0}/{1}.pangolin.csv".format(self.indir, self.ticket), "w+")
        pangolins = glob.glob("{0}/*.pangolin.csv".format(indir))
        with open(pangolins[0], "r") as f:
            header = f.readlines()[0]
        concat.write(header)
        for item in pangolins:
            single = open(item, "r")
            concat.write("\n".join(single.readlines()[1:]))
        concat.close()

    def create_concat_consensus(self):

        indir = "{0}/ncovIllumina_sequenceAnalysis_makeConsensus".format(self.indir)

        concat = open("{0}/{1}.consensus.fa".format(self.indir, self.ticket), "w+")
        for item in glob.glob("{0}/*.consensus.fa".format(indir)):
            single = open(item, "r")
            concat.write(single.read())
            concat.write("\n")
        concat.close()

    def create_fohm_csv(self):

        """Creates a summary file for FoHM for each region-lab-combination"""

        # Add header to summary files
        for regionlab in self.regionlabs:
            sumfile = os.path.join(
                self.indir,
                "{}_{}_komplettering.csv".format(regionlab, self.today),
            )
            with open(sumfile, "w") as summary:
                summary.write("provnummer,urvalskriterium,GISAID_accession\n")

        # Write sample information to corresponding summary file
        for record in self.caseinfo:
            region = record["region_code"]
            lab = record["lab_code"]
            sumfile = os.path.join(
                self.indir, "{}_{}_{}_komplettering.csv".format(region, lab, self.today)
            )
            with open(sumfile, "a") as out:
                summary = csv.writer(out)
                summary.writerow(
                    [
                        record["Customer_ID_sample"],
                        record["selection_criteria"],
                    ]
                )

    def create_sarscov2_resultfile(self):
        """Write summary csv report of Artic and Pangolin results"""

        ticket = self.ticket
        today = self.today
        if self.articdata == dict():
            print("No artic results loaded. Quitting sarscov2_resultfile")
            sys.exit(-1)
        indir = self.indir

        summaryfile = os.path.join(
            indir, "sars-cov-2_{}_results.csv".format(ticket, today)
        )
        with open(summaryfile, mode="w") as out:
            summary = csv.writer(out)
            summary.writerow(
                [
                    "Sample",
                    "Selection",
                    "Ticket",
                    "%N_bases",
                    "%10X_coverage",
                    "QC_pass",
                    "Lineage",
                    "PangoLEARN_version",
                    "VOC",
                ]
            )
            for sample, data in self.articdata.items():
                selection = "-"
                row = [
                    sample,
                    data["selection_criteria"],
                    ticket,
                    data["pct_n_bases"],
                    data["pct_10X_bases"],
                    data["qc"],
                    data["lineage"],
                    data["pangoLEARN_version"],
                    data["VOC"],
                ]
#        with open(summaryfile, mode="w") as out:
#            summary = csv.writer(out)
#            summary.writerow(
#                [
#                    "Sample",
#                    "Selection",
#                    "Region Code",
#                    "Ticket",
#                    "%N_bases",
#                    "%10X_coverage",
#                    "QC_pass",
#                    "Lineage",
#                    "PangoLEARN_version",
#                    "VOC",
#                    "Mutations",
#                ]
#            )
#            for sample, data in self.articdata.items():
#                selection = "-"
#                row = [
#                    sample,
#                    data["selection_criteria"],
#                    data["region_code"],
#                    ticket,
#                    data["pct_n_bases"],
#                    data["pct_10X_bases"],
#                    data["qc"],
#                    data["lineage"],
#                    data["pangoLEARN_version"],
#                    data["VOC"],
#                    data["VOC_aa"],
#                ]
                summary.writerow(row)

    def create_sarscov2_variantfile(self):
        """Write variant csv report of identified variants
        I am literally just variant_summary.csv but with sample names"""

        indir = self.indir
        ticket = self.ticket
        today = self.today

        varRep = glob.glob(os.path.join(indir, "*variant_summary.csv"))[0]
        varout = os.path.join(
            indir, "sars-cov-2_{}_variants.csv".format(ticket, today)
        )
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



    def load_lookup_dict(self):
        """ Loads articdata with data from various sources. Atm, artic output and the case          config input file """
        self.load_artic_results()
        self.load_case_config()

    def load_case_config(self):
        """ Appends additional data to articdata dictionary """
        casekeys = self.caseinfo[0].keys()

        packing = dict(zip(casekeys, "-"*len(casekeys)))

        #Packs with keys. Time consuming but not really
        for k, v in self.articdata.items():
            self.articdata[k].update(packing)
        #Writes caseconfig data where relevant
        for entry in self.caseinfo:
            k = entry['Customer_ID_sample']
            if k in self.articdata.keys():
                self.articdata[k].update(entry)


    def load_artic_results(self):
        """Parse artic output directory for analysis results. Returns dictionary data object        """
        indir = self.indir
        voc_pos = range(475, 486)
        muts = pandas.read_csv("{0}/standalone/spike_mutations.csv".format(WD), sep=",")
        # Magical unpacking into single list
        voc_pos_aa = sum(muts.values.tolist(), [])

        classifications = pandas.read_csv("{0}/standalone/classifications.csv".format(WD), sep=",")
        voc_strains = { 'lineage':'','spike':'','class':''}
        voc_strains['lineage'] = classifications['lineage'].tolist()
        voc_strains['spike'] = classifications['spike'].tolist()
        voc_strains['class'] = classifications['class'].tolist()

        #voc_pos_aa = get_json("{0}/standalone/voc_strains.json".format(WD))['voc_pos_aa']
        #voc_strains = get_json("{0}/standalone/voc_strains.json".format(WD))['voc_strains']

        artic_data = dict()
        var_all = dict()
        var_voc = dict()

        # Files of interest. ONLY ADD TO END OF THIS LIST
        files = [
            "*qc.csv",
            "*variant_summary.csv",
            "*pangolin.csv",
        ]
        paths = list()
        for f in files:
            try:
                hits = glob.glob(os.path.join(indir, f))
                if len(hits) == 0:
                    raise Exception("File not found")
                if len(hits) > 1:
                    print(
                        "Multiple hits for {0}/{1}, picking {2}".format(
                            indir, f, hits[0]
                        )
                    )
                paths.append(hits[0])
            except Exception as e:
                print("Unable to find {0} in {1} ({2})".format(f, indir, e))
                sys.exit(-1)

        # Parse qc report data
        with open(paths[0]) as f:
            content = csv.reader(f)
            next(content)
            for line in content:
                sample = line[0].split("_")[-1]
                if float(line[2]) > 95:
                    qc_flag = "TRUE"
                else:
                    qc_flag = "FALSE"
                artic_data[sample] = {
                    "pct_n_bases": line[1],
                    "pct_10X_bases": line[2],
                    "longest_no_N_run": line[3],
                    "num_aligned_reads": line[4],
                    "artic_qc": line[7],
                    "qc": qc_flag,
                }
        # Parse Pangolin report data
        with open(paths[2]) as f:
            content = csv.reader(f)
            next(content)
            for line in content:
                sample = line[0].split(".")[0].split("_")[-1]
                lineage = line[1]

                artic_data[sample].update(
                    {
                        "lineage": lineage,
                        "pangolin_probability": line[2],
                        "pangoLEARN_version": line[3],
                        "pangolin_qc": line[4],
                    }
                )

        # Parse Variant report data
        if os.stat(paths[1]).st_size != 0:
            with open(paths[1]) as f:
                content = csv.reader(f)
                next(content)
                for line in content:
                    sample = line[0].split("_")[-1]
                    variant = line[2]
                    pos = int(re.findall(r"\d+", variant)[0])
                    if (pos in voc_pos) or (variant in voc_pos_aa):
                        append_dict(var_voc, sample, variant)
                    append_dict(var_all, sample, variant)
        # Add variant data to results
        if var_voc:
            for sample in artic_data.keys():
                if sample in var_voc.keys():
                    artic_data[sample].update({"VOC_aa": ";".join(var_voc[sample])})
                else:
                    artic_data[sample].update({"VOC_aa": "-"})
        else:
            for sample in artic_data.keys():
                artic_data[sample].update({"VOC_aa": "-"})
        if var_all:
            for sample in artic_data.keys():
                if sample in var_all.keys():
                    if len(var_all[sample]) > 1:
                        artic_data[sample].update(
                            {"variants": ";".join(var_all[sample])}
                        )
                    else:
                        artic_data[sample].update({"variants": var_all[sample]})
                else:
                    artic_data[sample].update({"variants": "-"})


        #Classification
        for key, vals in artic_data.items():
            #Packing
            artic_data[key].update( {"VOC": "No" } )

            #Check for lineage
            if artic_data[key]["lineage"] == "None":
                artic_data[key].update( {"VOC":"-"})
            elif artic_data[key]["lineage"] in voc_strains['lineage']:
                index = voc_strains['lineage'].index(artic_data[key]['lineage'])
                #Check for spike
                if pandas.isna(voc_strains['spike'][index]) or voc_strains['spike'][index] in artic_data[key]['VOC_aa']:
                    artic_data[key].update( {"VOC":voc_strains['class'][index]} )



        self.articdata.update(artic_data)

    def create_deliveryfile(self):

        """Create deliverables file"""

        deliv = {"files": []}
        delivfile = "{}/{}_deliverables.yaml".format(self.indir, self.case)

        ## Per Case
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
        # MultiQC json (vogue)
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
 
        ## Artic Summary report
        #deliv["files"].append(
        #    {
        #        "format": "csv",
        #        "id": self.case,
        #        "path": "{}/{}.typing_summary.csv".format(self.indir, self.ticket),
        #        "path_index": "~",
        #        "step": "report",
        #        "tag": "artic-sum",
        #    }
        #)
        ## Artic Variant report
        #deliv["files"].append(
        #    {
        #        "format": "csv",
        #        "id": self.case,
        #        "path": "{}/{}.variant_summary.csv".format(self.indir, self.ticket),
        #        "path_index": "~",
        #        "step": "report",
        #        "tag": "artic-var",
        #    }
        #)
        ## Artic QC report
        #deliv["files"].append(
        #    {
        #        "format": "csv",
        #        "id": self.case,
        #        "path": "{}/{}.qc.csv".format(self.indir, self.ticket),
        #        "path_index": "~",
        #        "step": "result_aggregation",
        #        "tag": "artic-qc",
        #    }
        #)

        # Artic Json (Vogue) data
        deliv["files"].append(
            {
                "format": "json",
                "id": self.case,
                "path": "{}/{}_artic.json".format(self.indir, self.ticket),
                "path_index": "~",
                "step": "result_aggregation",
                "tag": "artic-json",
            }
        )
        # Provided CG CASE info from StatusDB
        deliv["files"].append(
            {
                "format": "json",
                "id": self.case,
                "path": "{}".format(self.casefile),
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
                "path": self.config_artic,
                "path_index": "~",
                "step": "runinfo",
                "tag": "runtime-settings",
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

        # Region split
        for regionlab in self.regionlabs:
            rl = regionlab
            # Region split Pangolin typing
#            deliv["files"].append(
#                {
#                    "format": "csv",
#                    "id": self.case,
#                    "path": "{}/ncovIllumina_sequenceAnalysis_makeConsensus/"
#                    "{}_{}_pangolin_classification.txt".format(
#                        self.indir, rl, self.today
#                    ),
#                    "path_index": "~",
#                    "step": "typing",
#                    "tag": "SARS-CoV-2-type",
#                }
#            )
            # Region split FoHM delivery file
            deliv["files"].append(
                {
                    "format": "csv",
                    "id": self.case,
                    "path": os.path.join(
                        self.indir, "{}_{}_komplettering.csv".format(rl, self.today)
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

            ## Commenting these to save space in CG. Can be reenabled dynamically

            ## Alignment (bam, sorted)
            #deliv["files"].append(
            #    {
            #        "format": "bam",
            #        "id": sampleID,
            #        "path": "{}/ncovIllumina_sequenceAnalysis_readMapping/{}.sorted.bam".format(
            #            self.indir, base_sample
            #        ),
            #        "path_index": "~",
            #        "step": "alignment",
            #        "tag": "reference-alignment-sorted",
            #    }
            #)
            ## Variants (csq-vcf)
            #deliv["files"].append(
            #    {
            #        "format": "vcf",
            #        "id": sampleID,
            #        "path": "{}/ncovIllumina_Genotyping_typeVariants/vcf/{}.vcf".format(
            #            self.indir, base_sample
            #        ),
            #        "path_index": "~",
            #        "step": "genotyping",
            #        "tag": "variants",
            #    }
            #)
            ## Variants (tsv)
            #deliv["files"].append(
            #    {
            #        "format": "tsv",
            #        "id": sampleID,
            #        "path": "{}/ncovIllumina_sequenceAnalysis_callVariants/{}.variants.tsv".format(
            #            self.indir, base_sample
            #        ),
            #        "path_index": "~",
            #        "step": "variant-calling",
            #        "tag": "variants",
            #    }
            #)

        with open(delivfile, "w") as out:
            yaml.dump(deliv, out)
