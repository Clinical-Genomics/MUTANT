""" This script creates a custom report from the QC output of the ARTIC pipeline,
    applied on SARS-CoV-2 samples 
    By: @talnor """
import os
import sys
import glob
from datetime import date
import csv
import re

if len(sys.argv) != 3:
    sys.exit(-1)

indir = sys.argv[1]
ticket = sys.argv[2]
voc_strain = ["A.23.1", "A.27", "A.28",
              "B.1.1.7", "B.1.351", "B.1.525", "B.1.526", "B.1.1.318", "B.1.427", "B.1.429",
              "C.16", "P.1", "P.2", "R.1", "R.2"]
voc_pos = range(475, 486)
voc_pos_aa = ["L452R"]

def create_summary(indir, voc_strain, voc_pos, voc_pos_aa):
    qcRep = glob.glob(os.path.join(indir, "*qc.csv"))[0]
    varRep = glob.glob(os.path.join(indir, "*variant_summary.csv"))[0]
    pangolinRep = glob.glob(os.path.join(indir, "ncovIllumina_sequenceAnalysis_makeConsensus/*pangolin.csv"))[0]
    res = dict()
    var = dict()
    # Parse qc report data
    with open(qcRep) as f:
        content = csv.reader(f)
        next(content)
        for line in content:
            sample = line[0].split("_")[2]
            res[sample] = [line[1], line[2], line[7]]
    # Parse Pangolin report data
    with open(pangolinRep) as f:
        content = csv.reader(f)
        next(content)
        for line in content:
            sample = line[0].split("_")[3].split(".")[0]
            res[sample].extend([line[1], line[3]])
    # Parse Variant report data
    with open(varRep) as f:
        content = csv.reader(f)
        next(content)
        for line in content:
            sample = line[0].split("_")[2]
            variant = line[2]
            pos = int(re.findall(r'\d+', variant)[0])
            if (pos in voc_pos) or (variant in voc_pos_aa):
                if sample in var.keys():
                    var[sample].append(variant)
                else:
                    var[sample] = [variant]

    # Write to summary report
    today = date.today().strftime("%Y%m%d")
    summaryfile = os.path.join(indir, "sars-cov-2_{}_results_{}.csv".format(ticket, today))
    with open(summaryfile, mode='w') as out:
        summary = csv.writer(out)
        summary.writerow(["Sample", "Selection", "Ticket", "%N_bases", "%covered_bases", "QC_pass", "Lineage",
                          "PangoLEARN_version", "Variants", "VOC"])
        for sample, data in res.items():
            selection = "-"
            data = [sample, selection, ticket] + data
            if sample in var.keys():
                data = data + [";".join(var[sample])]
            else:
                data = data + ["-"]
            if data[6] in voc_strain:
                data = data + ["Yes"]
            elif data[6] == "None":
                data = data + ["-"]
            else:
                data = data + ["No"]
            summary.writerow(data)
    # Write to variant report
    varfile = os.path.join(indir, "sars-cov-2_{}_variants_{}.csv".format(ticket, today))
    with open(varRep) as f, open(varfile, mode='w') as out:
        variants = f.readlines()
        varsummary = csv.writer(out)
        varsummary.writerow(variants[0].strip().split(","))
        for line in variants[1:]:
            line = line.strip().split(",")
            varsummary.writerow([line[0].split("_")[2]] + line[1:])
    print(varfile)    
    print(summaryfile)

create_summary(indir, voc_strain, voc_pos, voc_pos_aa)
