"""Parse SARS-CoV-2 results"""

import os
import sys
import json
import pandas
import glob
import csv
import re
from mutant import WD
from mutant.constants.artic import MULTIQC_TO_VOGUE, ILLUMINA_FILES_CASE
from mutant.modules.generic_parser import append_dict, parse_classifications


def get_results_paths(indir, case, ticket) -> dict:
    """Get paths for all reports and output files from GMS-Artic and MUTANT"""

    # Case paths
    path_dict = dict()
    case_paths = dict()
    for stepname, filepath in ILLUMINA_FILES_CASE.items():
        fullpath = filepath.format(resdir=indir, case=case, ticket=ticket)
        if "*" in fullpath:
            fullpath = glob.glob(fullpath)[0]
        case_paths[stepname] = fullpath
    path_dict[case] = case_paths
    return path_dict


def get_multiqc_template_data(record, field, step, direction="") -> dict:
    """Get default values for a multiqc metric"""

    info = MULTIQC_TO_VOGUE[step]
    sample_template = {
        "header": "~",
        "id": record["CG_ID_sample"],
        "input": info["input"].format(
            key="_".join(
                [
                    record["region_code"],
                    record["lab_code"],
                    record["Customer_ID_sample"],
                ]
            ),
            direction=direction,
        ),
        "name": info["fields"][field],
        "step": info["step"],
        "value": None,
    }
    return sample_template


def get_multiqc_sample_keys(data, sample, step, info) -> dict:
    """Get sample keys from multiqc.json data for an analysis step"""

    # Get sample keys
    sample_keys = {}
    # Get single key for sample
    if info["format"] == "single":
        for key in data["report_saved_raw_data"][step].keys():
            keycont = key.split("_")
            if sample in keycont:
                sample_keys["0"] = key
    # Get paired keys for sample
    elif info["format"] == "paired":
        for key in data["report_saved_raw_data"][step].keys():
            keycont = key.split("_")
            # Get keys containing sample name and "_1" or "_2" suffix
            if (sample in keycont) and (keycont[-1] in ["1", "2"]):
                suffix = keycont[-1]
                sample_keys[suffix] = key
    return sample_keys


def get_multiqc_metric_value(data, sample_keys, step, field, direction="0") -> str:
    """Get multiqc json data for metric value"""

    return data["report_saved_raw_data"][step][sample_keys[direction]][field]


def get_vogue_multiqc_data(multiqc, caseinfo) -> list:
    """Parse multiqc.json data for vogue delivery"""

    case_data = list()
    with open(multiqc) as f:
        data = json.load(f)
        for record in caseinfo:
            # sample_data = list()
            for step, info in MULTIQC_TO_VOGUE.items():
                # Get sample key naming in multiqc file
                if record["sequencing_qc_pass"]:
                    sample_keys = get_multiqc_sample_keys(
                        data, record["Customer_ID_sample"], step, info
                    )
                # Parse selected metrics in multiqc json
                for field, field_name in info["fields"].items():
                    # Add default values
                    if info["format"] == "single":
                        # Get template data
                        field_data = get_multiqc_template_data(record, field, step)
                        # Update value if analyzed
                        if record["sequencing_qc_pass"] and len(sample_keys) > 0:
                            field_data["value"] = get_multiqc_metric_value(
                                data, sample_keys, step, field
                            )
                        case_data.append(field_data)
                    else:
                        for direction in [1, 2]:
                            # Get template data
                            field_data = get_multiqc_template_data(
                                record, field, step, str(direction)
                            )
                            # Update value if analyzed
                            if record["sequencing_qc_pass"] and len(sample_keys) > 0:
                                field_data["value"] = get_multiqc_metric_value(
                                    data, sample_keys, step, field, str(direction)
                                )
                            case_data.append(field_data)
    return case_data


def get_artic_results(indir) -> dict:
    """Parse artic output directory for analysis results. Returns dictionary data object"""

    voc_pos = range(475, 486)
    muts = pandas.read_csv("{0}/standalone/spike_mutations.csv".format(WD), sep=",")
    # Magical unpacking into single list
    voc_pos_aa = sum(muts.values.tolist(), [])

    classifications_path = "{0}/standalone/classifications.csv".format(WD)
    voc_strains: dict = parse_classifications(csv_path=classifications_path)

    artic_data = dict()
    var_all = dict()
    var_voc = dict()

    # Files of interest. ONLY ADD TO END OF THIS LIST
    files = [
        "*qc.csv",
        "*variant_summary.csv",
    ]
    paths = list()
    for f in files:
        try:
            hits = glob.glob(os.path.join(indir, f))
            if len(hits) == 0:
                raise Exception("File not found")
            if len(hits) > 1:
                print(
                    "Multiple hits for {0}/{1}, picking {2}".format(indir, f, hits[0])
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
            if float(line[2]) >= 90:
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

    # Parse Pangolin report data
    pangodir = "{0}/ncovIllumina_sequenceAnalysis_pangolinTyping".format(indir)
    pangolins = glob.glob("{0}/*.pangolin.csv".format(pangodir))
    for path in pangolins:
        with open(path) as f:
            content = csv.reader(f)
            next(content)  # Skip header
            for line in content:
                sample = line[0].split(".")[0].split("_")[-1]
                artic_data[sample].update(
                    {
                        "lineage": line[1],
                        "pangolin_probability": line[3],
                        "pangolin_version": line[-7],
                        "pangolin_qc": line[-3],
                    }
                )

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
                    artic_data[sample].update({"variants": ";".join(var_all[sample])})
                else:
                    artic_data[sample].update({"variants": var_all[sample]})
            else:
                artic_data[sample].update({"variants": "-"})

    # Classification
    for sample, vals in artic_data.items():
        # Packing
        artic_data[sample].update({"VOC": "No"})

        # Check for lineage
        if artic_data[sample]["lineage"] == "None":
            artic_data[sample].update({"VOC": "-"})
        elif artic_data[sample]["lineage"] in voc_strains["lineage"]:
            index = voc_strains["lineage"].index(artic_data[sample]["lineage"])
            # Check for spike
            if (
                pandas.isna(voc_strains["spike"][index])
                or voc_strains["spike"][index] in artic_data[sample]["variants"]
            ):
                # Add variant class
                artic_data[sample].update({"VOC": voc_strains["class"][index]})
    return artic_data
