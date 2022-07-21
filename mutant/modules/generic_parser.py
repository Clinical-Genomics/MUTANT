import os
import sys
import json
import click
import pandas
import csv
import glob
import re

from mutant import WD


def get_artic_results(indir, nanopore) -> dict:
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
    if nanopore:
        pangodir = "{0}/articNcovNanopore_sequenceAnalysisMedaka_pangolinTyping".format(
            indir
        )
    else:
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
                        "pangolin_data_version": line[8],
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


def append_dict(dictionary, key, item) -> dict:
    if key in dictionary.keys():
        dictionary[key].append(item)
    else:
        dictionary[key] = [item]
    return dictionary


def get_json(config) -> dict:
    """Read json file"""
    if os.path.exists(config):
        try:
            with open(config) as json_file:
                data = json.load(json_file)
        except Exception as e:
            click.echo(
                "Unable to read provided json file: {}. Exiting..".format(config)
            )
            click.echo(e)
            sys.exit(-1)
    else:
        click.echo("Could not find supplied config: {}. Exiting..".format(config))
        sys.exit(-1)
    return data


def get_sarscov2_config(config) -> dict:
    """Parse SARS-CoV-2 sample config"""
    caseinfo = get_json(config)
    for i in range(len(caseinfo)):
        caseinfo[i]["region_code"] = caseinfo[i]["region_code"].replace(" ", "_")
        caseinfo[i]["lab_code"] = caseinfo[i]["lab_code"].replace(" ", "_")
        selection_criteria = caseinfo[i]["selection_criteria"]
        first_character = selection_criteria[:1]
        if first_character.isnumeric():
            caseinfo[i]["selection_criteria"] = (
                caseinfo[i]["selection_criteria"].split(".")[1].strip()
            )
    return caseinfo


def read_filelines(infile) -> list:
    try:
        with open(infile, "r") as f:
            contents = f.readlines()
    except Exception as e:
        click.echo("Unable to read file: {}. Exiting..".format(infile))
        click.echo(e)
        sys.exit(-1)
    return contents


def parse_classifications(csv_path: str) -> dict:
    """Parse classifications.csv, which contains info about VOC/VOI"""
    classifications = pandas.read_csv(csv_path, sep=",")
    voc_strains = {"lineage": "", "spike": "", "class": ""}
    voc_strains["lineage"] = classifications["lineage"].tolist()
    voc_strains["spike"] = classifications["spike"].tolist()
    voc_strains["class"] = classifications["class"].tolist()
    return voc_strains
