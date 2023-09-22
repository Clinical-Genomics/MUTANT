import glob
import os
import sys
import json
import click
import pandas

from mutant.constants.artic import ILLUMINA_FILES_CASE, NANOPORE_FILES_CASE


def get_results_paths(indir, case, ticket, nanopore) -> dict:
    """Get paths for all reports and output files from GMS-Artic and MUTANT"""

    # Case paths
    path_dict = dict()
    case_paths = dict()
    if nanopore:
        files_to_iterate = NANOPORE_FILES_CASE
    else:
        files_to_iterate = ILLUMINA_FILES_CASE
    for stepname, filepath in files_to_iterate.items():
        fullpath = filepath.format(resdir=indir, case=case, ticket=ticket)
        if "*" in fullpath:
            fullpath = glob.glob(fullpath)[0]
        case_paths[stepname] = fullpath
    path_dict[case] = case_paths
    return path_dict


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
            click.echo(f"Unable to read provided json file: {config}. Exiting..")
            click.echo(e)
            sys.exit(-1)
    else:
        click.echo(f"Could not find supplied config: {config}. Exiting..")
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
        click.echo(f"Unable to read file: {infile}. Exiting..")
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
