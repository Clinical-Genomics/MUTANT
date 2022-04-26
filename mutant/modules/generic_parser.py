import os
import sys
import json
import click
import pandas


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
