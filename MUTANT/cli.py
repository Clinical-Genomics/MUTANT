"""This is the main entry point of MUTANT.
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env pythons
import os
import sys
import click
import json

from MUTANT import __version__, preset_config, logger, wd

default_sampleinfo = {
    "CG_ID_project": "XXX0000",
    "CG_ID_sample": "XXX0000A1",
    "Customer_ID_project": "100100",
    "Customer_ID_sample": "10XY123456",
    "Customer_ID": "cust000",
    "application_tag": "SOMTIN100",
    "date_arrival": "0001-01-01 00:00:00",
    "date_libprep": "0001-01-01 00:00:00",
    "date_sequencing": "0001-01-01 00:00:00",
    "method_libprep": "Not in LIMS",
    "method_sequencing": "Not in LIMS",
    "organism": "SARS-CoV-2",
    "priority": "standard",
    "reference": "None",
}

if preset_config == "":
    click.echo(
        "ERROR - No properly set-up config under neither envvar MUTANT_CONFIG nor ~/.MUTANT/config.json. Exiting."
    )
    sys.exit(-1)

def set_cli_config(config):
    if config != "":
        if os.path.exists(config):
            try:
                with open(os.path.abspath(config), "r") as conf:
                    ctx.obj["config"] = json.load(conf)
                ctx.obj["config"]["config_path"] = os.path.abspath(config)
            except Exception as e:
                pass

def review_sampleinfo(pfile):
    """Reviews sample info. Returns loaded json object"""

    try:
        with open(pfile) as json_file:
            data = json.load(json_file)
    except Exception as e:
        click.echo("Unable to read provided sample info file as json. Exiting..")
        sys.exit(-1)

    if isinstance(data, list):
        for entry in data:
            for k, v in default_sampleinfo.items():
                if not k in entry:
                    click.echo(
                        "WARNING - Parameter {} needs to be provided in sample json. Formatting example: ({})".format(
                            k, v
                        )
                    )
    else:
        for k, v in default_sampleinfo.items():
            if not k in data:
                click.echo(
                    "WARNING - Parameter {} needs to be provided in sample json. Formatting example: ({})".format(
                        k, v
                    )
                )
    return data

def done():
    click.echo("INFO - Execution finished!")
    logger.debug("INFO - Execution finished!")

@click.group()
@click.version_option(__version__)
@click.pass_context
def root(ctx):
    """Microbial Utility Toolbox And wrapper for data traNsmission and Transformation """
    ctx.obj = {}
    ctx.obj["config"] = preset_config
    ctx.obj["log"] = logger

@root.group()
@click.pass_context
def analyse(ctx):
    """Run microbial analysis"""
    pass

@analyse.command()
@click.argument("sampleinfo_file")
@click.option("--input", help="Full path to input folder", default="")
@click.option("--config", help="MUTANT config to override default", default="")
@click.pass_context
def SARS_CoV_2(
    ctx, sampleinfo_file, input, config
):
    """Run SARS-CoV-2 analysis"""

    # Check configs
    set_cli_config(config)
    sampleinfo = review_sampleinfo(sampleinfo_file)

    # Check input data
    pool = []
    if not os.path.isdir(input):
        click.echo("ERROR - Sequence data folder {} does not exist.".format(input))
        ctx.abort()
    for subfolder in os.listdir(input):
        if os.path.isdir("{}/{}".format(input, subfolder)):
            pool.append(subfolder)
    run_settings = {
        "input": input,
        "pool": pool,
    }

    # Run
    """run_creator = Job_Creator(
        config=ctx.obj["config"],
        log=ctx.obj["log"],
        sampleinfo=sampleinfo,
        run_settings=run_settings,
    )"""

    done()