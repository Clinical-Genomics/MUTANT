""" CLI access point for various MICROBIAL tools 
    By: Isak Sylvin @sylvinite & C:/O """


#!/usr/bin/env python

import os
import sys
import re
import json
import click
import subprocess
from datetime import datetime
from mutant import version, log
from mutant.assets.deliver.deliver_artic import DeliverSC2
from mutant.assets.run.run_artic import RunSC2
from mutant.assets.utils.parse import get_sarscov2_config

#File work directory
WD = os.path.dirname(os.path.realpath(__file__))
TIMESTAMP = datetime.now().strftime("%y%m%d-%H%M%S")

@click.group()
@click.version_option(version)
@click.pass_context
def root(ctx):
    """ Microbial Utility Toolbox And wrapper for data traNsmission and Transformation """
    ctx.obj = {}

@root.group()
@click.pass_context
def analyse(ctx):
    pass

@analyse.command()
@click.argument("input_folder")
@click.option("--config_artic", help="Custom artic configuration file", default="{}/config/hasta/artic.json".format(WD))
@click.option("--config_case", help="Provided config for the case", default="")
@click.option("--config_mutant", help="General configuration file for MUTANT", default="")
@click.option("--outdir", help="Output folder to override general configutations", default="")
@click.option("--profiles", help="Execution profiles, comma-separated", default="singularity,slurm")
@click.pass_context
def sarscov2(ctx, input_folder, config_artic, config_case, config_mutant, outdir, profiles):

    # Set base for output files (Move this section)
    if config_case != "":
        caseinfo = get_sarscov2_config(config_case)
        caseID = caseinfo[0]["case_ID"]
    else:
        caseID = "artic"
    prefix = "{}_{}".format(caseID, TIMESTAMP)

    # Run
    run = RunSC2(
        input_folder=input_folder,
        config_artic=config_artic,
        caseID=caseID,
        prefix=prefix,
        profiles=profiles,
        timestamp=TIMESTAMP,
        WD=WD
    )
    resdir = run.get_results_dir(config_mutant, outdir)
    run.run_case(resdir)

    # Deliverables
    if config_case != "":
        delivery = DeliverSC2(
            caseinfo=config_case,
            resdir=os.path.abspath(resdir),
            config_artic=config_artic,
            timestamp=TIMESTAMP
        )
        delivery.rename_deliverables()
        delivery.create_deliveryfile()
        delivery.create_fohm_csv()


@analyse.command()
@click.pass_context
def jasen(ctx):
    pass

@root.group()
@click.pass_context
def toolbox(ctx):
    pass

@toolbox.group()
@click.pass_context
def sarscov2(ctx):
    pass


@sarscov2.command()
@click.argument("input_folder")
@click.option("--config_artic", help="Custom artic configuration file", default="{}/config/hasta/artic.json".format(WD))
@click.option("--config_case", help="Provided config for the case", required=True)
@click.pass_context
def cgmodifications(ctx, input_folder, config_artic, config_case):
    """Applies all cg modifications as a batch"""


    # Deliverables
    if config_case != "":
        delivery = DeliverSC2(
            caseinfo=config_case,
            resdir=os.path.abspath(input_folder),
            config_artic=config_artic,
            timestamp=TIMESTAMP
        )
        delivery.rename_deliverables()
        delivery.create_deliveryfile()
        delivery.create_fohm_csv()

@sarscov2.command()
@click.argument("input_folder")
@click.option("--config_artic", help="Custom artic configuration file", default="{}/config/hasta/artic.json".format(WD))
@click.option("--config_case", help="Provided config for the case", required=True)
@click.pass_context
def rename(ctx, input_folder, config_artic, config_case):
    """Renames covid output to CG standard"""

    # Deliverables
    if config_case != "":
        delivery = DeliverSC2(
            caseinfo=config_case,
            resdir=os.path.abspath(input_folder),
            config_artic=config_artic,
            timestamp=TIMESTAMP
        )
        delivery.rename_deliverables()

@sarscov2.command()
@click.argument("input_folder")
@click.option("--config_artic", help="Custom artic configuration file", default="{}/config/hasta/artic.json".format(WD))
@click.option("--config_case", help="Provided config for the case", required=True)
@click.pass_context
def deliveryfile(ctx, input_folder, config_artic, config_case):
    """Generates CG specific delivery file"""

    # Deliverables
    if config_case != "":
        delivery = DeliverSC2(
            caseinfo=config_case,
            resdir=os.path.abspath(input_folder),
            config_artic=config_artic,
            timestamp=TIMESTAMP
        )
        delivery.create_deliveryfile()

@sarscov2.command()
@click.argument("input_folder")
@click.option("--config_artic", help="Custom artic configuration file", default="{}/config/hasta/artic.json".format(WD))
@click.option("--config_case", help="Provided config for the case", required=True)
@click.pass_context
def fohmfile(ctx, input_folder, config_artic, config_case):
    """Generates FoHM demanded delivery file"""


    # Deliverables
    if config_case != "":
        delivery = DeliverSC2(
            caseinfo=config_case,
            resdir=os.path.abspath(input_folder),
            config_artic=config_artic,
            timestamp=TIMESTAMP
        )
        delivery.create_fohm_csv()

@toolbox.command()
@click.argument("input_folder")
@click.argument("ticket_number")
def ArticReport(input_folder, ticket_number):
    """Report for QC output of the ARTIC pipeline"""
    cmd = "python {0}/assets/artic_report.py {1} {2}".format(WD, input_folder, ticket_number)
    log.debug("Command ran: {}".format(cmd))
    proc = subprocess.Popen(cmd.split())
    out, err = proc.communicate()
