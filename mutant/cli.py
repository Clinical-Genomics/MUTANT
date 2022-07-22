""" CLI access point for various MICROBIAL tools 
    By: Isak Sylvin @sylvinite & C:/O """

# !/usr/bin/env python

import os
import subprocess

import click

from mutant import version, log, WD, TIMESTAMP
from mutant.modules.artic_illumina.start import RunSC2
from mutant.modules.artic_nanopore.reformat import reformat_fastq_folder
from mutant.modules.artic_nanopore.report import ReportPrinterNanopore
from mutant.modules.generic_parser import get_json
from mutant.modules.artic_illumina.report import ReportSC2
from mutant.modules.artic_illumina.delivery import DeliverySC2


@click.group()
@click.version_option(version)
@click.pass_context
def root(ctx):
    """Microbial Utility Toolbox And wrapper for data traNsmission and Transformation"""
    ctx.obj = {}


@root.group()
@click.pass_context
def analyse(ctx):
    pass


@analyse.command()
@click.argument("input_folder")
@click.option(
    "--config_artic",
    help="Custom artic configuration file",
)
@click.option("--config_case", help="Provided config for the case")
@click.option(
    "--config_mutant",
    help="General configuration file for MUTANT",
    default="{}/config/hasta/mutant.json".format(WD),
)
@click.option(
    "--outdir", help="Output folder to override general configutations", default=""
)
@click.option("-n", "--nanopore", is_flag=True)
@click.option(
    "--profiles",
    help="Execution profiles, comma-separated",
    default="singularity,slurm",
)
@click.pass_context
def sarscov2(
    ctx,
    input_folder,
    config_artic,
    config_case,
    config_mutant,
    outdir,
    profiles,
    nanopore,
):

    if not config_case:
        log.debug("No case config given, mutant needs a case config to run analysis")

    caseinfo = get_json(config_case)
    caseID = caseinfo[0]["case_ID"]
    prefix = "{}_{}".format(caseID, TIMESTAMP)

    primer_json = caseinfo[0]["primer"].replace(" ", "_")
    if not config_artic:
        primer_config = "{}/config/hasta/" + primer_json.lower() + ".json"
        config_artic = primer_config.format(WD)

    if "nanopore" in primer_json.lower():
        nanopore = True

    # Run
    run = RunSC2(
        input_folder=input_folder,
        caseID=caseID,
        config_artic=config_artic,
        prefix=prefix,
        profiles=profiles,
        timestamp=TIMESTAMP,
        WD=WD,
    )

    resdir = run.get_results_dir(config_mutant, outdir)
    if nanopore:
        barcode_to_sampleid: dict = reformat_fastq_folder(
            fastq_folder=os.path.abspath(input_folder)
        )
    run.run_case(resdir, nanopore)

    if nanopore:
        if config_case != "":
            nanopore_report = ReportPrinterNanopore(
                caseinfo=config_case,
                indir=os.path.abspath(resdir),
                barcode_to_sampleid=barcode_to_sampleid,
                config_artic=config_artic,
            )
            nanopore_report.create_all_nanopore_files()

    else:
        # Report
        if config_case != "":
            report = ReportSC2(
                caseinfo=config_case,
                indir=os.path.abspath(resdir),
                fastq_dir=os.path.abspath(input_folder),
                config_artic=config_artic,
                timestamp=TIMESTAMP,
            )
            report.create_all_files()

        # Deliverables
        if config_case != "":
            delivery = DeliverySC2(
                caseinfo=config_case,
                indir=os.path.abspath(resdir),
            )

            delivery.rename_deliverables()


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
@click.option(
    "--config_artic",
    help="Custom artic configuration file",
    default="{}/config/hasta/default_config.json".format(WD),
)
@click.option("--fastq_folder", help="Sequence data folder for the case", required=True)
@click.option("--config_case", help="Provided config for the case", required=True)
@click.pass_context
def postproc(ctx, input_folder, config_artic, fastq_folder, config_case):
    """Applies all cg post-processing of the sarscov2 pipeline"""

    # Reports
    if config_case != "":
        report = ReportSC2(
            caseinfo=config_case,
            indir=os.path.abspath(input_folder),
            config_artic=config_artic,
            fastq_dir=os.path.abspath(fastq_folder),
            timestamp=TIMESTAMP,
        )

        report.create_all_files()

    # Deliverables
    if config_case != "":
        delivery = DeliverySC2(
            caseinfo=config_case,
            indir=os.path.abspath(input_folder),
        )

        delivery.rename_deliverables()


@sarscov2.command()
@click.argument("input_folder")
@click.option(
    "--config_artic",
    help="Custom artic configuration file",
    default="{}/config/hasta/default_config.json".format(WD),
)
@click.option("--config_case", help="Provided config for the case", required=True)
@click.pass_context
def rename(ctx, input_folder, config_artic, config_case):
    """Renames sarcov2 pipeline output to CG standard"""

    # Delivery
    if config_case != "":
        delivery = DeliverySC2(
            caseinfo=config_case,
            indir=os.path.abspath(input_folder),
        )
        delivery.rename_deliverables()


@toolbox.command()
@click.option("--input_folder", help="Folder with fastq to concatenate", required=True)
@click.option(
    "--app_tag", help="Application tag", required=False, default="CONCATENATE"
)
@click.option(
    "--date",
    help="Date to add to the concatenated file name, e.g. order date",
    required=False,
    default=None,
)
@click.pass_context
def concatenate(ctx, input_folder, app_tag, date):
    """Concatenates fastq files if needed"""
    if date:
        cmd = "python {0}/standalone/concatenate.py --input_folder {1} --app_tag {2} --date {3}".format(
            WD, input_folder, app_tag, date
        )
    else:
        cmd = "python {0}/standalone/concatenate.py --input_folder {1} --app_tag {2}".format(
            WD, input_folder, app_tag
        )
    log.debug("Command ran: {}".format(cmd))
    proc = subprocess.Popen(cmd.split())
    out, err = proc.communicate()


@toolbox.command()
@click.pass_context
def create_images(ctx):
    """Builds the sarscov2 pipeline images"""
    bdir = os.getcwd()
    os.chdir("{0}/externals/gms-artic".format(WD))
    cmd = "bash scripts/build_singularity_containers.sh && chmod 0777 *.sif"
    log.debug("Command ran: {}".format(cmd))
    proc = subprocess.Popen(cmd.split())
    out, err = proc.communicate()
    os.chdir(bdir)
