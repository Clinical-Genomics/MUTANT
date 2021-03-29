""" This script renames the output files for the ARTIC pipeline and Pangolin analysis,
    and creates a deliverables file for Clinical Genomics Infrastructure

By: Isak Sylvin & Tanja Normark
"""

from datetime import date

import csv
import glob
import os
import json
import yaml


def get_json_data(config):
    if os.path.exists(config):
        """Get sample information as json object"""
        try:
            with open(config) as json_file:
                data = json.load(json_file)
        except Exception as e:
            click.echo("Unable to read provided json file: {}. Exiting..".format(config))
            click.echo(e)
            sys.exit(-1)
    else:
        click.echo("Could not find supplied config: {}. Exiting..".format(config))
        sys.exit(-1)
    return data

class DeliverSC2:

    def __init__(self, caseinfo, resdir, config_artic, timestamp):
        self.casefile = caseinfo
        caseinfo = get_json_data(caseinfo)

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
        self.indir = resdir
        self.config_artic = config_artic
        self.time = timestamp
        today = date.today().strftime("%Y%m%d")
        self.today = today

    def rename_deliverables(self):
        """Rename result files for delivery: fastq, consensus files, vcf and pangolin"""

        for sampleinfo in self.caseinfo:
            sample = sampleinfo["CG_ID_sample"]
            region = sampleinfo["region_code"]
            lab = sampleinfo["lab_code"]

            #rename makeConsensus
            prefix = "{0}/ncovIllumina_sequenceAnalysis_makeConsensus".format(self.indir)
            for item in glob.glob("{0}/*{1}*".format(prefix, sample)):
                newpath = "{0}/{1}_{2}_{3}.consensus.fasta".format(prefix, region, lab, sampleinfo["Customer_ID_sample"])
                os.rename(item, newpath)

            #rename typeVariants
            prefix = "{0}/ncovIllumina_Genotyping_typeVariants/vcf".format(self.indir)
            for item in glob.glob("{0}/*{1}*.csq.vcf".format(prefix, sample)):
                newpath = "{0}/{1}_{2}_{3}.vcf".format(prefix, region, lab, sampleinfo["Customer_ID_sample"])
                os.rename(item, newpath)

            #rename callVariants
            prefix = "{0}/ncovIllumina_sequenceAnalysis_callVariants".format(self.indir)
            for item in glob.glob("{0}/*{1}*.variants.tsv".format(prefix, sample)):
                newpath = "{0}/{1}_{2}_{3}.variants.tsv".format(prefix, region, lab, sampleinfo["Customer_ID_sample"])
                os.rename(item, newpath)

            #rename ReadMapping
            prefix = "{0}/ncovIllumina_sequenceAnalysis_readMapping".format(self.indir)
            for item in glob.glob("{0}/*{1}*.sorted.bam".format(prefix, sample)):
                newpath = "{0}/{1}.sorted.bam".format(prefix, sampleinfo["Customer_ID_sample"])
                os.rename(item, newpath)

            #rename core
            core_suffix = ['.qc.csv','.pangolin.csv','.typing_summary.csv','.variant_summary.csv']
            for thing in core_suffix:
                hit = glob.glob("{0}/*{1}".format(self.indir, thing))
                if len(hit) == 1:
                    hit=hit[0]
                    os.rename(hit, "{0}/{1}{2}".format(self.indir, self.ticket, thing))
           
    def create_deliveryfile(self):

        """Create deliverables file"""

        deliv = dict()
        deliv['files'] = list()
        delivfile = "{}/{}_deliverables.yaml".format(self.indir, self.case)

        # Project-wide

        # Summary report
        deliv['files'].append({'format': 'csv', 'id': self.case,
                               'path': "{}/{}.typing_summary.csv".format(self.indir, self.ticket),
                               'path_index': '~', 'step': 'report', 'tag': 'SARS-CoV-2-sum'})
        # Variant report
        deliv['files'].append({'format': 'csv', 'id': self.case,
                               'path': "{}/{}.variant_summary.csv".format(self.indir, self.ticket),
                               'path_index': '~', 'step': 'report', 'tag': 'SARS-CoV-2-var'})
        # QC report
        deliv['files'].append({'format': 'csv', 'id': self.case,
                               'path': "{}/{}.qc.csv".format(self.indir, self.ticket),
                               'path_index': '~', 'step': 'result_aggregation', 'tag': 'SARS-CoV-2-qc'})
        # Json (vogue) data
        deliv['files'].append({'format': 'json', 'id': self.case,
                               'path': "{}/{}.json".format(self.indir, self.case),
                               'path_index': '~', 'step': 'result_aggregation', 'tag': 'SARS-CoV-2-json'})
        # Sampleinfo/Case
        deliv['files'].append({'format': 'json', 'id': self.case,
                               'path': "{}".format(self.casefile),
                               'path_index': '~', 'step': 'runinfo', 'tag': 'sampleinfo'})
        # Settings dump
        deliv['files'].append({'format': 'txt', 'id': self.case,
                               'path': self.config_artic,
                               'path_index': '~', 'step': 'runinfo', 'tag': 'runtime-settings'})
        # Execution log
        deliv['files'].append({'format': 'txt', 'id': self.case,
                               'path': "{}/nextflow.log".format(self.indir),
                               'path_index': '~', 'step': 'runinfo', 'tag': 'logfile'})

        # Region-lab-wide

        for regionlab in self.regionlabs:
            rl = regionlab.replace(' ', "_")
            # Pangolin reports
            deliv['files'].append({'format': 'csv', 'id': self.case,
                                   'path': "{}/ncovIllumina_sequenceAnalysis_makeConsensus/"
                                           "{}_{}_pangolin_classification.txt".format(self.indir, rl,
                                                                                      self.today),
                                   'path_index': '~', 'step': 'typing', 'tag': 'SARS-CoV-2-type'})
            # FoHM delivery file
            deliv['files'].append({'format': 'csv', 'id': self.case,
                                   'path': os.path.join(self.indir,
                                                        "{}_{}_komplettering.csv".format(rl, self.today)),
                                   'path_index': '~', 'step': 'report', 'tag': 'SARS-CoV-2-info'})

        # Sample-wide

        for record in self.caseinfo:
            sampleID = record["CG_ID_sample"]
            sample = record["Customer_ID_sample"]
            region = record["region_code"].replace(' ', '_')
            lab = record["lab_code"].replace(' ', '_')
            base_sample = "{0}_{1}_{2}".format(region, lab, sample)
            # Concat reads forwards
            deliv['files'].append({'format': 'fastq', 'id': sampleID,
                                   'path': "{0}/ncovIllumina_sequenceConcatination/{1}_1.fastq.gz".format(self.indir, base_sample),
                                   'path_index': '~', 'step': 'concatination', 'tag': 'forward-reads'})
            # Concat reads reverse
            deliv['files'].append({'format': 'fastq', 'id': sampleID,
                                   'path': "{0}/ncovIllumina_sequenceConcatination/{1}_2.fastq.gz".format(self.indir, base_sample),
                                   'path_index': '~', 'step': 'concatination', 'tag': 'reverse-reads'})
            # Consensus files
            deliv['files'].append({'format': 'fasta', 'id': sampleID,
                                   'path': "{}/ncovIllumina_sequenceAnalysis_makeConsensus/{}.consensus.fasta"
                                   .format(self.indir, base_sample),
                                   'path_index': '~', 'step': 'consensus', 'tag': 'consensus'})
            # Alignment (bam, sorted)
            deliv['files'].append({'format': 'bam', 'id': sampleID,
                                   'path': "{}/ncovIllumina_sequenceAnalysis_readMapping/{}.sorted.bam".format(
                                       self.indir, sample),
                                   'path_index': '~', 'step': 'alignment', 'tag': 'reference-alignment-sorted'})
            # Variants (csq-vcf)
            deliv['files'].append({'format': 'vcf', 'id': sampleID,
                                   'path': "{}/ncovIllumina_Genotyping_typeVariants/vcf/{}.vcf".format(self.indir,
                                                                                                       base_sample),
                                   'path_index': '~', 'step': 'genotyping', 'tag': 'variants'})
            # Variants (tsv)
            deliv['files'].append({'format': 'tsv', 'id': sampleID,
                                   'path': "{}/ncovIllumina_sequenceAnalysis_callVariants/{}.variants.tsv".format(
                                       self.indir, base_sample),
                                   'path_index': '~', 'step': 'variant-calling', 'tag': 'variants'})
        with open(delivfile, 'w') as out:
            yaml.dump(deliv, out)

    def create_fohm_csv(self):

        """Creates a summary file for FoHM for each region-lab-combination"""

        # Add header to summary files
        for regionlab in self.regionlabs:
            sumfile = os.path.join(self.indir, "{}_{}_komplettering.csv".format(regionlab.replace(' ', '_'), self.today))
            with open(sumfile, "w") as summary:
                summary.write("provnummer,urvalskriterium,GISAID_accession\n")

        # Write sample information to corresponding summary file
        for record in self.caseinfo:
            region = record["region_code"].replace(' ', '_')
            lab = record["lab_code"].replace(' ', '_')
            sumfile = os.path.join(self.indir, "{}_{}_{}_komplettering.csv".format(region, lab, self.today))
            with open(sumfile, "a") as out:
                summary = csv.writer(out)
                summary.writerow([record["Customer_ID_sample"], record["selection_criteria"]] )


