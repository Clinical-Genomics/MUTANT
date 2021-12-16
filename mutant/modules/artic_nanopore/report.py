""" Using a dict as input, this class will print a report covering the
    information requested by the sarscov2-customers at Clinical Genomics
"""
from collections import abc
from pathlib import Path

from mutant.modules.generic_parser import get_sarscov2_config

class ReportPrinterNanopore:
    def __init__(self, caseinfo: str, indir: Path):
        self.casefile = caseinfo
        caseinfo = get_sarscov2_config(caseinfo)
        self.case = caseinfo[0]["case_ID"]
        self.ticket = caseinfo[0]["Customer_ID_project"]
        self.indir = indir

    def print_report(self, result: dict) -> None:
        samples = result.keys()
        for sample in samples:
            analysis_data = result[sample].keys()
            for data_point in analysis_data:
                message = "{0} for sample {1} has value {2}".format(data_point, sample, result[sample][data_point])
                print(message)
