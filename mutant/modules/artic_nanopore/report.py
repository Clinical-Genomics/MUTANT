""" Using a dict as input, this class will print a report covering the
    information requested by the sarscov2-customers at Clinical Genomics
"""
from mutant.modules.generic_parser import (
    get_sarscov2_config,
)

class ReportPrinterNanopore:
    def __init__(self, indir, caseinfo):
        self.casefile = caseinfo
        caseinfo = get_sarscov2_config(caseinfo)
        self.case = caseinfo[0]["case_ID"]
        self.ticket = caseinfo[0]["Customer_ID_project"]
        self.project = caseinfo[0]["Customer_ID_project"]
        self.indir = indir

    def print_report(self, result: dict) -> None:
        pass
