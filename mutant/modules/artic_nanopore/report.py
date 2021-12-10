""" Using a dict as input, this class will print a report covering the
    information requested by the sarscov2-customers at Clinical Genomics
"""
from mutant.modules.generic_parser import (
    get_sarscov2_config,
)


class ReportSC2:
    def __init__(self, caseinfo, indir):
        caseinfo = get_sarscov2_config(caseinfo)
        self.caseinfo = caseinfo
        self.indir = indir



