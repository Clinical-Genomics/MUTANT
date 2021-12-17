""" Using a dict as input, this class will print a report covering the
    information requested by the sarscov2-customers at Clinical Genomics
"""

from mutant.modules.generic_parser import get_sarscov2_config

class ReportPrinterNanopore:
    def __init__(self, caseinfo: str, indir: str):
        self.casefile = caseinfo
        caseinfo = get_sarscov2_config(caseinfo)
        self.case = caseinfo[0]["case_ID"]
        self.ticket = caseinfo[0]["Customer_ID_project"]
        self.indir = indir

    def print_report(self, result: dict) -> None:
        result_file = "/".join([self.indir, "sars-cov-2_results.csv"])
        with open(result_file, "a") as file_to_append:
            file_to_append.write("Sample,Lineage")
            samples = result.keys()
            for sample in samples:
                line_to_append = "{0},{1}".format(sample, result[sample]["pangolin_type"])
                file_to_append.write(line_to_append)
