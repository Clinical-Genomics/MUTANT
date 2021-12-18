""" Using a dict as input, this class will print a report covering the
    information requested by the sarscov2-customers at Clinical Genomics
"""

from mutant.modules.generic_parser import get_sarscov2_config

class ReportPrinterNanopore:
    def __init__(self, caseinfo: str, indir: str):
        self.casefile = caseinfo
        self.caseinfo = get_sarscov2_config(caseinfo)
        self.case = self.caseinfo[0]["case_ID"]
        self.ticket = self.caseinfo[0]["Customer_ID_project"]
        self.indir = indir

    def print_report(self, result: dict) -> None:
        """Append results from the analysis to a report"""
        file_name_report = "_".join(["sars-cov-2", self.ticket, "results.csv"])
        result_file = "/".join([self.indir, file_name_report])
        with open(result_file, "a") as file_to_append:
            file_to_append.write("Sample,Lineage\n")
            samples = result.keys()
            for sample in samples:
                line_to_append = "{0},{1}{2}".format(sample, result[sample]["pangolin_type"], "\n")
                file_to_append.write(line_to_append)
