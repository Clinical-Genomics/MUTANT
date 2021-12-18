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
        file_name_report = "_".join(["sars-cov-2", str(self.ticket), "results.csv"])
        result_file = "/".join([self.indir, file_name_report])
        with open(result_file, "a") as file_to_append:
            file_to_append.write("Sample,Selection,Region Code,Ticket,%N_bases,%10X_coverage,QC_pass,Lineage,PangoLEARN_version\n")
            samples = result.keys()
            for sample in samples:
                line_to_append = "{0},{1},{2},{3},{4},{5},{6},{7},{8}{9}".format(
                    sample,
                    result[sample]["selection_criteria"],
                    result[sample]["region_code"],
                    self.ticket,
                    result[sample]["fraction_n_bases"],
                    result[sample]["pct_10x_coverage"],
                    result[sample]["qc_pass"],
                    result[sample]["pangolin_type"],
                    result[sample]["pangolearn_version"],
                    "\n"
                )
                file_to_append.write(line_to_append)
        file_to_append.close()
