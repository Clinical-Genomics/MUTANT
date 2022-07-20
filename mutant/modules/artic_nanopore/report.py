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

    def create_all_nanopore_files(self, result: dict, variants: list):
        self.print_report(result=result)
        self.print_variants(variants=variants)

    def print_variants(self, variants: list) -> None:
        """Append data to the variant report"""
        file_name_report = "_".join(["sars-cov-2", str(self.ticket), "variants.csv"])
        variants_file = "/".join([self.indir, file_name_report])
        with open(variants_file, "a") as file_to_append:
            header_results = ",".join(
                [
                    "sampleID",
                    "gene",
                    "aa_var",
                    "dna_var",
                ]
            )
            file_to_append.write(header_results)
            for line in variants:
                file_to_append.write(line)
        file_to_append.close()

    def print_report(self, result: dict) -> None:
        """Append results from the analysis to a report"""
        file_name_report = "_".join(["sars-cov-2", str(self.ticket), "results.csv"])
        result_file = "/".join([self.indir, file_name_report])
        with open(result_file, "a") as file_to_append:
            header_results = ",".join(
                [
                    "Sample",
                    "Selection",
                    "Region Code",
                    "Ticket",
                    "%N_bases",
                    "%10X_coverage",
                    "QC_pass",
                    "Lineage",
                    "PangoLEARN_version",
                    "VOC",
                    "Mutations\n",
                ]
            )
            file_to_append.write(header_results)
            samples = result.keys()
            for sample in samples:
                line_to_append = (
                    "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10}{11}".format(
                        sample,
                        result[sample]["selection_criteria"],
                        result[sample]["region_code"],
                        self.ticket,
                        result[sample]["fraction_n_bases"],
                        result[sample]["pct_10x_coverage"],
                        result[sample]["qc_pass"],
                        result[sample]["pangolin_type"],
                        result[sample]["pangolearn_version"],
                        result[sample]["voc"],
                        result[sample]["mutations"],
                        "\n",
                    )
                )
                file_to_append.write(line_to_append)
        file_to_append.close()
