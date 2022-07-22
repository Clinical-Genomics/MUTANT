""" Using a dict as input, this class will print a report covering the
    information requested by the sarscov2-customers at Clinical Genomics
"""
import glob
from datetime import date

from mutant.modules.artic_nanopore.parser import collect_results, collect_variants
from mutant.modules.generic_parser import get_sarscov2_config, read_filelines
from mutant.modules.generic_reporter import GenericReporter


class ReportPrinterNanopore:
    def __init__(self, caseinfo: str, indir: str, barcode_to_sampleid: dict):
        self.casefile = caseinfo
        self.caseinfo = get_sarscov2_config(caseinfo)
        self.case = self.caseinfo[0]["case_ID"]
        self.ticket = self.caseinfo[0]["Customer_ID_project"]
        self.indir = indir
        self.barcode_to_sampleid = barcode_to_sampleid

    def create_all_nanopore_files(self):
        result: dict = collect_results(
            resdir=self.indir,
            barcode_to_sampleid=self.barcode_to_sampleid,
            caseinfo=self.casefile,
        )
        variants: list = collect_variants(
            resdir=self.indir, barcode_to_sampleid=self.barcode_to_sampleid
        )
        self.print_report(result=result)
        self.print_variants(variants=variants)
        generic_reporter = GenericReporter(
                indir=self.indir,
                nanopore=True,
            )
        generic_reporter.create_trailblazer_config()
        self.create_concat_pangolin()
        self.create_concat_pangolin_fohm(result=result)

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
                    "dna_var\n",
                ]
            )
            file_to_append.write(header_results)
            for line in variants:
                file_to_append.write(line + "\n")
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

    def extract_barcode_from_pangolin_csv(self, line: str) -> str:
        """ line in format noblecat_220721-191417_barcode01/ARTIC/medaka_MN908947.3 """
        split_on_slash = line.split("/")
        splid_on_underscore = split_on_slash[0].split("_")
        return splid_on_underscore[2]

    def create_concat_pangolin(self) -> None:
        """Concatenate nanopore pangolin results"""
        indir = "{0}/articNcovNanopore_sequenceAnalysisMedaka_pangolinTyping".format(self.indir)
        concatfile = "{0}/{1}.pangolin.csv".format(self.indir, self.ticket)
        pango_csvs = glob.glob("{0}/*.pangolin.csv".format(indir))

        header = read_filelines(pango_csvs[0])[0]
        with open(concatfile, "w") as concat:
            concat.write(header)
            # Parse sample pangolin data
            for csv in pango_csvs:
                data: list = read_filelines(csv)[1:]
                for line in data:
                    split_on_comma = line.split(",")
                    barcode = self.extract_barcode_from_pangolin_csv(line=split_on_comma[0])
                    split_on_comma[0] = self.barcode_to_sampleid[barcode]
                    concatenated_line = ""
                    for section in split_on_comma:
                        concatenated_line = ",".join([concatenated_line, section])
                    formatted_line = concatenated_line[1:]
                    concat.write(formatted_line)

    def create_concat_pangolin_fohm(self, result: dict) -> None:
        """Prints concatenated pangolin csv for samples passing QC"""
        concatenated_pangolin = "{0}/{1}.pangolin.csv".format(self.indir, self.ticket)
        pango_fohm = "{0}/{1}_{2}_pangolin_classification_format4.txt".format(
            self.indir, self.ticket, str(date.today())
        )
        with open(concatenated_pangolin) as f:
            lines = f.readlines()

        with open(pango_fohm, "w") as passed_qc:
            passed_qc.write(lines[0])
            for line in lines[1:]:
                split_on_comma = line.split(",")
                sample_id = split_on_comma[0]
                if result[sample_id]["qc_pass"] == "TRUE":
                    passed_qc.write(line)
