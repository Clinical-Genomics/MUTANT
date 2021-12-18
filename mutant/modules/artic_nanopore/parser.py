"""Parse SARS-CoV-2 nanopore results"""
import glob
import os

from mutant.modules.generic_parser import get_sarscov2_config

QC_PASS_THRESHOLD_COVERAGE_10X_OR_HIGHER = 95

class ParserNanopore:
    def __init__(self, caseinfo: str):
        self.caseinfo = caseinfo

    def get_line(self, filename: str, line_index_of_interest: int) -> str:
        """Return a certain line of a given file"""
        opened_file = open(filename)
        for i, line in enumerate(opened_file):
            if i == line_index_of_interest:
                return line

    def get_cust_sample_id(self, line_to_parse: str, barcode_translation: dict) -> str:
        """Return the customer ID of a sample"""
        split_on_slash = line_to_parse.split("/")
        sample_folder = split_on_slash[0]
        split_on_underscore = sample_folder.split("_")
        barcode = split_on_underscore[-1]
        cust_sample_id = barcode_translation[barcode]
        return cust_sample_id

    def translate_barcodes(self, parsed_config: dict) -> dict:
        """Builds a dict where barcodes point at customer sample id"""
        barcode_to_sample = {}
        for sample in parsed_config:
            barcode_to_sample[sample["barcode"]] = sample["Customer_ID_sample"]
        return barcode_to_sample

    def get_data_from_config(self, parsed_config: dict) -> dict:
        """Collect data for selection criteria and region"""
        data_to_report = {}
        for sample in parsed_config:
            cust_sample_id = sample["Customer_ID_sample"]
            data_to_report[cust_sample_id] = {}
            data_to_report[cust_sample_id]["selection_criteria"] = sample["selection_criteria"]
            data_to_report[cust_sample_id]["region_code"] = sample["region_code"]
        return data_to_report

    def get_fraction_n(self, input_file: str) -> float:
        """Calculates fraction of N bases in a fasta file"""
        with open(input_file, "r") as fasta:
            total_bases = 0
            total_N = 0
            for line in fasta:
                stripped_line = line.strip()
                total_bases += len(stripped_line)
                total_N += stripped_line.count('N')
            percentage_N_two_decimals = round((total_N/total_bases)*100, 2)
            return percentage_N_two_decimals

    def parse_assembly(self, results: dict, resdir: str, barcode_to_sample: dict) -> dict:
        """Collects data by parsing the assembly"""
        base_path = "/".join([resdir, "articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka"])
        for filename in os.listdir(base_path):
            if filename.endswith(".consensus.fasta"):
                abs_path = os.path.join(base_path, filename)
                first_line: str = self.get_line(filename=abs_path, line_index_of_interest=0)
                cust_sample_id: str = self.get_cust_sample_id(line_to_parse=first_line, barcode_translation=barcode_to_sample)
                fraction_N: float = self.get_fraction_n(input_file=abs_path)
                results[cust_sample_id]["fraction_n_bases"] = fraction_N
        return results

    def get_depth_files_paths(self, barcode: str, base_path: str) -> list:
        """Returns the paths to coverage statistics for a certain barcode"""
        wild_card_path = "*".join([base_path, barcode, "depths"])
        depth_files = glob.glob(wild_card_path)
        return depth_files

    def count_bases_w_10x_cov_or_more(self, coverage_stats: list) -> int:
        """Counts how many times a number is higher than 10 in a list"""
        numbers_higher_than_10 = 0
        for number in coverage_stats:
            if number >= 10:
                numbers_higher_than_10 += 1
        return numbers_higher_than_10

    def initiate_coverage_stats_list(self, file_path: str) -> list:
        coverage_stats = []
        with open(file_path, "r") as file1:
            for line in file1:
                stripped_line = line.strip()
                columns: list = stripped_line.split("\t")
                coverage_stats.append(int(columns[3]))
        file1.close()
        return coverage_stats

    def calculate_coverage(self, results: dict, resdir: str, barcode_to_sample: dict) -> dict:
        """Collects data for the fraction of each assembly that got 10x coverage or more"""
        base_path = "/".join([resdir, "articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka/"])
        for barcode in barcode_to_sample:
            coverage_files_paths: list = self.get_depth_files_paths(barcode=barcode, base_path=base_path)
            coverage_stats: list = self.initiate_coverage_stats_list(file_path=coverage_files_paths[0])
            with open(coverage_files_paths[1], "r") as file2:
                for line in file2:
                    stripped_line = line.strip()
                    columns: list = stripped_line.split("\t")
                    coverage_stats[int(columns[2])] += int(columns[3])
            file2.close()
            bases_w_10x_cov_or_more: int = self.count_bases_w_10x_cov_or_more(coverage_stats=coverage_stats)
            percentage_equal_or_greater_than_10 = round((bases_w_10x_cov_or_more / len(coverage_stats)) * 100, 2)
            results[barcode_to_sample[barcode]]["pct_10x_coverage"] = percentage_equal_or_greater_than_10
            if percentage_equal_or_greater_than_10 >= QC_PASS_THRESHOLD_COVERAGE_10X_OR_HIGHER:
                results[barcode_to_sample[barcode]]["qc_pass"] = "TRUE"
            else:
                results[barcode_to_sample[barcode]]["qc_pass"] = "FALSE"
        return results

    def get_pangolin_type(self, raw_pangolin_result: str) -> str:
        """Return the pangolin type of a sample"""
        split_on_comma = raw_pangolin_result.split(",")
        lineage = split_on_comma[1]
        return lineage

    def get_pangoLEARN_version(self, raw_pangolin_result: str) -> str:
        """Return the pangoLEARN_version of a sample"""
        split_on_comma = raw_pangolin_result.split(",")
        pango_learn_version = split_on_comma[9]
        return pango_learn_version

    def parse_pangolin(self, results: dict, barcode_to_sample: dict, resdir: str) -> dict:
        """Collect data for pangolin types"""
        base_path = "/".join([resdir, "articNcovNanopore_sequenceAnalysisMedaka_pangolinTyping"])
        for filename in os.listdir(base_path):
            abs_path = os.path.join(base_path, filename)
            second_line: str = self.get_line(filename=abs_path, line_index_of_interest=1)
            cust_sample_id: str = self.get_cust_sample_id(line_to_parse=second_line, barcode_translation=barcode_to_sample)
            pangolin_type: str = self.get_pangolin_type(raw_pangolin_result=second_line)
            results[cust_sample_id]["pangolin_type"] = pangolin_type
            pangoLEARN_version: str = self.get_pangoLEARN_version(raw_pangolin_result=second_line)
            results[cust_sample_id]["pangolearn_version"] = pangoLEARN_version
        return results

#    def parse_classifications(self, results: dict) -> dict:
#        pass

#    def parse_mutations(self, results: dict) -> dict:
#        pass

    def collect_results(self, resdir: str) -> dict:
        """Build a dictionary with data for the report"""
        parsed_config: dict = get_sarscov2_config(config=self.caseinfo)
        barcode_to_sample: dict = self.translate_barcodes(parsed_config=parsed_config)
        results: dict = self.get_data_from_config(parsed_config=parsed_config)
        results: dict = self.parse_assembly(results=results, resdir=resdir, barcode_to_sample=barcode_to_sample)
        results: dict = self.calculate_coverage(results=results, resdir=resdir, barcode_to_sample=barcode_to_sample)
        results: dict = self.parse_pangolin(results=results, barcode_to_sample=barcode_to_sample, resdir=resdir)
        #results: dict = self.parse_classifications(results=results)
        #results: dict = self.parse_mutations(results=results)
        return results
