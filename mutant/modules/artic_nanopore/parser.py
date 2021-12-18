"""Parse SARS-CoV-2 nanopore results"""
import os

from mutant.modules.generic_parser import get_sarscov2_config


class ParserNanopore:
    def __init__(self, caseinfo: str):
        self.caseinfo = caseinfo

#    def get_line(self, filename: str, line_index_of_interest: int) -> str:
#        """Return a certain line of a given file"""
#        opened_file = open(filename)
#        for i, line in enumerate(opened_file):
#            if i == line_index_of_interest:
#                return line

    def get_second_line(self, filename: str) -> str:
        """Return the second line of the file"""
        opened_file = open(filename)
        second_line = ""
        for i, line in enumerate(opened_file):
            if i == 1:
                second_line = line
        return second_line

    def get_first_line(self, filename: str) -> str:
        """Return the first line of the file"""
        opened_file = open(filename)
        first_line = ""
        for i, line in enumerate(opened_file):
            if i == 0:
                first_line = line
        return first_line

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
                total_bases =+ len(stripped_line)
                total_N =+ stripped_line.count('N')
            fraction_N_two_decimals = float("{:.2f}".format(total_N/total_bases))
            return fraction_N_two_decimals

    def parse_assembly(self, results: dict, resdir: str, barcode_to_sample: dict) -> dict:
        """Collects data by parsing the assembly"""
        base_path = "/".join([resdir, "articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka"])
        for filename in os.listdir(base_path):
            abs_path = os.path.join(base_path, filename)
            first_line: str = self.get_first_line(filename=abs_path)
            cust_sample_id: str = self.get_cust_sample_id(line_to_parse=first_line, barcode_translation=barcode_to_sample)
            fraction_N: float = self.get_fraction_n(input_file=abs_path)
            results[cust_sample_id]["fraction_n_bases"] = fraction_N
        return results

#    def calculate_coverage(self, results: dict) -> dict:
#        pass

    def get_pangolin_type(self, raw_pangolin_result: str) -> str:
        """Return the pangolin type of a sample"""
        split_on_comma = raw_pangolin_result.split(",")
        lineage = split_on_comma[1]
        return lineage

    def parse_pangolin(self, results: dict, barcode_to_sample: dict, resdir: str) -> dict:
        """Collect data for pangolin types"""
        base_path = "/".join([resdir, "articNcovNanopore_sequenceAnalysisMedaka_pangolinTyping"])
        for filename in os.listdir(base_path):
            abs_path = os.path.join(base_path, filename)
            second_line: str = self.get_second_line(filename=abs_path)
            cust_sample_id: str = self.get_cust_sample_id(line_to_parse=second_line, barcode_translation=barcode_to_sample)
            pangolin_type: str = self.get_pangolin_type(raw_pangolin_result=second_line)
            results[cust_sample_id]["pangolin_type"] = pangolin_type
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
        #results: dict = self.calculate_coverage(results=results)
        results: dict = self.parse_pangolin(results=results, barcode_to_sample=barcode_to_sample, resdir=resdir)
        #results: dict = self.parse_classifications(results=results)
        #results: dict = self.parse_mutations(results=results)
        return results
