"""Parse SARS-CoV-2 nanopore results"""
import os

from pathlib import Path
from mutant.modules.generic_parser import get_sarscov2_config


class ParserNanopore:
    def __init__(self, caseinfo: str, indir: Path):
        self.caseinfo = caseinfo
        self.indir = indir

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

#    def parse_assembly(self, results: dict) -> dict:
#        pass

#    def calculate_coverage(self, results: dict) -> dict:
#        pass

    def get_second_line(self, filename: str) -> str:
        """Return the second line of the file"""
        opened_file = open(filename)
        second_line = ""
        for i, line in enumerate(opened_file):
            if i == 1:
                second_line = line
        return second_line

    def get_cust_sample_id(self, raw_pangolin_result: str, barcode_translation: dict) -> str:
        """Return the customer ID of a sample"""
        split_on_slash = raw_pangolin_result.split("/")
        sample_folder = split_on_slash[0]
        split_on_underscore = sample_folder.split("_")
        barcode = split_on_underscore[-1]
        cust_sample_id = barcode_translation[barcode]
        return cust_sample_id

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
            cust_sample_id: str = self.get_cust_sample_id(raw_pangolin_result=second_line, barcode_translation=barcode_to_sample)
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
        #results: dict = self.parse_assembly(results=results)
        #results: dict = self.calculate_coverage(results=results)
        results: dict = self.parse_pangolin(results=results, barcode_to_sample=barcode_to_sample, resdir=resdir)
        #results: dict = self.parse_classifications(results=results)
        #results: dict = self.parse_mutations(results=results)
        return results
