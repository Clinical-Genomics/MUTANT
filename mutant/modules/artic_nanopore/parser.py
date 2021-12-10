"""Parse SARS-CoV-2 nanopore results"""
from mutant.modules.generic_parser import (
    get_sarscov2_config,
)


class ParserNanopore:
    def __init__(self, caseinfo, indir):
        caseinfo = get_sarscov2_config(caseinfo)
        self.caseinfo = caseinfo
        self.indir = indir

    def collect_results(self) -> dict:
        results: dict = self.parse_config()
        results: dict = self.parse_assembly(results=results)
        results: dict = self.calculate_coverage(results=results)
        results: dict = self.parse_pangolin(results=results)
        results: dict = self.parse_classifications(results=results)
        results: dict = self.parse_mutations(results=results)
        return results

    def parse_config(self) -> dict:
        pass

    def parse_assembly(self, results: dict) -> dict:
        pass

    def calculate_coverage(self, results: dict) -> dict:
        pass

    def parse_pangolin(self, results: dict) -> dict:
        pass

    def parse_classifications(self, results: dict) -> dict:
        pass

    def parse_mutations(self, results: dict) -> dict:
        pass
