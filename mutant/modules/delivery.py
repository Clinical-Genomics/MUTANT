""" This class modifies existing files to fit the Clinical Genomics infrastructure.
    Specifically it acts on the sarscov2 pipeline.
"""

import glob
import os

from mutant.modules.generic_parser import get_sarscov2_config


class DeliverySC2:
    def __init__(
        self,
        caseinfo: str,
        indir: str,
        nanopore: bool,
        timestamp: str,
        barcode_to_sampleid: dict,
    ):
        self.casefile = caseinfo
        caseinfo = get_sarscov2_config(caseinfo)
        self.caseinfo = caseinfo
        self.case = caseinfo[0]["case_ID"]
        self.ticket = caseinfo[0]["Customer_ID_project"]
        self.indir = indir
        self.nanopore = nanopore
        self.timestamp = timestamp
        self.barcode_to_sampleid = barcode_to_sampleid

    def map_sampleid_to_barcode(self) -> dict:
        sampleid_to_barcode = {}
        for barcode in self.barcode_to_sampleid:
            sampleid = self.barcode_to_sampleid[barcode]
            sampleid_to_barcode[sampleid] = barcode
        return sampleid_to_barcode

    def translate_basename_to_nanopore_format(self, sampleid_to_barcode: dict) -> dict:
        illumina_to_nanopore_base = {}
        for sampleinfo in self.caseinfo:
            illumina_base = "{0}_{1}_{2}".format(
                sampleinfo["region_code"],
                sampleinfo["lab_code"],
                sampleinfo["Customer_ID_sample"],
            )
            nanopore_base = "{0}_{1}_{2}".format(
                self.case,
                self.timestamp,
                sampleid_to_barcode[sampleinfo["Customer_ID_sample"]],
            )
            illumina_to_nanopore_base[illumina_base] = nanopore_base
        return illumina_to_nanopore_base

    def rename_deliverables(self):
        """Rename result files for delivery: fastq, consensus files, vcf and pangolin"""

        illumina_to_nanopore_base = {}
        if self.nanopore:
            sampleid_to_barcode = self.map_sampleid_to_barcode()
            illumina_to_nanopore_base: dict = (
                self.translate_basename_to_nanopore_format(
                    sampleid_to_barcode=sampleid_to_barcode
                )
            )

        # Rename sample folders
        if self.nanopore:
            fastq_folder = "{0}/../fastq".format(self.indir)
            for barcode in os.listdir(fastq_folder):
                barcode_abspath = fastq_folder + "/" + barcode
                sampleid_abspath = (
                    fastq_folder + "/" + self.barcode_to_sampleid[barcode]
                )
                os.rename(barcode_abspath, sampleid_abspath)

        # Rename sample files
        for sampleinfo in self.caseinfo:
            region = sampleinfo["region_code"]
            lab = sampleinfo["lab_code"]
            base_sample = "{0}_{1}_{2}".format(
                region, lab, sampleinfo["Customer_ID_sample"]
            )
            if not self.nanopore:
                if not sampleinfo["sequencing_qc_pass"]:
                    continue

            # rename makeConsensus
            if self.nanopore:
                consensus_path = f"{self.indir}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka"
                consensus_files = f"{consensus_path}/{illumina_to_nanopore_base[base_sample]}.consensus.fasta"
            else:
                consensus_path = f"{self.indir}/ncovIllumina_sequenceAnalysis_makeConsensus"
                consensus_files = f"{consensus_path}/{base_sample}.*.fa"
            for item in glob.glob(consensus_files):
                newpath = f"{consensus_path}/{base_sample}.consensus.fasta"
                try:
                    with open(item, "r") as old_consensus_io:
                        with open(newpath, "w") as new_consensus_io:
                            new_consensus_io.write(
                                f">{base_sample}\n{old_consensus_io.readlines()[1]}"
                            )

                except Exception as e:
                    print(e)
                    pass

            # rename typeVariants
            if self.nanopore:
                vcf_path = f"{self.indir}/articNcovNanopore_Genotyping_typeVariants/vcf"
                vcf_files = f"{vcf_path}/{illumina_to_nanopore_base[base_sample]}.csq.vcf"
            else:
                vcf_path = f"{self.indir}/ncovIllumina_Genotyping_typeVariants/vcf"
                vcf_files = f"{vcf_path}/{base_sample}.csq.vcf"
            for item in glob.glob(vcf_files):
                newpath = f"{vcf_path}/{base_sample}.vcf"
                try:
                    os.symlink(item, newpath)
                except Exception as e:
                    pass

        ## Rename case files

        # rename multiqc
        if self.nanopore:
            hit = glob.glob(
                f"{self.indir}/QCStats/articNcovNanopore_sequenceAnalysisMedaka_multiqcNanopore/*_multiqc.html"
            )
        else:
            hit = glob.glob(
                f"{self.indir}/QCStats/ncovIllumina_sequenceAnalysis_multiqc/*_multiqc.html"
            )
        if len(hit) == 1:
            hit = hit[0]
            try:
                os.symlink(hit, f"{self.indir}/{self.ticket}_multiqc.html")
            except Exception as e:
                pass

        # rename multiqc json
        if self.nanopore:
            hit = glob.glob(
                f"{self.indir}/QCStats/articNcovNanopore_sequenceAnalysisMedaka_multiqcNanopore/*_multiqc_data/multiqc_data.json"
            )
        else:
            hit = glob.glob(
                f"{self.indir}/QCStats/ncovIllumina_sequenceAnalysis_multiqc/*_multiqc_data/multiqc_data.json"
            )
        if len(hit) == 1:
            hit = hit[0]
            try:
                os.symlink(hit, f"{self.indir}/{self.ticket}_multiqc.json")
            except Exception as e:
                pass

        core_suffix = [
            ".qc.csv",
            ".typing_summary.csv",
            ".variant_summary.csv",
        ]
        for thing in core_suffix:
            hit = glob.glob(f"{self.indir}/*{thing}")
            if len(hit) == 1:
                hit = hit[0]
                try:
                    os.symlink(hit, f"{self.indir}/{self.ticket}{thing}")
                except Exception as e:
                    pass
