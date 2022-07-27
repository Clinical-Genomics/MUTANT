""" This class modifies existing files to fit the Clinical Genomics infrastructure.
    Specifically it acts on the sarscov2 pipeline.

    By: Isak Sylvin & Tanja Normark

"""

import glob
import os

from mutant.modules.generic_parser import get_sarscov2_config


class DeliverySC2:
    def __init__(self, caseinfo: str, indir: str, nanopore: bool, timestamp: str, barcode_to_sampleid: dict):
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
        illunmina_to_nanopore_base = {}
        for sampleinfo in self.caseinfo:
            illumina_base = "{0}_{1}_{2}".format(
                sampleinfo["region_code"], sampleinfo["lab_code"], sampleinfo["Customer_ID_sample"]
            )
            nanopore_base = "{0}_{1}_{2}".format(
                self.case, self.timestamp, sampleid_to_barcode[sampleinfo["Customer_ID_sample"]]
            )
            illunmina_to_nanopore_base[illumina_base] = nanopore_base
        return illunmina_to_nanopore_base

    def rename_deliverables(self):
        """Rename result files for delivery: fastq, consensus files, vcf and pangolin"""

        sampleid_to_barcode = {}
        illunmina_to_nanopore_base = {}
        if self.nanopore:
            sampleid_to_barcode = self.map_sampleid_to_barcode()
            illunmina_to_nanopore_base: dict = self.translate_basename_to_nanopore_format(sampleid_to_barcode=sampleid_to_barcode)

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
                consensus_path = "{0}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka".format(
                    self.indir
                )
                consensus_file = "{0}/{1}.consensus.fasta".format(consensus_path, illunmina_to_nanopore_base[base_sample])
            else:
                consensus_path = "{0}/ncovIllumina_sequenceAnalysis_makeConsensus".format(
                    self.indir
                )
                consensus_file = "{0}/{1}.*.fa".format(consensus_path, base_sample)
            for item in glob.glob(consensus_files):
                newpath = "{0}/{1}.consensus.fasta".format(consensus_path, base_sample)
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
                vcf_path = "{0}/articNcovNanopore_Genotyping_typeVariants/vcf".format(self.indir)
                vcf_files = "{0}/*.csq.vcf".format(vcf_path)
            else:
                vcf_path = "{0}/ncovIllumina_Genotyping_typeVariants/vcf".format(self.indir)
                vcf_files = "{0}/{1}.csq.vcf".format(vcf_path, base_sample)
            for item in glob.glob(vcf_files):
                newpath = "{0}/{1}.vcf".format(vcf_path, base_sample)
                try:
                    os.symlink(item, newpath)
                except Exception as e:
                    pass

        ## Rename case files

        # rename multiqc
        hit = glob.glob(
            "{}/QCStats/ncovIllumina_sequenceAnalysis_multiqc/*_multiqc.html".format(
                self.indir
            )
        )
        if len(hit) == 1:
            hit = hit[0]
            try:
                os.symlink(hit, "{}/{}_multiqc.html".format(self.indir, self.ticket))
            except Exception as e:
                pass

        # rename multiqc json
        hit = glob.glob(
            "{}/QCStats/ncovIllumina_sequenceAnalysis_multiqc/*_multiqc_data/multiqc_data.json".format(
                self.indir
            )
        )
        if len(hit) == 1:
            hit = hit[0]
            try:
                os.symlink(hit, "{}/{}_multiqc.json".format(self.indir, self.ticket))
            except Exception as e:
                pass

        core_suffix = [
            ".qc.csv",
            ".typing_summary.csv",
            ".variant_summary.csv",
        ]
        for thing in core_suffix:
            hit = glob.glob("{0}/*{1}".format(self.indir, thing))
            if len(hit) == 1:
                hit = hit[0]
                try:
                    os.symlink(hit, "{0}/{1}{2}".format(self.indir, self.ticket, thing))
                except Exception as e:
                    pass
