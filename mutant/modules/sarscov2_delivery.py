""" This class modifies existing files to fit the Clinical Genomics infrastructure.
    Specifically it acts on the sarscov2 pipeline.

    By: Isak Sylvin & Tanja Normark

"""

import glob
import os

from mutant.modules.generic_parser import get_sarscov2_config


class DeliverySC2:
    def __init__(self, caseinfo, indir):
        self.casefile = caseinfo
        caseinfo = get_sarscov2_config(caseinfo)
        self.caseinfo = caseinfo
        self.case = caseinfo[0]["case_ID"]
        self.ticket = caseinfo[0]["Customer_ID_project"]
        self.project = caseinfo[0]["Customer_ID_project"]
        self.indir = indir



    def rename_deliverables(self):
        """Rename result files for delivery: fastq, consensus files, vcf and pangolin"""

        # Rename sample files

        for sampleinfo in self.caseinfo:
            sample = sampleinfo["CG_ID_sample"]
            region = sampleinfo["region_code"]
            lab = sampleinfo["lab_code"]
            base_sample = "{0}_{1}_{2}".format(
                region, lab, sampleinfo["Customer_ID_sample"]
            )
            if not sampleinfo["sequencing_qc_pass"]:
                continue

            # rename makeConsensus
            prefix = "{0}/ncovIllumina_sequenceAnalysis_makeConsensus".format(
                self.indir
            )
            for item in glob.glob("{0}/{1}.*".format(prefix, base_sample)):
                newpath = "{0}/{1}.consensus.fasta".format(prefix, base_sample)
                try:
                    with open(item, "r") as old_consensus_io:
                        with open(newpath, "w") as new_consensus_io:
                            new_consensus_io.write(f">{base_sample}\n")
                            new_consensus_io.write(old_consensus_io.readlines()[1])

                except Exception as e:
                    pass

            # rename typeVariants
            prefix = "{0}/ncovIllumina_Genotyping_typeVariants/vcf".format(self.indir)
            for item in glob.glob("{0}/{1}.csq.vcf".format(prefix, base_sample)):
                newpath = "{0}/{1}.vcf".format(prefix, base_sample)
                try:
                    os.symlink(item, newpath)
                except Exception as e:
                    pass

        ## Rename case files

        # rename multiqc
        hit = glob.glob("{}/multiqc/*_multiqc.html".format(self.indir))
        if len(hit) == 1:
            hit = hit[0]
            try:
                os.symlink(hit, "{}/{}_multiqc.html".format(self.indir, self.ticket))
            except Exception as e:
                pass

        # rename multiqc json
        hit = glob.glob(
            "{}/multiqc/*_multiqc_data/multiqc_data.json".format(self.indir)
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
