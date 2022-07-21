from pathlib import Path

import yaml


class GenericReporter:
    def __init__(self, indir: str, nanopore: bool):
        self.indir = indir
        self.nanopore = nanopore


    def get_finished_slurm_ids(self) -> list:
        """Get slurm IDs"""

        trace_file_path = Path(self.indir, "pipeline_info", "execution_trace.txt")
        slurm_id_list = []
        with open(trace_file_path, "r") as trace_file_contents:
            for line in trace_file_contents:
                slurm_id = line.split()[2]
                try:
                    slurm_id_list.append(int(slurm_id))
                except Exception:
                    continue
        return slurm_id_list

    def create_trailblazer_config(self) -> None:
        """Create Trailblazer config file"""

        trailblazer_config_path = Path(self.indir, "trailblazer_config.yaml")
        finished_slurm_ids = self.get_finished_slurm_ids()
        if not finished_slurm_ids:
            return
        with open(trailblazer_config_path, "w") as trailblazer_config_file:
            yaml.dump(data={"jobs": finished_slurm_ids}, stream=trailblazer_config_file)

    def create_concat_pangolin(self):
        """Concatenate pangolin results"""

        indir = "{0}/ncovIllumina_sequenceAnalysis_pangolinTyping".format(self.indir)
        concatfile = "{0}/{1}.pangolin.csv".format(self.indir, self.ticket)
        pangolins = glob.glob("{0}/*.pangolin.csv".format(indir))
        # Copy header
        header = read_filelines(pangolins[0])[0]
        with open(concatfile, "w") as concat:
            concat.write(header)
            # Parse sample pangolin data
            for pango in pangolins:
                data = read_filelines(pango)[1:]
                for line in data:
                    # Use sample name at taxon field
                    taxon_regex = "(\w+)_(\w+)_(\w+)_(?P<name>\w+).(\S+)"
                    sample, subs = re.subn(taxon_regex, r"\g<name>", line.split(",")[0])
                    if subs == 0:
                        print(
                            "Unable to rename taxon - using original: {}".format(pango)
                        )
                    else:
                        line = sample + line[line.find(","):]
                    concat.write(line)
