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
