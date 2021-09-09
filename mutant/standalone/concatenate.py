""" This script concatenates all forward fastq files into one file and
    all reverse fastq files into one file for each sample. The files
    that are used for concatenation will be removed afterwards if they
    got more than 1 inode.
    By: @henningonsbring """

import os
import re
import subprocess
import sys

from argparse import ArgumentParser
from pathlib import Path


PREFIX_TO_CONCATENATE = ["MWG", "MWL", "MWM", "MWR", "MWX", "CONCATENATE"]
parser = ArgumentParser()

parser.add_argument(
    "-i",
    "--input_folder",
    dest="input_folder",
    help="Folder with fastq to concatenate",
    metavar="<PATH>",
    required=True,
    type=str,
)
parser.add_argument(
    "-a",
    "--app_tag",
    dest="app_tag",
    help="Application tag",
    metavar="<STRING>",
    required=False,
    type=str,
    default="CONCATENATE",
)
parser.add_argument(
    "-d",
    "--date",
    dest="date",
    help="Date to add to the concatenated file name, e.g. order date",
    metavar="<DATE>",
    required=False,
    type=str,
    default="",
)

args = parser.parse_args()
should_concatenate = False

for prefix in PREFIX_TO_CONCATENATE:
    if args.app_tag.startswith(prefix):
        print(
            "Apptag %s identified, data generated with this application tag should be concatenated"
            % (args.app_tag)
        )
        should_concatenate = True

if should_concatenate == False:
    print(
        "Data with application tag %s should not be concatenated, skipping concatenation"
        % (args.app_tag)
    )
    sys.exit(-1)

for dir_name in os.listdir(args.input_folder):
    dir_path = os.path.join(args.input_folder, dir_name)
    if len(os.listdir(dir_path)) == 0:
        print("Empty folder found: %s" % (dir_path))
        cmd = "rmdir " + str(dir_path)
        os.system(cmd)
        print("Running command: %s" % (cmd))
        continue
    if not os.path.isdir(dir_path):
        continue
    for read_direction in [1, 2]:
        same_direction = []
        total_size = 0
        for file in os.listdir(dir_path):
            abs_path_file = os.path.join(dir_path, file)
            direction_string = ".+_R" + str(read_direction) + "_[0-9]+.fastq.gz"
            direction_pattern = re.compile(direction_string)
            if direction_pattern.match(abs_path_file):
                same_direction.append(os.path.abspath(abs_path_file))
                total_size = total_size + Path(abs_path_file).stat().st_size
        same_direction.sort()
        cmd = "cat"
        for i in range(len(same_direction)):
            cmd = cmd + " " + same_direction[i]
        if args.date:
            output = (
                dir_path
                + "/"
                + str(args.date)
                + "_"
                + dir_name
                + "_"
                + str(read_direction)
                + ".fastq.gz"
            )
        else:
            output = dir_path + "/" + dir_name + "_" + str(read_direction) + ".fastq.gz"
        cmd = cmd + " > " + output
        print("Running command: %s" % (cmd))
        os.system(cmd)
        concatenated_size = Path(output).stat().st_size
        if total_size == concatenated_size:
            print(
                "QC PASSED: Total size for files used in concatenation match the size of the concatenated file"
            )
            for file in same_direction:
                inode_check_cmd = "stat -c %h " + file
                n_inodes = subprocess.getoutput(inode_check_cmd)
                if int(n_inodes) > 1:
                    print("Removing file: %s" % (file))
                    os.remove(file)
                else:
                    print(
                        "WARNING %s only got 1 inode, file will not be removed" % (file)
                    )
        else:
            print("WARNING data lost in concatenation")
