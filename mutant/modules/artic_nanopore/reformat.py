"""Reformat SARS-CoV-2 nanopore data organisation"""
import gzip
import os


def get_barcode_from_custid(sampleid_abspath: str) -> str:
    fastq_files: list = os.listdir(sampleid_abspath)
    with gzip.open(fastq_files[0], "rt") as f:
        for line in f:
            nanopore_run_info = str(line).split(" ")
            barcode_info = nanopore_run_info[8].split("=")
            return barcode_info[1]


def reformat_fastq_folder(fastq_folder: str) -> dict:
    """Changes the name of folders with custid to barcode and
    returns a dict where barcode points on custid"""
    barcode_to_sampleid = {}
    external_sampleids: list = os.listdir(fastq_folder)
    for sampleid in external_sampleids:
        sampleid_abspath = fastq_folder + "/" + sampleid
        barcode: str = get_barcode_from_custid(sampleid_abspath=sampleid_abspath)
        barcode_abspath = fastq_folder + "/" + barcode
        os.rename(sampleid_abspath, barcode_abspath)
        barcode_to_sampleid[barcode] = sampleid
    return barcode_to_sampleid
