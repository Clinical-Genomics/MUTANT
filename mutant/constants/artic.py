"""Constants for GMS-Artic and MUTANT"""

# GMS-Artic and MUTANT output files
ARTIC_FILES_CASE = {
    "multiqc-json": "{resdir}/QCStats/ncovIllumina_sequenceAnalysis_multiqc/*_multiqc_data/multiqc_data.json",
    "multiqc-html": "{resdir}/QCStats/ncovIllumina_sequenceAnalysis_multiqc/*_multiqc.html",
    "vogue-metrics": "{resdir}/{case}_metrics_deliverables.yaml",
    "results-file": "{resdir}/sars-cov-2_{ticket}_results.csv",
}

# Multiqc metrics to report to Vogue
MULTIQC_TO_VOGUE = {
    "multiqc_picard_AlignmentSummaryMetrics": {
        "format": "single",
        "step": "AlignmentSummaryMetrics",
        "input": "{key}.sorted.bam",
        "fields": {
            "TOTAL_READS": "TOTAL_READS",
            "PCT_PF_READS_ALIGNED": "PCT_PF_READS_ALIGNED",
            "MEAN_READ_LENGTH": "MEAN_READ_LENGTH",
        },
    },
    "multiqc_picard_insertSize": {
        "format": "single",
        "step": "CollectInsertSizeMetrics",
        "input": "{key}.mapped.primertrimmed.sorted.bam",
        "fields": {"MEDIAN_INSERT_SIZE": "MEDIAN_INSERT_SIZE"},
    },
    "multiqc_picard_wgsmetrics": {
        "format": "single",
        "step": "CollectWgsMetrics",
        "input": "{key}.mapped.primertrimmed.sorted.bam",
        "fields": {
            "MEDIAN_COVERAGE": "MEDIAN_COVERAGE",
            "PCT_10X": "PCT_10X",
            "PCT_30X": "PCT_30X",
            "PCT_50X": "PCT_50X",
            "PCT_100X": "PCT_100X",
        },
    },
    "multiqc_cutadapt": {
        "format": "paired",
        "step": "cutadapt",
        "input": "{key}_{direction}.fastq.gz",
        "fields": {"percent_trimmed": "percent_trimmed"},
    },
    "multiqc_general_stats": {
        "format": "paired",
        "step": "fastQC",
        "input": "{key}_{direction}.fastq.gz",
        "fields": {
            "FastQC_mqc-generalstats-fastqc-percent_duplicates": "percent_duplicates",
            "FastQC_mqc-generalstats-fastqc-percent_gc": "percent_gc",
            "FastQC_mqc-generalstats-fastqc-total_sequences": "raw_total_sequences",
        },
    },
}
