"""Parse SARS-CoV-2 results"""


import json
from mutant.constants.artic import MULTIQC_TO_VOGUE


def get_multiqc_template_data(record, field, step, direction="") -> dict:
    """Get default values for a multiqc metric"""

    info = MULTIQC_TO_VOGUE[step]
    sample_template = {
        "header": "~",
        "id": record["CG_ID_sample"],
        "input": info["input"].format(
            key="_".join(
                [
                    record["region_code"],
                    record["lab_code"],
                    record["Customer_ID_sample"],
                ]
            ),
            direction=direction,
        ),
        "name": info["fields"][field],
        "step": info["step"],
        "value": None,
    }
    return sample_template


def get_multiqc_sample_keys(data, sample, step, info) -> dict:
    """Get sample keys from multiqc.json data for an analysis step"""

    # Get sample keys
    sample_keys = {}
    # Get single key for sample
    if info["format"] == "single":
        for key in data["report_saved_raw_data"][step].keys():
            keycont = key.split("_")
            if sample in keycont:
                sample_keys["0"] = key
    # Get paired keys for sample
    elif info["format"] == "paired":
        for key in data["report_saved_raw_data"][step].keys():
            keycont = key.split("_")
            # Get keys containing sample name and "_1" or "_2" suffix
            if (sample in keycont) and (keycont[-1] in ["1", "2"]):
                suffix = keycont[-1]
                sample_keys[suffix] = key
    return sample_keys


def get_multiqc_metric_value(data, sample_keys, step, field, direction="0") -> str:
    """Get multiqc json data for metric value"""

    return data["report_saved_raw_data"][step][sample_keys[direction]][field]


def get_vogue_multiqc_data(multiqc, caseinfo) -> list:
    """Parse multiqc.json data for vogue delivery"""

    case_data = list()
    with open(multiqc) as f:
        data = json.load(f)
        for record in caseinfo:
            # sample_data = list()
            for step, info in MULTIQC_TO_VOGUE.items():
                # Get sample key naming in multiqc file
                if record["sequencing_qc_pass"]:
                    sample_keys = get_multiqc_sample_keys(
                        data, record["Customer_ID_sample"], step, info
                    )
                # Parse selected metrics in multiqc json
                for field, field_name in info["fields"].items():
                    # Add default values
                    if info["format"] == "single":
                        # Get template data
                        field_data = get_multiqc_template_data(record, field, step)
                        # Update value if analyzed
                        if record["sequencing_qc_pass"] and len(sample_keys) > 0:
                            field_data["value"] = get_multiqc_metric_value(
                                data, sample_keys, step, field
                            )
                        case_data.append(field_data)
                    else:
                        for direction in [1, 2]:
                            # Get template data
                            field_data = get_multiqc_template_data(
                                record, field, step, str(direction)
                            )
                            # Update value if analyzed
                            if record["sequencing_qc_pass"] and len(sample_keys) > 0:
                                field_data["value"] = get_multiqc_metric_value(
                                    data, sample_keys, step, field, str(direction)
                                )
                            case_data.append(field_data)
    return case_data
