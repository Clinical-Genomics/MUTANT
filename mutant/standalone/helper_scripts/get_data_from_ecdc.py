"""Parse SARS-CoV-2 classifications from ECDC"""
import os
import sys
import csv
import requests
from bs4 import BeautifulSoup

# For future implementations, get data from https://github.com/erikalmecdc/ecdc_virology

# Output directory
outdir = sys.argv[1]

# Request from URL
ecdc_url = "https://www.ecdc.europa.eu/en/covid-19/variants-concern/"
html_text = requests.get(ecdc_url).text
soup = BeautifulSoup(html_text, features="lxml")

# Initialize
classtype = ["VOC", "VOI", "VUM", "De-escalated"]
all_mutations = list()

with open(os.path.join(outdir, "classifications.csv"), "w") as class_outfile:
    # Set up output
    classifications = csv.writer(class_outfile, lineterminator="\n")
    classifications.writerow(["lineage", "spike", "class"])

    # Scrape variants from each of the 4 tables of different classes
    for i in range(4):
        class_lineage = list()

        # Parse lineage data
        for row in soup.findAll("table")[i].tbody.findAll("tr"):
            # Parse mutations of interest
            spike_muts = (
                row.findAll("td")[3].getText().strip().replace(" ", "").split(",")
            )
            for mut in spike_muts:
                # Remove comments in parenthesis
                if "(" in mut:
                    mut = mut.split("(")[0]
                if mut not in all_mutations:
                  all_mutations.append(mut)

            # Get data from variant column
            lineage = row.findAll("td")[1].getText().strip().replace(" ", "")
            # Remove comments in parenthesis
            if "(" in lineage:
                lineage = lineage.split("(")[0]
            # Split multiple synonymous variants
            if "/" in lineage:
                class_lineage.extend(lineage.split("/"))
            # Add single variant
            else:
                class_lineage.append(lineage)

        # Output formatted lineages
        for lineage in class_lineage:
            spike = ""
            if "+" in lineage:
                spike = lineage.split("+")[1]
                lineage = lineage.split("+")[0]
            classifications.writerow([lineage, spike, classtype[i]])

with open(os.path.join(outdir, "spike_mutations.csv"), "w") as mutations_outfile:
    mutations = csv.writer(mutations_outfile, lineterminator="\n")
    mutations.writerow(["Spike_mutations_of_interest"])
    for mut in all_mutations:
        mutations.writerow([mut])
