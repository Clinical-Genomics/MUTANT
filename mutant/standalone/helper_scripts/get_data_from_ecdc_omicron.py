"""Parse SARS-CoV-2 classifications from ECDC
Adaptation of get_data_from_ecdc.py for updated table structure at ECDC with manual handling of omicron
"""
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
classtype = ["VOC", None, "VOI", "VUM", "De-escalated"]
all_mutations = list()

with open(os.path.join(outdir, "classifications.csv"), "w") as class_outfile:
    # Set up output
    classifications = csv.writer(class_outfile, lineterminator="\n")
    classifications.writerow(["lineage", "spike", "class"])

    # Scrape variants from each of the 4 tables of different classes
    for i in [0, 2, 3, 4]:
        class_lineage = list()

        # Parse lineage data
        for row in soup.findAll("table")[i].tbody.findAll("tr"):
            # Skip malformated Omicron row
            if row.findAll("td")[0].getText().strip() == "Omicron":
                continue
            # Skip malformated rows
            if len(row.findAll("td")) != 9:
                continue
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

        # Add omicron data
        if i==0:
            class_lineage.append("B.1.1.529")

        # Output formatted lineages
        for lineage in class_lineage:
            spike = ""
            if "+" in lineage:
                spike = lineage.split("+")[1]
                lineage = lineage.split("+")[0]
            classifications.writerow([lineage, spike, classtype[i]])

# Add omicron mutations
all_omicron_mutations = "A67V, Δ69-70, T95I, G142D, Δ143-145, Δ211-212, ins214EPE, G339D, S371L, S373P, S375F, K417N," \
                       "N440K, G446S, S477N, T478K, E484A, Q493R, G496S, Q498R, N501Y, Y505H, T547K, D614G, H655Y," \
                       "N679K, P681H, N764K, D796Y, N856K, Q954H, N969K, L981F".replace(" ", "").split(",")
for mut in all_omicron_mutations:
    if mut not in all_mutations:
        all_mutations.append(mut)

with open(os.path.join(outdir, "spike_mutations.csv"), "w") as mutations_outfile:
    mutations = csv.writer(mutations_outfile, lineterminator="\n")
    mutations.writerow(["Spike_mutations_of_interest"])
    for mut in all_mutations:
        mutations.writerow([mut])
