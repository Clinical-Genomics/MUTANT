#!/usr/bin/env bash

set -e

INSTANCE=$1
BRANCH=$2

OLDDIR=$(pwd)

if [ "$INSTANCE" != "production" ] && [ "$INSTANCE" != "stage" ]; then
    echo "Error: Please provide argument: production or stage"
    exit 1
fi

#INSTANCE='stage' or INSTANCE='production'
cd /home/proj/${INSTANCE}/mutant/MUTANT
git fetch
git checkout $BRANCH
git pull origin $BRANCH
INITIAL="$(echo $INSTANCE | head -c 1)"
source activate ${INITIAL^}_mutant
pip install -e .
source deactivate
cd mutant/externals/gms-artic
git fetch
git checkout master
git pull origin master
cd ${OLDDIR}

# Pull latest image
singularity pull --force /home/proj/${INSTANCE}/mutant/MUTANT/mutant/externals/gms-artic/artic-ncov2019-illumina.sif docker://clinicalgenomics/artic-ncov2019-illumina:latest
