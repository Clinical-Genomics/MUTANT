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
git submodule init
git submodule update
#git submodule foreach 'git fetch origin; git checkout $(git rev-parse --abbrev-ref HEAD); git reset --hard origin/$(git rev-parse --abbrev-ref HEAD); git submodule update --recursive; git clean -dfx'

# Pull latest image
singularity pull --force /home/proj/${INSTANCE}/mutant/MUTANT/mutant/externals/gms-artic/artic-ncov2019-illumina.sif docker://genomicmedicinesweden/gms-artic-illumina:0c599b9
singularity pull --force /home/proj/${INSTANCE}/mutant/MUTANT/mutant/externals/gms-artic/artic-ncov2019-pangolin.sif docker://genomicmedicinesweden/gms-artic-pangolin:latest
chmod u=rwx,g=rx,o=rx mutant/externals/gms-artic/artic-ncov2019-illumina.sif
chmod u=rwx,g=rx,o=rx mutant/externals/gms-artic/artic-ncov2019-pangolin.sif
