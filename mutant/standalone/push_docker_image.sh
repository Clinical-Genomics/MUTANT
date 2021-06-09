#!/usr/bin/env bash

set -e

VERSION=$1

echo "Creating tags for artic-ncov2019-illumina:$VERSION."
docker tag artic-ncov2019-illumina:$VERSION clinicalgenomics/artic-ncov2019-illumina:$VERSION
docker tag artic-ncov2019-illumina:$VERSION clinicalgenomics/artic-ncov2019-illumina:latest
echo "Pushing artic-ncov2019-illumina:$VERSION to clinicalgenomics Dockerhub."
docker push clinicalgenomics/artic-ncov2019-illumina:$VERSION
docker push clinicalgenomics/artic-ncov2019-illumina:latest

echo "Done!"
