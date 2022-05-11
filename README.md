# MUTANT
[Current Version](https://github.com/Clinical-Genomics/MUTANT/blob/main/mutant/__init__.py#L5)

Microbial Utility Toolbox And wrapper for data traNsmission and Transformation

## Installation

### Prereq
* Install Java 8 or above (`java --version`)
* Install [docker](https://www.docker.com/).
* Alternatively, install singularity: 
  * On linux: `conda install -c conda-forge singularity`
  * On Mac/Windows follow the installation instructions [here](https://sylabs.io/guides/3.7/admin-guide/installation.html#installation-on-windows-or-mac).
    * See information below on using vagrant for running singularity.

### Setup
* `git clone --recurse-submodules https://github.com/Clinical-Genomics/MUTANT.git`
* `cd MUTANT`
* `source setup.sh D_mutant`
* `source activate D_mutant`
* Pull the latest image from Dockerhub or build images as described below.

## Development
* Put your changes in a branch. Make a Pull Request of said branch.
* Run the `self-test`.
* Fill in everything you can in the PR template.
* Ping colleague and ask them to review your PR & suggest changes.
* Integrate said changes. If the changes are only aestetic, simply mention you're done without rerunning any tests.
* Merge your changes and version bump. 
* Deploy new version in production.


### Self-test
* `source activate D_mutant`
* `cd MUTANT`
* `mutant analyse sarscov2 tests/testdata/fasta_files --profiles local,docker --config_artic mutant/config/local/default_config.json --config_mutant mutant/config/local/mutant.json --config_case tests/testdata/MIC3109_artic.json` 
* Wait for pipeline completion (~3m). Check results in `./results/` 

Or install MUTANT under S_mutant and run `cg workflow mutant start frankhusky` with results in `/home/proj/stage/mutant/cases/frankhusky`

### Version bumping
MUTANTs versioning is bumped manually post PR merge, using [semver](https://semver.org/) standards on [this](https://github.com/Clinical-Genomics/MUTANT/blob/main/mutant/__init__.py#L3) variable.

## Deploying new version
Deploy new version in **production**:
* `bash mutant/standalone/deploy_hasta_update.sh production main`

## Docker containers
Containers are available on [dockerhub](https://hub.docker.com/r/clinicalgenomics/artic-ncov2019-illumina). 
The containers can be built locally (not possible on hasta) using the instructions below.

### Build docker containers
1. Start docker.
2. `cd mutant/externals/gms-artic/`
3. `docker build -f environments/illumina/Dockerfile -t artic-ncov2019-illumina:<version> .`

### Push docker containers to Dockerhub
The following script will push `artic-ncov2019-illumina:<version>` to Dockerhub.
* `docker login --username=<>`
* `bash mutant/standalone/push_docker_image.sh <version>`

### Pull containers from Dockerhub
* Automatically during analysis: Set container to `clinicalgenomics/artic-ncov2019-illumina`
* or pull singularity image manually: `singularity pull <MUTANTDIR>/mutant/externals/gms-artic/artic-ncov2019-illumina.sif docker://clinicalgenomics/artic-ncov2019-illumina`
* or pull docker image manually: `docker pull clinicalgenomics/artic-ncov2019-illumina`

## Singularity containers

### Build singularity containers
1. If needed, set up vagrant vm as below and perform `Setup`.
2. `mutant toolbox create-images`

#### Basics in using vagrant
* Set up the Vagrant virtual machine in an **empty directory**. This directory will be shared to the path `/vagrant` in the vm.
```
export VM=sylabs/singularity-ce-3.8-ubuntu-bionic64
vagrant init $VM
vagrant up
vagrant ssh
```
* Install conda.
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh -b && miniconda3/bin/conda init bash && source ~/.bashrc
```
**Basic vagrant commands**:
* Configure the vm using the `Vagrantfile`
* Exit the vm: `ctrl + d`
* Destroy the vm: `vagrant destroy`
* Manage your boxes: `vagrant box --help`
