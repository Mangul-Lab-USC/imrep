# ImReP

![](https://github.com/mandricigor/imrep/blob/master/figure_imrep.png)

## Quick Start 
ImReP is a method for rapid and accurate profiling of the adaptive immune repertoires from regular RNA-Seq data
 
Download ImReP using 
```bash
git clone https://github.com/mandricigor/imrep.git
```

Environment Setup
Create the imrep Conda Environment:
```bash
conda create -n imrep python=2.7
```
This command creates a new conda environment named imrep with Python version 2.7.

Activate the imrep Environment:
```bash
conda activate imrep
```
Activate the newly created environment to use Python 2.7 and install the necessary packages.

Installation and Setup
Build and Install the Python Package: Navigate to the suffix_tree-2.1 directory and run the following commands:
```bash
cd suffix_tree-2.1
python2 setup.py build
python2 setup.py install
```
These commands build and install the Python package in the current directory.

Run the Installation Script: Execute the install.sh script in the parent directory:
```bash
cd ..
sh install.sh
```
This script performs additional setup tasks required for the environment.

Install Additional Packages: Install the pysam and biopython packages in the imrep2 environment:
```bash
conda install pysam
conda install biopython
```
These packages are necessary for the imrep.py script to function correctly.

Running the imrep.py Script
Display Help Message: To understand the usage of the imrep.py script, run:
```bash
python2 imrep.py -h
```
This command displays the help message for the script.

Process a FASTQ File: Navigate to the example directory and run the imrep.py script with a FASTQ file named toyExample.fastq, outputting to output1:
```bash
cd example
python2 ../imrep.py --fastq toyExample.fastq output1
```
This command processes the FASTQ file and generates output in the specified directory.

Additional Information
Python Version: This environment is specifically configured for Python 2.7, as required by the imrep.py script.
Conda Environment: The imrep environment is isolated from other conda environments, ensuring that the dependencies of this project do not interfere with other projects.
Dependencies: The environment includes Python 2.7, pysam, and biopython, along with any other dependencies required by the imrep.py script and related tasks.


Find ImReP analysis in _toyExample_ directory. Learn more [here](https://github.com/mandricigor/imrep/wiki/Quick-Start) 


## Docker Usage
Docker enables you to run ImRep without needing to explicitly install and manage its dependencies. You can learn more about Docker [here](https://docs.docker.com/get-started/overview/).

First, build the docker container from the current working directory:
```bash
git clone https://github.com/Mangul-Lab-USC/imrep.git
# Make sure you are in the correct directory...
cd imrep
# build the docker container and give it a <tag> `-t`
docker build -t imrep .
```
Run the docker container like so, and pass in your intended ImRep `<options>`.
```bash
docker run --rm -v ${PWD}:/imrep -ti imrep <options>
```
In this case, we mount a Docker [volume](https://docs.docker.com/storage/volumes/) `-v` to recieve the output generated from within the container locally.

Here is a full example with Docker:
```bash
docker build -t imrep .
docker run --rm -v ${PWD}:/imrep -ti imrep --bam example/toyExample.bam example/toyExample.cdr3
```
# ImReP Tutorial

Use the sidebar on the right to navigate ImReP tutorial. Get started with a toy example of 200 RNA-Seq reads (_example/toyExample.bam_) distributed with ImRep package


# What is ImReP

ImReP is a method to quantify individual immune response based on a recombination landscape of genes encoding B and T cell receptors (BCR and TCR).  ImReP is able to efficiently extract TCR and BCR reads from the RNA-Seq data and assemble clonotypes (defined as clones with identical CDR3 amino acid sequences)  and detect corresponding V(D)J recombinations. Using CAST clustering technique, ImReP is able to correct assembled clonotypes for PCR and sequencing errors.

# Publication

Mangul, Serghei, et al. ["Profiling adaptive immune repertoires across multiple human tissues by RNA Sequencing"](http://www.biorxiv.org/content/early/2017/03/25/089235) bioRxiv (2016): 089235. 

# Presentations

Slides from my talk at ISMB 2017 (HiTSeq 2017) are available [here](https://sergheimangul.files.wordpress.com/2017/07/hitseq_2017_public.pdf)


## The Atlas of Immune Repertoires (TheAIR)

Using ImReP we have created the Atlas of Immunoglobilin Repertoires (TheAIR). The AIR is a collection of CDR3 from Immunoglobulin (Ig) receptor repetoires. CDR3s are assembled from 8,555 human RNA-seq samples across 544 individuals from 53 tissues from the Genotype-Tissue Expression (GTEx v6) project using ImReP.The AIR has one of the largest collection of CDR3 sequences (n=3.6 million) and tissue types (n=53). TheAIR is freely available at : https://smangul1.github.io/TheAIR/


# Contact 

This software was developed by [Igor Mandric](https://github.com/mandricigor) and [Serghei Mangul](https://github.com/Mangul-Lab-USC/). Please do not hesitate to contact us (mandric.igor@gmail.com, mangul@usc.edu) if you have any comments, suggestions, or clarification requests regarding the tutorial or if you would like to contribute to this resource.





