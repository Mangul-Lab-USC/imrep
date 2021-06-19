# ImReP

![](https://github.com/mandricigor/imrep/blob/master/figure_imrep.png)

## Quick Start 
ImReP is a method for rapid and accurate profiling of the adaptive immune repertoires from regular RNA-Seq data
 
Download ImReP using 
```
git clone https://github.com/mandricigor/imrep.git
```

Install ImReP from the base directory

```
cd imrep
./install.sh 
```
Run ImReP analysis by a single command for the BAM file with mapped and unmapped reads (preferred). BAM file needs to be indexed (.bai file). Forgot to save unmapped reads, we got you covered. Learn more [here](https://github.com/mandricigor/imrep/wiki/Forgot-to-save-unmapped-reads%3F)

```
python imrep.py --bam example/toyExample.bam example/toyExample.cdr3
```

Note : ImReP is written in Python2.7. So if your cluster default Python version is not 2.7, then you will have to first load it and then compile the above command.

Find ImReP analysis in _toyExample_ directory. Learn more [here](https://github.com/mandricigor/imrep/wiki/Quick-Start) 


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





