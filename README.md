[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/advntr/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/advntr/badges/downloads.svg)](https://anaconda.org/bioconda/advntr)
[![Build Status](https://travis-ci.org/mehrdadbakhtiari/adVNTR.svg?branch=master)](https://travis-ci.org/mehrdadbakhtiari/adVNTR)

adVNTR - A tool for genotyping VNTRs
------------------------------------
[adVNTR](https://github.com/mehrdadbakhtiari/adVNTR/) is a tool for genotyping Variable Number Tandem Repeats (VNTR)
from sequence data. It works with both NGS short reads (Illumina HiSeq) and SMRT reads (PacBio) and finds
diploid repeating counts for VNTRs and identifies possible mutations in the VNTR sequences.

Installation
------------
If you are using the conda packaging manager (e.g. [miniconda](https://conda.io/miniconda.html) or anaconda),
you can install adVNTR from the [bioconda  channel](https://bioconda.github.io/):

    conda config --add channels bioconda
    conda install advntr

adVNTR could be invoked from command line with ``advntr``

Alternatively, you can install dependencies and [install the adVNTR from source](http://advntr.readthedocs.io/en/latest/installation.html#install-from-source-not-recommended).


Data Requirements
-----------------
* To run adVNTR on trained VNTR models:
    - Download [vntr_data.zip](https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data.zip) and extract it inside the project directory.

Alternatively, you can add model for custom VNTR. See [Add Custom VNTR](http://advntr.readthedocs.io/en/latest/tutorial.html#add-custom-vntr-label) for more information.

Execution:
----------
Use following command to see the help for running the tool.

    advntr --help

The program outputs the RU count genotypes of VNTRs. To specify a single VNTR by its ID use ``--vntr_id <id>`` option.
The list of some known VNTRs and their ID is available at [Disease-linked-VNTRs page](https://github.com/mehrdadbakhtiari/adVNTR/wiki/Disease-linked-VNTRs) in wiki.

Demo 1: input in BAM format
---------------------------
* ``--alignment_file`` specifies the alignment file containing mapped and unmapped reads:

```sh
    advntr genotype --alignment_file aligned_illumina_reads.bam --working_directory ./log_dir/
```

* With ``--pacbio``, adVNTR assumes the alignment file contains PacBio sequencing data:

```sh
    advntr genotype --alignment_file aligned_pacbio_reads.bam --working_directory ./log_dir/ --pacbio
```

* Use ``--frameshift`` to find the possible frameshifts in VNTR:

```sh
    advntr genotype --alignment_file aligned_illumina_reads.bam --working_directory ./log_dir/ --frameshift
```

Citation:
---------
Bakhtiari, M., Shleizer-Burko, S., Gymrek, M., Bansal, V. and Bafna, V., 2017. [Targeted Genotyping of Variable Number Tandem Repeats with adVNTR](https://doi.org/10.1101/221754/). bioRxiv, p.221754.
