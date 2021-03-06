Installation
============
In order to use adVNTR, it is recommended to (1) install adVNTR using conda packaging manager and (2) download the
predefined models for human genome from :ref:`data-requirements` section. However, you can install it from the source
and/or use custom models.

.. _install-with-conda:

Install adVNTR with conda
-------------------------
If you are using the conda packaging manager (e.g. `miniconda <https://conda.io/miniconda.html>`_ or anaconda),
you can install adVNTR from the `bioconda  channel <https://bioconda.github.io/>`_:

::

    conda config --add channels bioconda
    conda install advntr

adVNTR could be invoked from command line with ``advntr``


.. _data-requirements:

Data Requirements
-----------------
* To run adVNTR on trained VNTR models:
    - Download `vntr_data.zip <https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data.zip>`_ and extract it inside the project directory.
Alternatively, you can add model for custom VNTR. See :ref:`add-custom-vntr-label` for more information.

Execution:
----------
Use following command to see the help for running the tool.

::

    advntr --help

The program outputs the RU count genotypes of VNTRs. To specify a single VNTR by its ID use ``--vntr_id <id>`` option.
The list of some known VNTRs and their ID is available at `Disease-linked-VNTRs page <https://github.com/mehrdadbakhtiari/adVNTR/wiki/Disease-linked-VNTRs>`_ in wiki.

See the demo execution here or :ref:`quickstart` page to see an example data set with step-by-step genotyping commands.

Demo: input in BAM format
-------------------------
* ``--alignment_file`` specifies the alignment file containing mapped and unmapped reads:

::

    advntr genotype --alignment_file aligned_illumina_reads.bam --working_directory ./log_dir/

* With ``--pacbio``, adVNTR assumes the alignment file contains PacBio sequencing data:

::

    advntr genotype --alignment_file aligned_pacbio_reads.bam --working_directory ./log_dir/ --pacbio

* Use ``--frameshift`` to find the possible frameshifts in VNTR:

::

    advntr genotype --alignment_file aligned_illumina_reads.bam --working_directory ./log_dir/ --frameshift


Install from source (Not recommended)
-------------------------------------
You can also install adVNTR from the source instead of using conda.

Dependencies
^^^^^^^^^^^^
1. Following libraries are required
    -   ``python2.7``
    -   ``python-pip``
    -   ``python-tk``
    -   ``libz-dev``
    -   ``samtools``
    -   ``muscle``

You can install these requirement in Ubuntu Linux by running ``sudo apt-get install python2.7 python-pip python-tk libz-dev samtools muscle``

2. Following python2.7 packages are required:
    -   ``biopython``
    -   ``pysam`` version 0.9.1.4 or above
    -   ``cython``
    -   ``networkx`` version 1.11
    -   ``scipy``
    -   ``joblib``
    -   ``scikit-learn``

You can install required python libraries by running ``pip install -r requirements.txt``

3. In addition, ``ncbi-blast`` version 2.2.29 or above is required:
    - The easiest way to install ``ncbi-blast`` is to download `BLAST+ executables <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`_ and add them in the path.

To Install
^^^^^^^^^^
To get the latest version and install it locally, run:

::

    git clone https://github.com/mehrdadbakhtiari/adVNTR
    cd adVNTR
    python setup.py install

adVNTR could be invoked from command line with ``advntr``
