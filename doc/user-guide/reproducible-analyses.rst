*****************************************
Setting up a fully Reproducible Analysis
*****************************************

EOS now has the capability to run a fully reproducible analysis, using `snakemake <https://snakemake.readthedocs.io/en/stable>`_.
This means that the analysis file describes fully all the steps needed, including the relative dependencies, and can be reproduced by anyone with that analysis file and the same EOS version.
