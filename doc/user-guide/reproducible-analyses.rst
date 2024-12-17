*****************************************
Setting up a fully reproducible analysis
*****************************************

EOS now has the capability to run a fully reproducible analysis, using `snakemake <https://snakemake.readthedocs.io/en/stable>`_.
This means that the analysis file describes fully all the steps needed, including the relative dependencies, and can be easily reproduced by anyone with that analysis file and the same EOS version.

Below is a Snakefile that can be used with minimal modifications - the only change you need to make is to set the base directory for all the EOS output, via the ``basedir`` variable.

.. literalinclude:: ../../examples/Snakefile
   :language: default

If you have snakemake installed, you can run this Snakefile by running the ``snakemake`` command in the same directory as the Snakefile and your analysis file, which should be named ``analysis.yaml``.
See the `snakemake documentation <https://snakemake.readthedocs.io/en/stable>`_ for more information on how to use snakemake, in particular other features such as parallelising locally or batch running on a cluster.
