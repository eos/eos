**********************
Command-Line Interface
**********************

.. _cli:

.. note::

   The EOS command-line interface is completely optional and does not provide any functionality beyond the interactive Python interface.

Although using EOS within an interactive Jupyter notebook on your personal computer or laptop
is useful to prototype an analysis, this approach sometimes suffers from limited computing power.
To circumvent this problem, you can alternatively

  * use EOS in Jupyter interactively on a remote workstation computer via an SSH tunnel (see the `FAQ <../faq>`_);
  * use EOS on remote workstations or compute clusters via the command-line interface.

Working in the command-line interface requires that all analyses are defined within an analysis file.
For an example of an analysis file see the corresponding `section in the user guide <user-guide/analysis-file-format-example.html>`_.
For the description of the file format see the corresponding `section in the reference <reference/analysis-file-format.html>`_.
The following provide a documentation of the command-line interface, which is available through a single script called ``eos-analysis``.

.. argparse::
   :filename: ../src/scripts/eos-analysis
   :func: _parser
   :prog: eos-analysis
   :nodescription:
   :nodefault:

   The ``eos-analysis`` script provides several subcommands that

    * inspect the analysis file;
    * sample from a posterior density with a variety of methods, including Metroplis-Hastings, Population Monte Carlos, and nested sampling;
    * perform auxiliary tasks on intermediate results.

   The output of these commands are stored on disk as directories filled with YAML files
   (for descriptions and small numerical datasets) and Numpy datafiles (for samples).
   The datafiles can be access with the classes documented as part of the :obj:`eos.data` module.
