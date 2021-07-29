##############
Advanced Usage
##############

******************************
The EOS Command-Line Interface
******************************

Although using EOS within an interactive Jupyter notebook on your personal computer or laptop
is useful to prototype an analysis, this approach sometimes suffers from limited computing power.
To circumvent this problem, you can alternatively

  * use EOS in Jupyter interactively on a remote workstation computer via an SSH tunnel (see the `FAQ <faq>`_);
  * use EOS on remote workstations or compute clusters via the command-line interface.

In the following we document the command-line interface and the file format used in conjunction with it.

.. note::

   The EOS command-line interface is completely optional and does not provide any means beyond the
   interactive Python interface.

The Analysis Description Format
===============================

EOS uses a YAML file to describe the individual steps of one or more statistical analyses.
At the top level, the format includes the following YAML keys:

 - ``priors`` (**mandatory**) --- The list of priors within the analysis.
 - ``likelihoods`` (**mandatory**) --- The list of likelihoods within the analysis..
 - ``posteriors`` (**mandatory**) --- The list of posteriors within the analysis.
 - ``predictions`` (**optional**) --- The list of theory predictions within the analysis.

Describing Priors
~~~~~~~~~~~~~~~~~

The ``priors`` key contains a list of *named* priors. Each prior has two mandatory keys:

  - ``name`` (**mandatory**) --- The unique name of this prior.
  - ``parameters`` (**mandatory**) --- The ordered list of parameters described by this prior.

The description of each individual parameter follows the prior description used in the
:class:`Analysis <eos.Analysis>` constructor.


Describing Likelihoods
~~~~~~~~~~~~~~~~~~~~~~

The ``likelihoods`` key contains a list of *named* likelihoods. Each likelihood has two mandatory keys:

  - ``name`` (**mandatory**) --- The unique name of this likelihood.
  - ``constraints`` (**mandatory**) --- The ordered list of EOS constraint names that comprise this likelihood.

Describing Posteriors
~~~~~~~~~~~~~~~~~~~~~

The ``posteriors`` key contains a list of *named* posteriors. Each posterior can contain several keys:

  - ``name`` (**mandatory**) --- The unique name of this posterior.
  - ``global_options`` (**optional**) --- A key/value map providing global options, i.e., options that apply to all observables used by this posterior.
  - ``prior`` (**mandatory**) --- The ordered list of named priors that are used as part of this posterior.
  - ``likelihood`` (**optional**) --- The ordered list of named likelihoods that are used as part of this posterior.
  - ``fixed_parameter`` (**optional**) --- A key/value map providing values for parameters that deviate from the default values.

Example
~~~~~~~

.. toggle-header::
   :header: Example `examples/cli/btopilnu.analysis <https://github.com/eos/eos/tree/master/examples/cli/btopilnu.analysis>`_

   .. literalinclude:: ../examples/cli/btopilnu.analysis
      :language: YAML

|

The Command-Line Interface
==========================

.. argparse::
   :filename: ../src/scripts/eos-analysis
   :func: _parser
   :prog: eos-analysis
   :nodescription:
   :nodefault:

   The ``eos-analysis`` script provides several subcommands that

    * inspect the analysis file;
    * sample from a posterior density with Monte Carlo methods;
    * perform auxiliary tasks on intermediate results.

   The output of these commands are stored on disk as directories filled with YAML files
   (for descriptions and small numerical datasets) and Numpy datafiles (for samples).
   The datafiles can be access with the classes documented as part of the :obj:`eos.data` module.
