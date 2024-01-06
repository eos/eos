*****************************************
Common Tasks in Phenomenological Analyses
*****************************************

A number of tasks are common to all phenomenological analyses.
These include

 * finding the best-fit point of a posterior distribution;
 * sampling from a posterior distribution;
 * and computing posterior predictive distributions.

EOS strives to provide a versatile interface to carry out these tasks.
All tasks are designed to be used in conjunction with other tasks.
To transfer information from one task to another, tasks store their (intermediate) results in various directories.
The root directory for this storage should be provided to the invocation of each task;
otherwise the tasks default to using the current working directory.

To ensure that a task can be used both from within a Jupyter notebook and from the command line,
every task is implemented as a function in the ``eos.tasks`` module.
We refer to the :ref:`EOS Python API documentation <api_eos_tasks>` for the full list of tasks and their arguments.
The command-line interface to the tasks is provided by the :command:`eos-analysis` command,
which is documented in the :ref:`Command-Line Interface <cli>` section.
