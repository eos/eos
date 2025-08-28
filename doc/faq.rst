##########################
Frequently Asked Questions
##########################

This is a list of Frequently Asked Questions about EOS. Feel free to suggest new entries!

************
How Do I ...
************

.. _faq-virtualenv:

... use EOS within a virtual Python environment?
================================================

We highly recommend to use EOS within a **virtual Python environment**. If you modify the EOS software, we recommend that you place your virtual environment folder within the EOS root directory.
To create the virtual environment, make sure you have the ``venv`` module installed, and then create a new virtual environment:

::

    cd eos
    python3 -m venv .venv


The new virtual environment will now be **activated**.
That means, the ``python`` (and ``pip``) command will now execute the local binary instead of a system-wide one.
At any later point, you can (re-)activate the virtual environment by running the following command:

::

    source .venv/bin/activate


While activated, you can **install** a binary distribution of EOS using ``pip install eoshep``.
If you plan to install EOS from source, as discussed in the :ref:`installation section<installation-from-source>`, you must set the prefix
to match the path to your virtual environment: ``$VIRTUAL_ENV``.

You most likely want to install **jupyter and notebook**.
To be able to access the virtual environment from the notebook interface, you also need to *install it* (details `here <https://ipython.readthedocs.io/en/stable/install/kernel_install.html#kernels-for-different-environments>`_).
While activated, execute the following:

::

    pip install jupyter notebook
    python -m ipykernel install --user --name venv-eos --display-name "python3 with EOS"


You can now **use jupyter** with ``python -m jupyter notebook`` and choose *python3 with EOS* as a kernel.

To **deactivate** the virtual environment, run

::

    deactivate


Of course, you can choose a **different directory** than ``eos/.venv`` to host your virtual environment. Please see
the ``venv`` `documentation <https://docs.python.org/3/tutorial/venv.html>`_ for more information.


... use EOS remotely?
=====================

Since EOS can be used from within a Jupyter notebook, you can use it on a remote Linux computer.
We assume that remote computer has the hostname ```REMOTE_HOST```. We further assume that you
can connect to that computer through a secure shell, via ```ssh REMOTE_HOST``` or some similar
command. Last but not least, we assume that you have installed EOS on the remote computer, e.g.
via the Python Package Index (PyPI): ```pip3 install eoshep```.
To use EOS on the remote computer through your local computer's browser, you need to slightly
change your connection details.
Use ```ssh -L 8888:localhost:8888 REMOTE_HOST``` to connect. Now you can run ```jupyter-notebook --no-browser```,
which will provide you with some output that looks as follows

::

    To access the notebook, open this file in a browser:
        ...
    Or copy and paste one of these URLs:
        http://localhost:8888/?token=...
     or http://127.0.0.1:8888/?token=...


Copy one of these links, making sure to include the lengthy token, and navigate to it in the
web browser on your local computer. You will be able to use EOS within an interactive Python environment.

You can make this more persistent by configuring your SSH connection to ```REMOTE_HOST```. Please see the documentation
of your local computer's secure shell program.

... use and develop EOS on Windows?
===================================

EOS can be used on Windows through the Windows Subsystem for Linux (WSL),
which is offered and developed by Microsoft themselves.
The installation of WSL is described in detail `here <https://docs.microsoft.com/en-us/windows/wsl/install>`__.
The installation provides you with a Linux container in a virtual machine.
From here on, you can follow the instructions for the installation method of your choosing described for EOS under Linux.
If you would like to make adaptions to the code for a project using EOS and/or contribute to its development,
we recommend using the Visual Studio Code editor as your development environment.
It supports connecting natively to the WSL container if you follow these `steps <https://code.visualstudio.com/docs/remote/wsl>`__.

.. _faq-check-installation:

... check the installation?
===========================

You can test with a single command that EOS is available as a python module:

.. code-block:: bash

    python3 -c "import eos; print(eos.__doc__)"

Please verify that the command prints a description of the EOS module:

.. code-block:: none

    EOS is a software package that addresses several use cases in the field of high-energy flavor physics (HEP). [...]

If you built and installed from the source code
and the import fails, i.e. you see ``ModuleNotFoundError: No module named 'eos'``,
the module was either not installed with ``make install`` or the ``PYTHONPATH`` is not configured correctly.

If you installed via ``pip``
and the output description is incorrect or ``None``,
you have another module installed that is also named ``eos``.
Mind that the package name for the ``pip3 install`` command is ``eoshep``.

.. _faq-tags:

... choose which ``tag`` to use for :math:`b \to s \ell \ell` observables
=========================================================================

The contributions from non-local operators in :math:`b \to s \ell \ell` observables
can be estimated with different theoretical approaches, usually not valid in the same energy range,
based on different input parameters and different theoretical assumptions.
Depending on the analysis under consideration, the user can choose between these approaches by specifying an option called ``tag``.
As different ``tags`` hide completely different assumptions, we have made this option mandatory.

We collect here the different ``tags`` used in EOS and the corresponding references.

.. list-table::
   :widths: 15 55 30
   :header-rows: 1

   * - ``tag``
     - Description
     - Reference

   * - ``Naive``
     - Absence of any non-local contributions.
     -

   * - ``BFS2004``
     - Calculation of the non-local contribution based on the framework of QCD factorization
       and a perturbative treatment of the intermediate quark states.
       This approach is expected to yield sensible results only for :math:`q^2 \ll 4 m_c^2`.
     - [BFS:2001A], [BFS:2004A]

   * - ``GP2004``
     - Calculation of the non-local contributions in a local OPE. This approach is expected
       to yield sensible results only for :math:`q^2`-integrated observables above the open-charm threshold.
     - [GP:2004A]

   * - ``GvDV2020``
     - Parametrization based on unitarity and analyticity valid for small dilepton masses below the :math:`\psi(2S)` state.
     - [GvDV:2020A], [GRvDV:2022A]

.. _faq-verbose:

... get more information on what is going on
============================================

EOS verbosity is handeled by the `logging module <https://docs.python.org/3/library/logging.html>`__.
For the CLI, verbosity is set using the EOS_VERBOSITY environment variable and ranges from 0 (show only errors) to 5 (show all messages).
Similarly, verbosity can be set for an individual command using the ``-v`` or ``--verbose`` argument.
The same can be achive in a python notebook with the following trick:

.. code-block:: python

    import sys
    import logging

    eos.logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler(stream=sys.stdout)
    eos.logger.addHandler(handler)

Above, ``logging.DEBUG`` can be replaced by any logging level, which are listed `here <https://docs.python.org/3/library/logging.html#levels>`__.
