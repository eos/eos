##########################
Frequently Asked Questions
##########################

This is a list of Frequently Asked Questions about EOS. Feel free to suggest new entries!

************
How Do I ...
************

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

If you installed via `pip`
and the output description is incorrect or ``None``,
you have another module installed that is also named `eos`.
Mind that the package name for the ``pip3 install`` command is ``eoshep``.
