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
To use EOS on the remote computer through your local computers browser, you need to slightly
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
of your local computer's secure shell programm.

... use EOS on Windows?
=======================

Presently we do not support using EOS on Windows. We recommend that you install EOS on a
remote Linux machine (e.g., one of your department's workstations) and connect remotely
as discussed :ref:`above <... use EOS remotely?>`.
