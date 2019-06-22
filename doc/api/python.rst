##########
Python API
##########

**************
Module ``eos``
**************

.. autoclass:: eos.Analysis
   :members:

.. autoclass:: eos.BestFitPoint
   :members:

.. autoclass:: eos.GoodnessOfFit
   :members:

.. autoclass:: eos.Kinematics
   :members:

.. autoclass:: eos.Observable
   :members:

.. autoclass:: eos.Observables
   :members:

.. autoclass:: eos.Options
   :members:

.. autoclass:: eos.Parameters
   :members:

.. autoclass:: eos.QualifiedName
   :members:

*******************
Module ``eos.plot``
*******************

EOS provides a plotting framework based on `Matplotlib <https://matplotlib.org/>`_.
Plots can readily be created from a Python script, from within a Jupyter notebook,
or in the command line using the ``eos-plot`` script.
For all of these cases a description of the plot is required in the format described below.
For the command-line script ``eos-plot``, the Python dictionary describing the plots must be provided as a YAML file.

.. autoclass:: eos.plot.Plotter
   :members:

Plot description format
=======================

The input must be formatted as a dictionary containing the keys ``plot`` and ``contents``.
The ``plot`` key must be mapped to a dictionary; it describes the layout of the plot,
including axis labels, positioning of the key, and similar settings.
The ``contents`` key must be mapped to a list; it describes the contents of the plot,
expressed in terms of independent plot items. 

.. code-block::

   plot_desc = {
       'plot': {
           'x': { ... }, # description of the x axis
           'y': { ... }, # description of the y axis
           'legend': { ... } # description of the legend
       },
       'contents': [
           { ... }, # first plot item
           { ... }, # second plot item
       ]
   }
   eos.plot.Plotter(plot_desc, FILENAME).plot()

In the above, ``FILENAME`` is an optional argument naming the file into which the plot shall be placed.
The format is automatically determined based on the file name extension.

Plot layouting
--------------

By default plots lack any axis labels and units, and legend.

Axis descriptions can be provided through the following key/value pairs, which apply equally to the x and y axis:

``label`` : string : (may contain LaTeX commands)
    The axis' label.
``unit``  : string : (may contain LaTeX commands)
    The axis' unit, which will be appended to the axis' label in square brackets.
``range`` : list or tuple of two floating point numbers
    The tuple of (minimal, maximal) values, which will be displayed along the axis.

An exampe illustrating plot layouting follows:

.. code-block::

   plot_args = {
       'plot': {
           'x': { 'label': r'$q^2$', 'unit': r'$\textnormal{GeV}^2$', 'range': [0.0, 11.60] },
           'y': { 'label': r'$d\mathcal{B}/dq^2$',                    'range': [0.0,  5e-3] },
           'legend': { 'location': 'upper center' }
       },
       'contents': [
           ...
       ]
   }

Plot contents
-------------

Each item in a plot's contents is represented by a dictionary. The only mandatory key is the item's ``type``.
Every item can feature an optional ``name``, which will be used when notifying the user about pertinent information,
warnings or errors.

``type`` : string : (mandatory)
    The type of the item, from one of the following recognized item types:

    ``observable``
        See `Plotting observables`_.
    ``constraints``
        See `Plotting constraints`_.

``name`` : string : (optional)
    The name of the item, for convenience when reporting warnings and errors.


All item types recognize the following optional keys:

``alpha`` : float between 0.0 and 1.0
    The opacity of the plot item expressed as an alpha value. 0.0 means completely transparent, 1.0 means completely
    opaque.

``color`` : string containing any valid Matplotlib color specification
    The color of the plot item. Defaults to one of the colors in the Matplotlib default color cycler.

``label`` : string : (may contain LaTeX commands)
    The label that appears in the plot's legend for this content item.

Plotting observables
--------------------

Contents items of type ``observable`` are used to display one of the built-in :class:`observables <eos.Observable>`.
The following keys are mandatory:

``observable`` : :class:`qualified name <eos.QualifiedName>`
    The name of the observable that will be plotted. Must be identify one of the observables known to EOS.

``range`` : list or tuple of two floating point numbers
    The tuple of (minimal, maximal) values of the specified kinematic variable for which the observable will be evaluated.

Exactly one of the follow keys is mandatory, to specify either a kinematic variable or a parameter to which the x coordinate will be mapped:

``variable`` : string
    The name of the kinematic variable to which the x axis will be mapped.

``parameter`` : string representation of a :class:`qualified name <QualifiedName>`
    The name of the :class:`parameter <eos.Parameter>` to which the x axis will be mapped.

Example:

.. code-block::

   plot_args = {
       'plot': { ... },
       'contents': [
           {
               'label': r'$\ell=\mu$',
               'type': 'observable',
               'observable': 'B->Dlnu::dBR/dq2;l=mu',
               'kinematic': 'q2',
               'range': [0.02, 11.60],
           },
       ]
   }

Plotting constraints
--------------------

Contents items of type ``constraints`` are used to display one of the built-in :class:`experimental or theoretical constraints <eos.Constraint>`.
The following keys are mandatory:

``constraints`` : :class:`qualified name <eos.QualifiedName>` or list thereof
    The name or the list of names of the constraints that will be plotted. Must identify at least one of the constraint known to EOS.

``variable`` : string
    The name of the kinematic variable to which the x axis will be mapped.

When plotting multivariate constraints, the following key is also mandatory:

``observable``  :class:`qualified name <eos.QualifiedName>`
    The name of the observable whose constraints will be plotted. Must be identify one of the observables known to EOS.
    This is only mandatory in multivariate constraints, since these can constrain more than one observable simultaneously.

Example:

.. code-block::

   plot_args = {
       'plot': { ... },
       'contents': [
           {
               'label': r'Belle 2015 $\ell=e,\, q=d$',
               'type': 'constraint',
               'color': 'C0',
               'constraints': 'B^0->D^+e^-nu::BRs@Belle-2015A',
               'observable': 'B->Dlnu::BR',
               'variable': 'q2',
               'rescale-by-width': False
           }
       ]
   }
