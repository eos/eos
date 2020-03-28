##########
Python API
##########

**************
Module ``eos``
**************

EOS provides its basic functionality via the main ``eos`` module.

.. autoclass:: eos.Analysis
   :members:

   .. _eos-Analysis-prior-descriptions:

   Each prior description is a dictionary with the following mandatory elements:

   * **type** (*str*) -- The type specification of the prior. Must be one of ``uniform``, ``flat``, or ``gaussian``.
   * **parameter** (*str*) -- The name of the parameter for which the prior shall apply.
   * **min** (*float*) -- The lower boundary of the prior's support.
   * **max** (*float*) -- The upper boundary of the prior's support.

   A ``uniform`` or ``flat`` prior does not require any further description. A ``gaussian`` prior requires in addition
   providing the following two elements:

   * **central** (*float*) -- The median value of the parameter.
   * **sigma** (*float*, or *list*, *tuple* of *float*) -- The width of the 68% probability interval. If a list or tuple
     of two numbers is provided, the prior will by a asymmetric but continuous. The two values are then taken to be the
     distance to the lower and upper end of the 68% probability interval.


.. autoclass:: eos.BestFitPoint
   :members:

.. autoclass:: eos.GoodnessOfFit
   :members:

.. autoclass:: eos.Kinematics
   :members:
   :inherited-members:

.. autoclass:: eos.LogLikelihood
   :members:

.. autoclass:: eos.LogPosterior
   :members:

.. autoclass:: eos.LogPrior
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
Module ``eos.data``
*******************

EOS provides access to save and load the various (intermediate) results of analyses via the ``eos.data`` module.

.. autoclass:: eos.data.MarkovChain
   :members:

.. autoclass:: eos.data.MixtureDensity
   :members:

.. autoclass:: eos.data.PMCSampler
   :members:

*******************
Module ``eos.plot``
*******************

EOS provides a plotting framework based on `Matplotlib <https://matplotlib.org/>`_.
Plots can readily be created from a Python script, from within a Jupyter notebook,
or in the command line using the ``eos-plot`` script.
For all of these cases a description of the plot is required in the format described below.
For the command-line script ``eos-plot``, the Python dictionary describing the plots must be provided as a YAML file.

.. note::

   Import ``eos.plot`` before you do something like ``import matplotlib.pyplot as plt``,
   because the ``eos.plot`` module sets its default plot style and a matplotlib backend.
   All options (except the backend) can be overwritten by updating ``matplotlib.rcParams[...]``;
   see also the ``matplotlib`` documentation.
   Note that the default settings use LaTeX to create labels and math expressions,
   so for this to work latex needs to be available on your system.

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

 * ``label`` (*str*, may contain LaTeX commands) -- The axis' label.
 * ``unit``  (*str*, may contain LaTeX commands) -- The axis' unit, which will be appended to the axis' label in square brackets.
 * ``range`` (*list* or *tuple* of two *float*) -- The tuple of [minimal, maximal] values, which will be displayed along the axis.

An example illustrating plot layouting follows:

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

 * ``type`` (*str*, mandatory) -- The type of the plot item, from one of the following recognized item types:

    - ``observable`` -- See `Plotting observables`_ for the description of this type of plot item.
    - ``constraints`` -- See `Plotting constraints`_. for the description of this type of plot item.

 * ``name`` (*str*, optional) -- The name of the plot item, for convenience when reporting warnings and errors.

All item types recognize the following optional keys:

 * ``alpha`` (*float*, between 0.0 and 1.0) -- The opacity of the plot item expressed as an alpha value. 0.0 means completely transparent,
   1.0 means completely opaque.

 * ``color`` (*str*, containing any valid Matplotlib color specification) -- The color of the plot item.
   Defaults to one of the colors in the Matplotlib default color cycler.
 * ``label`` (*str*, may contain LaTeX commands) -- The label that appears in the plot's legend for this plot item.

Plotting observables
--------------------

Contents items of type ``observable`` are used to display one of the built-in `observables <../observables.html>`_.
The following keys are mandatory:

 * ``observable`` (:class:`QualifiedName <eos.QualifiedName>`) -- The name of the observable that will be plotted.
   Must identify one of the observables known to EOS; see `the complete list of observables <../observables.html>`_.
 * ``range`` (*list* or *tuple* of two *float*) --The tuple of [minimal, maximal] values of the specified kinematic variable
   for which the observable will be evaluated.

Exactly one of the following keys is mandatory, to specify either a kinematic variable or a parameter to which the x coordinate
will be mapped:

 * ``variable`` (*str*) -- The name of the kinematic variable to which the x axis will be mapped.
 * ``kinematic`` (*str*) -- Alias for ``variable``.
 * ``parameter`` (*str*) -- The name of the parameter to which the x axis will be mapped;
   see `the complete list of parameters <../parameters.html>`_.

Example:

.. code-block::

   plot_args = {
       'plot': { ... },
       'contents': [
           {
               'label': r'$\ell=\mu$',
               'type': 'observable',
               'observable': 'B->Dlnu::dBR/dq2;l=mu',
               'variable': 'q2',
               'range': [0.02, 11.60],
           },
       ]
   }

Plotting constraints
--------------------

Contents items of type ``constraints`` are used to display one of the built-in `experimental or theoretical constraints <../constraints.html>`_.
The following keys are mandatory:

 * ``constraints`` (:class:`QualifiedName <eos.QualifiedName>` or iterable thereof) -- The name or the list of names of the constraints
   that will be plotted. Must identify at least one of the constraints known to EOS; see `the complete list of constraints <../constraints.html>`_.
 * ``variable`` (*str*) -- The name of the kinematic variable to which the x axis will be mapped.

When plotting multivariate constraints, the following key is also mandatory:

 * ``observable`` (:class:`QualifiedName <eos.QualifiedName>`) -- The name of the observable whose constraints will be plotted.
   Must identify one of the observables known to EOS; see `the complete list of observables <../observables.html>`_.
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
