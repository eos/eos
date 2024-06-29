=========================
Figure Description Format
=========================

EOS provides a descriptive framework to produce production-quality figures, using YAML inputs.
Figure descriptions can be stored as standalone YAML files or be included within the ``figures`` section of an `Analysis file <./analysis-file-format.html>`_.
Each figure description specifies the type of figure, the plot or plots included in the figure, and the content items displayed in any given plot.


Figures
-------

All figure descriptions accept two top-level YAML keys:

 - ``name`` (**mandatory only when used within an analysis file**) --- The name of the figure. Must by a valid file name.
 - ``type`` (**mandatory**) --- The type of the figure.


Single-Plot Figure
^^^^^^^^^^^^^^^^^^

A figure of type ``single`` contains a singe plot. The plot description is stored under the ``plot`` key.


Grid Figure
^^^^^^^^^^^

A figure of type ``grid`` contains a rectangular grid of individual plots.
The shape of the grid is stores as ``(rows, columns)`` tuple in the ``shape`` key.
The individual plots are described through the elements of the ``plots`` list,
which is a row-major flattening of the grid.


Corner Figure
^^^^^^^^^^^^^

A corner figure is specified using the type ``corner``.
It accepts a list of one or more elements as its ``contents`` key,
with each element describing a data file produced by the `common tasks framework <./python.html#common-tasks>`_.
