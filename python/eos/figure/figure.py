# Copyright (c) 2023 Danny van Dyk
# Copyright (c) 2023 Philip Lueghausen
#
# This file is part of the EOS project. EOS is free software;
# you can redistribute it and/or modify it under the terms of the GNU General
# Public License version 2, as published by the Free Software Foundation.
#
# EOS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 59 Temple
# Place, Suite 330, Boston, MA  02111-1307  USA

from .plot import PlotFactory

import matplotlib.pyplot as plt
import numpy as np


class Figure:
    def __init__(self, figure, axes, descriptions:list):
        self.figure = figure
        self.axes = axes
        self.plots = []

        if not descriptions:
            raise ValueError(f"Argument 'descriptions' is empty list")

        for ax, plot_description in zip(axes, descriptions):
            if plot_description is None:
                self.plots.append(PlotFactory.make(ax, [], plot_type='empty'))
                continue

            if 'description' not in plot_description:
                raise ValueError(f"Plot does not contain element 'description'")
            if 'items' not in plot_description:
                raise ValueError(f"Plot does not contain element 'items'")

            contents = plot_description['items']
            plot_description = PlotFactory.sanitize(plot_description['description'])

            self.plots.append(PlotFactory.make(ax, contents, **plot_description))


    def draw(self):
        for plot in self.plots:
            plot.draw()


class SingleFigure(Figure):
    """Represent a figure with a single set of axes.

    Use the value of 'plot' as a set of arguments to construct a Figure.

    :param plot: The plot's data, as a dictionary containing the 'description' and the 'contents'.
    :type plot: dict
    """

    def __init__(self, description):
        figure, axes = plt.subplots()
        super().__init__(figure, [axes,], [description,])


class GridFigure(Figure):
    """Represent a figure with several axes arranged in a grid.

    :param plots: The list of plots' data, as a list of dictionaries each containing the 'description' and the 'contents' for a single plot.
    :type plots: list[dict] or iterable[dict]
    :param nrows: The number of rows in the grid.
    :type nrows: int
    :param ncols: The number of columns in the grid.
    :type ncols: int
    """

    def __init__(self, descriptions, nrows, ncols, **kwargs):
        figure, axes = plt.subplots(nrows, ncols, **kwargs)
        axes = axes.flatten('C') # row-major style
        super().__init__(figure, np.array(axes), np.array(descriptions))


def make(description:dict) -> Figure:
    """
    This function takes a figure description and handles the contruction of a publication-grade
    figure in a convenient way.
    The description contains the arguments used for initialization of the objects in the
    :py:mod:`eos.figure` module.

    A figure can contain one or more plots, where each plot containts one or more items such
    as a data series from the evaluation of an observable or a kernel density estimate to
    visualize a set of statistical samples.

    The description can be either a dictionary (`figure description`) describing a single plot
    or a list of dictionaries describing a grid of plots.

    +---------------+----------------------------+
    | `figure description` (dict)                |
    +===============+============================+
    | `description` | `plot description` (dict)  |
    +---------------+----------------------------+
    | `items`       | list of `item description` |
    +---------------+----------------------------+

    +-------------+----------------------------------+
    | `plot description` (dict)                      |
    +=============+==================================+
    | `plot_type` | valid plot type (str)            |
    +-------------+----------------------------------+
    | ...         | one or several keyword arguments |
    +-------------+----------------------------------+

    +-------------+----------------------------------+
    | `item description` (dict)                      |
    +=============+==================================+
    | `item_type` | valid item type (str)            |
    +-------------+----------------------------------+
    | ...         | one or several keyword arguments |
    +-------------+----------------------------------+


    :param description: A figure description as documented
    :type description: dict or list of dict
    :returns: EOS figure
    :rtype: eos.figure.Figure
    """

    if type(description) == dict:
        # Make figure with a single plot
        return SingleFigure(description)

    elif type(description) == list:
        # Mad grid figure with several plots arranged on a grid
        raise RuntimeError("Making a grid figure is not yet implemented")

    else:
        raise RuntimeError("Argument 'description' must be either of type dict or list of dict")
