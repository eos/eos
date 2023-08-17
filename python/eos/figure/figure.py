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
    def __init__(self,
                 figure,
                 axes:list,
                 plots:list,
                 title:str=None):
        self._figure = figure
        self._axes   = axes
        self._plots  = []

        if len(plots) == 0:
            raise ValueError(f"Argument 'plots' is empty list")

        for ax, p in zip(axes, plots):
            if p == {} or p is None:
                self._plots.append(PlotFactory.make(ax, [], plot_type='empty'))
                continue

            if 'description' not in p:
                raise ValueError(f"Plot does not contain element 'description'")
            if 'items' not in p:
                raise ValueError(f"Plot does not contain element 'items'")

            items = p['items']
            description = PlotFactory.sanitize(p['description'])

            self._plots.append(PlotFactory.make(ax, items, **description))

        if title is not None:
            self._figure.suptitle(title)


    def draw(self):
        for plot in self._plots:
            plot.draw()


class SingleFigure(Figure):
    """Represent a figure with a single set of axes.

    Use the value of 'plot' as a set of arguments to construct a Figure.

    :param plot: The plot's data, as a dictionary containing the 'description' and the 'contents'.
    :type plot: dict
    :param title: The title of the figure. (optional)
    :type title: str
    """

    def __init__(self,
                 plot=None,
                 title:str=None):
        figure, axes = plt.subplots()
        if title is not None:
            axes.set_title(title)

        super().__init__(figure, [axes], [plot], title=title)

    @property
    def plot(self):
        return self._plots[0]


    @plot.setter
    def plot(self, value):
        self._plots[0] = value


class GridFigure(Figure):
    """Represent a figure with multiple sets of axes arranged in a grid.

    :param plots: The list of plot data sets, as a list of dictionaries each containing the ``description`` and the ``contents`` for a single plot.
    :type plots: list[dict] or iterable[dict]
    :param nrows: The number of rows in the grid.
    :type nrows: int
    :param ncols: The number of columns in the grid.
    :type ncols: int
    :param title: The title of the figure. (optional)
    :type title: str
    """

    def __init__(self,
                 plots:list,
                 nrows:int,
                 ncols:int,
                 title:str=None):
        if len(plots) != nrows * ncols:
            raise ValueError(f"Number of plots ({len(plots)}) does not match number of axes ({nrows*ncols})")

        figure, axes = plt.subplots(nrows, ncols)
        axes = axes.flatten('C') # row-major style
        super().__init__(figure, [ax for ax in axes], [p for p in plots], title=title)


def make(**kwargs) -> Figure:
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

    if 'items' in kwargs:
        return SingleFigure(**kwargs)
    elif 'plots' in kwargs:
        return GridFigure(**kwargs)
    else:
        raise RuntimeError("Figure description does not match any known format")
