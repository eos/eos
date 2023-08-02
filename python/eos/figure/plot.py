# Copyright (c) 2023 Danny van Dyk
# Copyright (c) 2023 Philip Lüghausen
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

import eos
from .item import ItemFactory

from collections import OrderedDict
import matplotlib.pyplot as plt
import numpy as np

class Plot:
    """Represents a plot along a single set of axes.

    :param ax: Set of axes
    :type ax: matplotlib.pyplot.Axes or similar
    :param items: List of content items displayed in this plot
    :type items: list[dict] or iterable[dict]
    :param title: Title of the plot
    :type title: str or NoneType (optional)
    :param legend: Keyword arguments to pass to ``matplotlib.pyplot.legend``
    :type legend: dict
    """
    def __init__(self,
            ax,
            items:list,
            title=None,
            legend:dict=None,
        ):
        self.ax = ax
        self.ax.set_title(title)
        self.legend_kwargs = legend

        for item in items:
            if 'type' not in item:
                raise ValueError(f"Content item { item } does not contain element 'type'")

        self.content_items = []
        for item in items:
            self.content_items.append(ItemFactory.make(ax, item['type'], **item['description']))


    def draw(self):
        "Draw the plot including all items"
        # A child of Plot class can add items to ax before this method should be called.
        # Child classes must call super().draw() at the **end** of their draw() method.

        # Remove default margin used by matplotlib
        self.ax.margins(0.0)

        # Draw all items
        for item in self.content_items:
            item.draw(self.ax)

        self.ax.legend(**self.legend_kwargs)


class TwoDimPlot(Plot):
    r"""Represents a 2d plot with x and y axes.

    :param ax: Set of axes
    :type ax: matplotlib.pyplot.Axes or similar
    :param items: List of content items displayed in this plot
    :type items: list[dict] or iterable[dict]
    :param x_range: Range of the x axis (automatic if value is `[None, None]`)
    :type x_range: tuple of float (optional)
    :param y_range: Range of the y axis (automatic if value is `[None, None]`)
    :type y_range: tuple of float (optional)
    :param x_label: Label of the x axis
    :type x_label: str (optional)
    :param y_label: Label of the y axis
    :type y_label: str (optional)
    :param \**kwargs: Additional keyword arguments to pass to the base class :py:class:`eos.figure.plot.Plot`
    """
    def __init__(self,
            ax,
            items:list,
            x_range:tuple[float, float]=(None, None),
            y_range:tuple[float, float]=(None, None),
            x_label:str=None,
            y_label:str=None,
            **kwargs
        ):
        super().__init__(ax, items, **kwargs)

        self.x_range = x_range
        self.y_range = y_range

        self.ax.set_xlabel(x_label)
        self.ax.set_ylabel(y_label)

        # Ranges: if not set, will be handled by matplotlib
        if not self.x_range == (None, None):
            self.ax.set_xlim(self.x_range)
        if not self.y_range == (None, None):
            self.ax.set_ylim(self.y_range)


class EmptyPlot(Plot):
    "Can be used in a grid as empty space instead of a plot."

    def __init__(self, ax, contents:list, **kwargs):
        super().__init__(ax, contents, **kwargs)
        ax.set_axis_off()

    def draw(self):
        "Draw nothing"
        pass


class PlotFactory:
    # Also build the documentation based on ordered registry
    # Initializer is well-defined for python version >= 3.6
    registry = OrderedDict({
        '2D-plot': TwoDimPlot,
        'empty': EmptyPlot
    })

    @staticmethod
    def sanitize(kwargs):
        result = kwargs

        # Issue here: non-descriptive error messages about missing arguments.
        replacements = [
            ('plot_type', 'type'),
        ]

        for old_key, new_key in replacements:
            if old_key not in result:
                continue

            result[new_key] = result[old_key]
            del result[old_key]

        # print(f'Sanitize result = {result}')
        return result


    @staticmethod
    def make(ax, contents:list, type:str='2D-plot', **kwargs):
        if type not in PlotFactory.registry:
            raise ValueError(f'Unknown plot type: { type }')

        return PlotFactory.registry[type](ax, contents, **kwargs)
