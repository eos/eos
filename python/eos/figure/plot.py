# Copyright (c) 2023-2024 Danny van Dyk
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

from collections import OrderedDict
from dataclasses import dataclass, field
from eos.deserializable import Deserializable

from .item import ItemFactory

import copy as _copy
import matplotlib.pyplot as plt
import numpy as np
import yaml as _yaml

@dataclass
class TwoDimensionalPlot(Deserializable):
    r"""Represents a 1D plot along a single set of axes.

    :param items: List of content items displayed in this plot
    :type items: list[dict] or iterable[dict]
    :param title: Title of the plot
    :type title: str or NoneType (optional)
    :param legend: Keyword arguments to pass to ``matplotlib.pyplot.legend``
    :type legend: dict
    """

    items:list
    title:str=None
    legend:dict=None,

    def __post_init__(self):
        pass
        #for item in items:
        #    if 'type' not in item:
        #        raise ValueError(f"Content item { item } does not contain element 'type'")

        #self.content_items = []
        #for item in items:
        #    self.content_items.append(ItemFactory.make(ax, item['type'], **item['description']))

    def prepare(self):
        for item in self.items:
            item.prepare()

    def draw(self, ax):
        "Draw the plot including all items"
        # A child of Plot class can add items to ax before this method should be called.
        # Child classes must call super().draw() at the **end** of their draw() method.

        # Set title
        ax.set_title(self.title)

        # Remove default margin used by matplotlib
        ax.margins(0.0)

        # Draw all items
        for item in self.items:
            item.draw(ax)

    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs['items'] = [ItemFactory.from_dict(**i) for i in kwargs['items']]
        return Deserializable.make(cls, **_kwargs)


#class TwoDimensionalPlot(Plot):
#    r"""Represents a 2D plot with x and y axes.
#
#    :param items: List of content items displayed in this plot
#    :type items: list[dict] or iterable[dict]
#    :param x_range: Range of the x axis (automatic if value is `[None, None]`)
#    :type x_range: tuple of float (optional)
#    :param y_range: Range of the y axis (automatic if value is `[None, None]`)
#    :type y_range: tuple of float (optional)
#    :param x_label: Label of the x axis
#    :type x_label: str (optional)
#    :param y_label: Label of the y axis
#    :type y_label: str (optional)
#    :param \**kwargs: Additional keyword arguments to pass to the base class :py:class:`eos.figure.plot.Plot`
#    """
#    def __init__(self,
#            ax,
#            items:list,
#            x_range:tuple[float, float]=(None, None),
#            y_range:tuple[float, float]=(None, None),
#            x_label:str=None,
#            y_label:str=None,
#            **kwargs
#        ):
#        super().__init__(ax, items, **kwargs)
#
#        self.x_range = x_range
#        self.y_range = y_range
#
#        self.ax.set_xlabel(x_label)
#        self.ax.set_ylabel(y_label)
#
#        # Ranges: if not set, will be handled by matplotlib
#        if not self.x_range == (None, None):
#            self.ax.set_xlim(self.x_range)
#        if not self.y_range == (None, None):
#            self.ax.set_ylim(self.y_range)

@dataclass
class EmptyPlot(Deserializable):
    r"""Represents and empty plot.

    Can be used in a grid as empty space instead of a plot.
    """

    def __post_init__(self):
        pass

    def prepare(self):
        pass

    def draw(self, ax):
        ax.set_axis_off()


class PlotFactory:
    # Also build the documentation based on ordered registry
    # Initializer is well-defined for python version >= 3.6
    registry = OrderedDict({
        '2D':    TwoDimensionalPlot, # default
        'empty': EmptyPlot
    })

    @staticmethod
    def from_yaml(yaml_data:str):
        kwargs = _yaml.safe_load(yaml_data)
        return PlotFactory.from_dict(**kwargs)

    @staticmethod
    def from_dict(**kwargs):
        "Factory method to create a plot from a dictionary"

        if 'type' not in kwargs:
            plot_type = '2D'
        else:
            plot_type = kwargs.pop('type')

        if plot_type not in PlotFactory.registry:
            raise ValueError(f'Unknown plot type: {plot_type}')

        return PlotFactory.registry[plot_type].from_dict(**kwargs)

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
