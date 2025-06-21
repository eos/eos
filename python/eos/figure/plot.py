# Copyright (c) 2023-2025 Danny van Dyk
# Copyright (c) 2023      Philip Lüghausen
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

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from eos.analysis_file_description import AnalysisFileContext
from eos.deserializable import Deserializable

from .item import ItemFactory, ItemColorCycler

import copy as _copy
import eos
import inspect
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import yaml as _yaml


@dataclass
class Legend(Deserializable):
    r"""Represents the legend properties in a plot."""

    position:str=field(default='best')

    def __post_init__(self):
        pass

    def draw(self, ax):
        ax.legend(loc=self.position)


@dataclass
class Grid(Deserializable):
    r"""Represents the grid properties in a plot."""

    visible:bool=field(default=False)
    axis:str=field(default='both')

    def __post_init__(self):
        if self.axis not in ['both', 'x', 'y']:
            raise ValueError(f'Unknown axis: {self.axis}')

    def draw(self, ax):
        ax.grid(visible=self.visible, axis=self.axis)


@dataclass
class XTicks(Deserializable):
    r""""Represents the x axis ticks properties in a plot."""

    minor:bool=field(default=True)
    visible:bool=field(default=True)

    def __post_init__(self):
        pass

    def draw(self, ax):
        if not self.visible:
            ax.xaxis.set_major_formatter(plt.NullFormatter())
        else:
            ax.xaxis.set_major_locator(matplotlib.ticker.AutoLocator())
            if self.minor:
                ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            ax.xaxis.set_tick_params(bottom=self.visible, top=self.visible)


@dataclass
class XAxis(Deserializable):
    r""""Represents the x axis properties in a plot."""

    label:str=field(default=None)
    range:tuple[float, float]=field(default=None)
    ticks:XTicks=field(default_factory=XTicks)
    unit:str=None

    def __post_init__(self):
        if self.range is not None and len(self.range) != 2:
            raise ValueError("Range must be a tuple of two values (min, max)")

        if self.range is not None:
            self.range = tuple(float(x) for x in self.range)

    def draw(self, ax):
        if self.label is not None and self.unit is not None:
            ax.set_xlabel(f'{self.label} [{self.unit}]')
        elif self.label is not None:
            ax.set_xlabel(self.label)

        if self.range is not None:
            ax.set_xlim(self.range)
        self.ticks.draw(ax)

    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        if 'ticks' in kwargs:
            _kwargs['ticks'] = XTicks.from_dict(**kwargs['ticks'])
        return Deserializable.make(cls, **_kwargs)


@dataclass
class YTicks(Deserializable):
    r""""Represents the y axis ticks properties in a plot."""

    minor:bool=field(default=True)
    visible:bool=field(default=True)

    def __post_init__(self):
        pass

    def draw(self, ax):
        if not self.visible:
            ax.yaxis.set_major_formatter(plt.NullFormatter())
        else:
            ax.yaxis.set_major_locator(matplotlib.ticker.AutoLocator())
            if self.minor:
                ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            ax.yaxis.set_tick_params(bottom=self.visible, top=self.visible)

@dataclass
class YAxis(Deserializable):
    r""""Represents the y axis properties in a plot."""

    label:str=None
    range:tuple[float, float]=None
    ticks:YTicks=field(default_factory=YTicks)
    unit:str=None

    def __post_init__(self):
        if self.range is not None and len(self.range) != 2:
            raise ValueError("Range must be a tuple of two values (min, max)")

        if self.range is not None:
            self.range = tuple(float(y) for y in self.range)

    def draw(self, ax):
        if self.label is not None and self.unit is not None:
            ax.set_ylabel(f'{self.label} [{self.unit}]')
        elif self.label is not None:
            ax.set_ylabel(self.label)

        if self.range is not None:
            ax.set_ylim(self.range)
        self.ticks.draw(ax)

    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        if 'ticks' in kwargs:
            _kwargs['ticks'] = YTicks.from_dict(**kwargs['ticks'])
        return Deserializable.make(cls, **_kwargs)


class Plot(ABC, Deserializable):
    r"""Base class for plots to be drawn into a figure."""

    @abstractmethod
    def prepare(self, context:AnalysisFileContext=None):
        "Prepare the plot for drawing"
        raise NotImplementedError

    @abstractmethod
    def draw(self, ax):
        "Draw the plot on the axes"
        raise NotImplementedError


@dataclass(kw_only=True)
class TwoDimensionalPlot(Plot):
    """Draws a 2D plot along a single set of axes."""

    aspect:float|None=field(default=None)
    grid:Grid=field(default_factory=Grid)
    items:list
    legend:Legend=field(default=None)
    title:str=None
    xaxis:XAxis=field(default_factory=XAxis)
    yaxis:YAxis=field(default_factory=YAxis)

    _api_doc = inspect.cleandoc("""
    Drawing a 2D Plot Along a Single Set of Axes
    --------------------------------------------

    This plot's type is ``2D``, which is the default plot type. It produces a two-dimensional plot
    along a single set of axes. The plot can contain multiple plot items.

    The following keys are mandatory:

        * ``items``: A list of items to be drawn in the plot.

    The following keys are optional:

        * ``aspect``: A float, which represents the aspect ratio of the plot.
        * ``grid``: An object of type :class:`eos.figure.Grid``, which contains the grid properties.
        * ``legend``: An object of type :class:`eos.figure.Legend``, which contains the legend properties.
        * ``title``: A string, which contains the title of the plot.
        * ``xaxis``: An object of type :class:`eos.figure.XAxis`, which contains the x axis properties.
        * ``yaxis``: An object of type :class:`eos.figure.YAxis`, which contains the y axis properties.

    """)

    def __post_init__(self):
        pass

    def prepare(self, context:AnalysisFileContext=None):
        context = AnalysisFileContext() if context is None else context
        for item in self.items:
            item.prepare(context=context)

    def draw(self, ax):
        "Draw the plot including all items"
        # Set title
        ax.set_title(self.title)

        # Remove default margin used by matplotlib
        ax.margins(0.0)

        # Draw grid
        self.grid.draw(ax)

        # Set aspect ratio
        if self.aspect is not None:
            ax.set_box_aspect(self.aspect)

        # Handle axes
        self.xaxis.draw(ax)
        self.yaxis.draw(ax)

        # Draw all items
        for item in self.items:
            item.draw(ax)

        # Draw legend
        if self.legend is not None:
            self.legend.draw(ax)


    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        if 'grid' in kwargs:
            _kwargs['grid'] = Grid.from_dict(**kwargs['grid'])
        if 'legend' in kwargs:
            _kwargs['legend'] = Legend.from_dict(**kwargs['legend'])
        if 'xaxis' in kwargs:
            _kwargs['xaxis'] = XAxis.from_dict(**kwargs['xaxis'])
        if 'yaxis' in kwargs:
            _kwargs['yaxis'] = YAxis.from_dict(**kwargs['yaxis'])
        _kwargs['items'] = [ItemFactory.from_dict(**i) for i in kwargs['items']]
        return Deserializable.make(cls, **_kwargs)


@dataclass
class EmptyPlot(Plot):
    """Draws an empty plot.

    Can be used in a grid as empty space instead of a plot."""

    _api_doc = inspect.cleandoc("""
    Drawing an Empty Plot
    ---------------------

    This plot's type is ``empty``. It produces an empty plot, which can be used in a grid as empty space.

    """)

    def __post_init__(self):
        pass

    def prepare(self, context:AnalysisFileContext=None):
        pass

    def draw(self, ax):
        ax.set_axis_off()


class PlotFactory:
    # Also build the documentation based on ordered registry
    # Initializer is well-defined for python version >= 3.6
    registry = {
        '2D':    TwoDimensionalPlot, # default
        'empty': EmptyPlot
    }

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

        ItemColorCycler.reset()

        return PlotFactory.registry[plot_type].from_dict(**kwargs)
