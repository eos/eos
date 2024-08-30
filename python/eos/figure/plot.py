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

    visible:bool=field(default=True)

    def __post_init__(self):
        pass

    def draw(self, ax):
        if not self.visible:
            ax.xaxis.set_major_formatter(plt.NullFormatter())
        ax.xaxis.set_tick_params(bottom=self.visible, top=self.visible)


@dataclass
class XAxis(Deserializable):
    r""""Represents the x axis properties in a plot."""

    label:str=field(default=None)
    range:tuple[float, float]=field(default=None)
    ticks:XTicks=field(default_factory=XTicks)

    def __post_init__(self):
        pass

    def draw(self, ax):
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

    visible:bool=field(default=True)

    def __post_init__(self):
        pass

    def draw(self, ax):
        if not self.visible:
            ax.yaxis.set_major_formatter(plt.NullFormatter())
        ax.yaxis.set_tick_params(left=self.visible, right=self.visible)


@dataclass
class YAxis(Deserializable):
    r""""Represents the y axis properties in a plot."""

    label:str=None
    range:tuple[float, float]=None
    ticks:YTicks=field(default_factory=YTicks)

    def __post_init__(self):
        pass

    def draw(self, ax):
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


@dataclass
class TwoDimensionalPlot(Plot):
    r"""Represents a 2D plot along a single set of axes.

    :param items: List of content items displayed in this plot
    :type items: list[dict] or iterable[dict]
    :param title: Title of the plot
    :type title: str or NoneType (optional)
    :param legend: Keyword arguments to pass to ``matplotlib.pyplot.legend``
    :type legend: dict
    """

    items:list
    title:str=None
    grid:Grid=field(default_factory=Grid)
    legend:Legend=field(default=None)
    xaxis:XAxis=field(default_factory=XAxis)
    yaxis:YAxis=field(default_factory=YAxis)
    aspect:float=field(default=1.0)

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
        ax.set_box_aspect(self.aspect)

        # Draw all items
        for item in self.items:
            item.draw(ax)

        # Handle axes
        self.xaxis.draw(ax)
        self.yaxis.draw(ax)

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
    r"""Represents and empty plot.

    Can be used in a grid as empty space instead of a plot.
    """

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
