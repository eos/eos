# Copyright (c) 2023-2025 Danny van Dyk
# Copyright (c) 2023      Philip Lueghausen
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

from .plot import Plot, PlotFactory
from .common import Watermark
from .data import DataFile

import os
import copy as _copy
import inspect
import matplotlib.pyplot as plt
import numpy as np
import yaml as _yaml

@dataclass(kw_only=True)
class Figure(ABC, Deserializable):
    r"""Base class for figures to be drawn using matplotlib."""
    name:str=field(default=None)

    @abstractmethod
    def draw(self, context:AnalysisFileContext=None, output:str|None=None):
        "Draw the figure"
        raise NotImplementedError


@dataclass(kw_only=True)
class SingleFigure(Figure):
    """Produces a figure with a single plot."""

    type:str=field(repr=False, init=False, default='single')

    plot:Deserializable
    size:tuple[float, float]=field(default=(6.4, 4.8))
    watermark:Watermark=field(default_factory=Watermark)

    _api_doc = inspect.cleandoc("""
    Producing a Figure with one Set of Plots
    ----------------------------------------

    This figure's type is ``single``, which is the default type of figure. It displays a single plot.

    The following keys are mandatory:

        * ``plot`` (:class:`Plot <eos.figure.plot.Plot>`) -- The plot to be drawn in the figure. This can be any of the
            :class:`Plot <eos.figure.plot.Plot>` descendants.

    The following keys are optional:

        * ``size`` (*tuple[float, float]*) -- The size of the figure in inches. Defaults to (6.4, 4.8).


    """)

    def __post_init__(self):
        self._figure, self._ax = plt.subplots(figsize=self.size)

    def draw(self, context:AnalysisFileContext=None, output:str|None=None):
        context = AnalysisFileContext() if context is None else context
        self.plot.prepare(context)
        self.plot.draw(self._ax)
        self.watermark.draw(self._ax)
        if output is not None:
            self._figure.savefig(output, bbox_inches='tight')

    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs['plot'] = PlotFactory.from_dict(**kwargs['plot'])
        if 'watermark' in kwargs:
            _kwargs['watermark'] = Watermark.from_dict(**kwargs['watermark'])
        return Deserializable.make(cls, **_kwargs)


@dataclass(kw_only=True)
class Inset(Deserializable):
    """Represents the inset properties for an `InsetFigure`."""

    plot:Deserializable
    position:tuple[float, float]
    size:tuple[float, float]

    def __post_init__(self):
        pass

    def prepare(self, context, ax):
        """Prepare the inset plot for drawing."""
        self._inset_ax = ax.inset_axes([self.position[0], self.position[1], self.size[0], self.size[1]])
        self.plot.prepare(context)

    def draw(self, ax):
        self.plot.draw(self._inset_ax)
        ax.indicate_inset_zoom(self._inset_ax, edgecolor="black")

    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs['plot'] = PlotFactory.from_dict(**kwargs['plot'])
        return Deserializable.make(cls, **_kwargs)


@dataclass(kw_only=True)
class InsetFigure(Figure):
    """Produces an inset figure with a main plot and a smaller inset plot."""

    type:str=field(repr=False, init=False, default='inset')

    plot:Deserializable
    inset:Deserializable
    size:tuple[float, float]=field(default=(6.4, 4.8))
    watermark:Watermark=field(default_factory=Watermark)

    _api_doc = inspect.cleandoc("""
    Producing a Figure with an Inset Plot
    -------------------------------------

    This figure's type is `inset`. It display a main plot covering the full area of the figure, and a smaller inset plot in one corner.

    The following keys are mandatory:

        * ``plot`` (:class:`Plot <eos.figure.plot.Plot>`) -- The main plot to be drawn in the figure. This can be any of the
            :class:`Plot <eos.figure.plot.Plot>` descendants.
        * ``inset`` (:class:`Inset <eos.figure.common.Inset>`) -- The inset plot to be drawn in the figure. This should be an instance of
            :class:`Inset <eos.figure.common.Inset>`.

    The following keys are optional:

        * ``size`` (*tuple[float, float]*) -- The size of the figure in inches. Defaults to (6.4, 4.8).
    """)

    def __post_init__(self):
        self._figure, self._ax = plt.subplots(figsize=self.size)

    def draw(self, context:AnalysisFileContext=None, output:str|None=None):
        context = AnalysisFileContext() if context is None else context
        self.plot.prepare(context)
        self.plot.draw(self._ax)
        self.watermark.draw(self._ax)
        self.inset.prepare(context, self._ax)
        self.inset.draw(self._ax)
        if output is not None:
            self._figure.savefig(output, bbox_inches='tight')

    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs['plot'] = PlotFactory.from_dict(**kwargs['plot'])
        _kwargs['inset'] = Inset.from_dict(**kwargs['inset'])
        if 'watermark' in kwargs:
            _kwargs['watermark'] = Watermark.from_dict(**kwargs['watermark'])
        return Deserializable.make(cls, **_kwargs)


@dataclass(kw_only=True)
class GridFigure(Figure):
    """Produces a figure with a configurable number of plots, arranged in a grid."""

    type:str=field(repr=False, init=False, default='grid')

    plots:list[Plot]
    shape:tuple[int, int]

    _api_doc = inspect.cleandoc("""
    Producing a Figure with a Grid of Plots
    ---------------------------------------

    This figure's type is ``grid``. It produces a grid of plots.

    The following keys are mandatory:

        * ``plots`` (*list* of :class:`Plot <eos.figure.plot.Plot>`) -- The list of plots to be drawn in the figure. Each plot can be any of the
            :class:`Plot <eos.figure.plot.Plot>` descendants. The list of plots are assigned to the grid positions in row-major order, i.e. the
            first plot is assigned to the first row and first column, the second plot to the first row and second column, etc.

        * ``shape`` (*tuple[int, int]*) -- The shape of the figure's grid, specifying the number of rows and columns in that order.


    """)
    def __post_init__(self):
        nrow, ncol = self.shape
        self._figure, axes = plt.subplots(nrow, ncol, figsize=(3.0 * nrow, 3.0 * ncol))
        self._axes = axes.flatten('C') # flatten to row-major style

    def draw(self, context:AnalysisFileContext=None, output:str|None=None):
        context = AnalysisFileContext() if context is None else context
        for idx, plot in enumerate(self.plots):
            plot.prepare()
            plot.draw(self._axes[idx])

        self._figure.tight_layout()
        if output is not None:
            self._figure.savefig(output, bbox_inches='tight')

    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs['plots'] = [PlotFactory.from_dict(**p) for p in kwargs['plots']]
        return Deserializable.make(cls, **_kwargs)


@dataclass(kw_only=True)
class CornerFigure(Figure):
    r"""Produces a corner figure, i.e., a figure with a triangular arrangement of 1D and 2D marginal PDFs.

    Distributions of the variables are shown on the diagonal, while correlations are plotted on the lower-left corner."""

    type:str=field(repr=False, init=False, default='corner')

    contents:list[DataFile]
    variables:list[str]=None

    _api_doc = inspect.cleandoc("""
    Producing a Corner Figure
    -------------------------

    This figure's type is ``corner``. It produces a corner figure, i.e., a figure with a triangular arrangement of 1D and 2D marginal PDFs.

    The following keys are mandatory:
        * ``contents`` (*list* of :class:`DataFile <eos.figure.data.DataFile>`) -- The list of data files to be drawn in the figure.
          Each data file should be a dictionary containing the path to the data file, a label and optionally a color.

    The following keys are optional:
        * ``variables`` (*list[str]*) -- The list of variable names to be considered. Defaults to None, in which case all variables contained in the first data file are shown.

    """)

    def __post_init__(self):
        if not self.contents:
            raise ValueError("Contents must include at least one item to be plotted.")

    def prepare(self, context:AnalysisFileContext=None, output:str|None=None):
        context = AnalysisFileContext() if context is None else context
        for content in self.contents:
            content.prepare(context=context)

        if self.variables:
            self._variables = self.variables
        else:
            self._variables = list(self.contents[0].variables)

        # determine useful ranges empirically
        absmin, absmax = np.array([+np.inf] * len(self._variables)), np.array([-np.inf] * len(self._variables))
        for content in self.contents:
            indices = [content._datafile.lookup_table[v] for v in self._variables]
            cmin, cmax = content.empirical_range
            for idx, cidx in enumerate(indices):
                absmin[idx] = cmin[cidx] if cmin[cidx] < absmin[idx] else absmin[idx]
                absmax[idx] = cmax[cidx] if cmax[cidx] > absmax[idx] else absmax[idx]

        for idx, v in enumerate(self._variables):
            print(f"variable {v}: [{absmin[idx]}, {absmax[idx]}]")

        # Check that the variables of the data files match
        for content in self.contents:
            unknown_variables = set(self._variables) - set(content.variables)
            if len(unknown_variables) > 0:
                raise ValueError(f"Unknown variables requested from data file '{content.path}': {list(unknown_variables)}")

        self._labels = self.contents[0].labels(self._variables)

        plots = []
        size = len(self._variables)

        for i in range(size):
            for j in range(size):

                if i < j:
                    plots.append(PlotFactory.from_dict(**{
                        'type': 'empty'
                    }))

                elif i == j:
                    plots.append(PlotFactory.from_dict(**{
                        'xaxis': {
                            'ticks': { 'visible': True },
                            'label': self._labels[i],
                            'range': [ absmin[i], absmax[i] ]
                        }
                        if (i == size - 1) else
                        {
                            'ticks': { 'visible': False }
                        },
                        'yaxis': {
                            'ticks': { 'visible': False },
                        },
                        'grid': { 'visible': True, 'axis': 'x' },
                        'items': [
                            {
                                'type': 'kde1D', 'label': content.label,
                                'datafile': context.data_path(content.path),
                                'variable': self._variables[i],
                                'color': content.color,
                            }
                        for content in self.contents]
                    }))

                else:
                    plots.append(PlotFactory.from_dict(**{
                        'xaxis': {
                            'ticks': { 'visible': True },
                            'label': self._labels[j],
                            'range': [ absmin[j], absmax[j] ]
                        }
                        if (i == size - 1) else
                        {
                            'ticks': { 'visible': False }
                        },
                        'yaxis': {
                            'ticks': { 'visible': True },
                            'label': self._labels[i],
                            'range': [ absmin[i], absmax[i] ]
                        }
                        if (j == 0) else
                        {
                            'ticks': { 'visible': False }
                        },
                        'grid': { 'visible': True},
                        'items': [
                            {
                                'type': 'kde2D', 'label': content.label ,
                                'datafile': context.data_path(content.path),
                                'variables': [self._variables[j], self._variables[i]],
                                'color': content.color,
                            }
                        for content in self.contents]
                    }))

        self._figure = GridFigure(shape=(size, size), plots=plots)


    def draw(self, context:AnalysisFileContext=None, output:str|None=None):
        context = AnalysisFileContext() if context is None else context
        if not hasattr(self, '_figure'):
            self.prepare(context=context, output=output)

        self._figure.draw(context=context, output=output)


    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        if 'contents' in kwargs:
            _kwargs['contents'] = [DataFile.from_dict(**c) for c in kwargs['contents']]
        return Deserializable.make(cls, **_kwargs)


class FigureFactory:
    # Also build the documentation based on ordered registry
    # Initializer is well-defined for python version >= 3.6
    registry = {
        'single': SingleFigure, # default
        'inset':  InsetFigure,
        'grid':   GridFigure,
        'corner': CornerFigure,
    }

    @staticmethod
    def from_yaml(yaml_data:str):
        kwargs = _yaml.safe_load(yaml_data)
        return FigureFactory.from_dict(**kwargs)

    @staticmethod
    def from_dict(**kwargs):
        r"""Factory method to create a plot from a dictionary

        This method takes a dictionary description of a figure and uses matplotlib to render that
        figure conveniently and at a publication-grade level.

        The description contains the arguments used for initializing objects of classes within the
        :py:mod:`eos.figure` module.

        A figure can contain one or more individual plots. Each plot containts one or more plottable
        items, such as a data series from the evaluation of an observable or a kernel density estimate to
        visualize a set of statistical samples.

        :param description: A figure description as required by any of the descendant classes of :py:`eos.figure.Figure`.
        :type description: dict
        :returns: A descendant of the :py:class:`eos.figure.Figure` class
        :rtype: :py:class:`eos.figure.Figure`
        """

        if 'type' not in kwargs:
            figure_type = 'single'
        else:
            figure_type = kwargs.pop('type')

        if figure_type not in FigureFactory.registry:
            raise ValueError(f'Unknown figure type: {figure_type}')

        return FigureFactory.registry[figure_type].from_dict(**kwargs)
