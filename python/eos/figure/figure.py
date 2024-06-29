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

import copy as _copy
import matplotlib.pyplot as plt
import numpy as np
import yaml as _yaml

@dataclass(kw_only=True)
class Figure(ABC, Deserializable):
    r"""Base class for figures to be drawn using matplotlib."""
    name:str=field(default=None)

    @abstractmethod
    def draw(self, context:AnalysisFileContext=None):
        "Draw the figure"
        raise NotImplementedError


@dataclass(kw_only=True)
class SingleFigure(Figure):
    r"""Represents a figure with a single set of axes.

    :param plot: A representation of the various plot types created by :py:class:`eos.figure.PlotFactory`.
    :type plot: :py:class:`eos.figure.Plot`
    """
    type:str=field(repr=False, init=False, default='single')
    plot:Deserializable

    def __post_init__(self):
        self._figure, self._ax = plt.subplots(figsize=(6, 4))

    def draw(self, context:AnalysisFileContext=None):
        context = AnalysisFileContext() if context is None else context
        self.plot.prepare(context)
        self.plot.draw(self._ax)

    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs['plot'] = PlotFactory.from_dict(**kwargs['plot'])
        return Deserializable.make(cls, **_kwargs)


@dataclass(kw_only=True)
class GridFigure(Figure):
    r"""Represents a figure with several axes arranged in a grid.

    :param shape: The shape of the grid, i.e. the number of rows and columns.
    :type shape: tuple[int, int]
    :param plots: The flattened list of plots (row-major style).
    :type plots: list[Plot] or iterable[Plot]
    """
    type:str=field(repr=False, init=False, default='grid')
    shape:tuple[int, int]
    plots:list[Plot]

    def __post_init__(self):
        nrow, ncol = self.shape
        self._figure, axes = plt.subplots(nrow, ncol, figsize=(3.0 * nrow, 3.0 * ncol))
        self._axes = axes.flatten('C') # flatten to row-major style

    def draw(self, context:AnalysisFileContext=None):
        context = AnalysisFileContext() if context is None else context
        for idx, plot in enumerate(self.plots):
            plot.prepare()
            plot.draw(self._axes[idx])

        self._figure.tight_layout()

    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs['plots'] = [PlotFactory.from_dict(**p) for p in kwargs['plots']]
        return Deserializable.make(cls, **_kwargs)


class FigureFactory:
    # Also build the documentation based on ordered registry
    # Initializer is well-defined for python version >= 3.6
    registry = {
        'single': SingleFigure, # default
        'grid':   GridFigure,
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
