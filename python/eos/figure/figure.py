# Copyright (c) 2023-2024 Danny van Dyk
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

from collections import OrderedDict
from dataclasses import dataclass, field
from eos.deserializable import Deserializable

from .plot import Plot, PlotFactory

import eos
import os
import copy as _copy
import matplotlib.pyplot as plt
import numpy as np
import yaml as _yaml

class Figure(Deserializable):
    r"""Base class for figures to be drawn using matplotlib."""
    name:str=field(default=None)

    def draw(self):
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

    def draw(self):
        self.plot.prepare()
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

    def draw(self):
        for idx, plot in enumerate(self.plots):
            plot.prepare()
            plot.draw(self._axes[idx])

        self._figure.tight_layout()

    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs['plots'] = [PlotFactory.from_dict(**p) for p in kwargs['plots']]
        return Deserializable.make(cls, **_kwargs)


@dataclass(kw_only=True)
class DataFileContent(Deserializable):
    r""" Collects the relevant information for a data file to be plotted.

    :param datafilepath: Path to the eos.Prediction or eos.ImportanceSample file to be used.
    :type datafilepath: str
    :param label: Name of the data file that will appear in the plot's legend.
    :type label: str
    :param color: Color to be used in the plot.
    :type color: str (optional)
    """
    datafilepath:str
    label:str
    color:str=None

    def validate(self):
        os.path.exists(self.datafilepath) or eos.error(f"Data file '{self.datafilepath}' does not exist")

    def datafile(self):
        self.validate()
        abspath = os.path.abspath(self.datafilepath)
        name = os.path.split(abspath)[-1]

        if name == 'samples':
            # The datafile refers to importance (parameter) samples
            return eos.data.ImportanceSamples(abspath)
        elif name.startswith('pred-'):
            # The datafile refers to observable predictions
            return eos.data.Prediction(abspath)
        else:
            eos.error(f"Data file '{self.datafilepath}' has an unsupported format")
            raise NotImplementedError

    def variables(self):
        return self.datafile().lookup_table.keys()

    def variablelabels(self, variablelist):
        self.validate()
        abspath = os.path.abspath(self.datafilepath)
        name = os.path.split(abspath)[-1]

        if name == 'samples':
            # The datafile refers to importance (parameter) samples
            df = eos.data.ImportanceSamples(abspath)
            p = eos.Parameters()
            return [p[dist].latex() for dist in variablelist]
        elif name.startswith('pred-'):
            # The datafile refers to observable predictions
            df = eos.data.Prediction(abspath)
            o = eos.Observables()
            return ["$" + o[dist].latex() + "$" for dist in variablelist]
        else:
            eos.error(f"Data file '{self.datafilepath}' has an unsupported format")
            raise NotImplementedError


@dataclass(kw_only=True)
class CornerFigure(Figure):
    r"""Produces a grid figure representing the variables described by one or several datafiles.
    Distributions of the variables are shown on the diagonal, while correlations are plotted on the lower-left corner.

    :param contents: List of contents to be drawn. Contents should be dictionnaries containing the path to the data file, a label and optionnaly a color.
    :type contents: list[DataFileContent] or iterable[DataFileContent]
    :param variables: List of the names of the variables to be considered. Defaults to None, in this case all variables are shown.
    :type variables: list[str] (optional)
    """
    type:str=field(repr=False, init=False, default='corner')
    contents:list[DataFileContent]
    variables:list[str]=None

    def __post_init__(self):
        if not self.contents:
            raise ValueError("At least one content should be plotted.")


    def draw(self):
        # Check that the variables of the data files match
        for content in self.contents:
            if self.variables:
                if not set(self.variables).issubset(content.variables()):
                    raise ValueError(f"The provided variables are not contained in the data file {content.datafilepath}")
            elif content.variables() != self.contents[0].variables():
                raise ValueError("variables defaulted to None but the provided datafiles contain different variables")

        if self.variables:
            self._variables = self.variables
        else:
            self._variables = list(self.contents[0].variables())

        self._labels = self.contents[0].variablelabels(self._variables)

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
                        'xaxis': { 'ticks': { 'visible': True }, 'label': self._labels[i] } if (i == size - 1) else { 'ticks': { 'visible': False } },
                        'yaxis': { 'ticks': { 'visible': False } },
                        'grid': { 'visible': True, 'axis': 'x' },
                        'items': [
                            {
                                'type': 'kde1D', 'label': content.label,
                                'datafile': content.datafilepath,
                                'variable': self._variables[i],
                                'color': content.color,
                            }
                        for content in self.contents]
                    }))

                else:
                    plots.append(PlotFactory.from_dict(**{
                        'xaxis': { 'ticks': { 'visible': True }, 'label': self._labels[j] } if (i == size - 1) else { 'ticks': { 'visible': False } },
                        'yaxis': { 'ticks': { 'visible': True }, 'label': self._labels[i] } if (j == 0) else { 'ticks': { 'visible': False } },
                        'grid': { 'visible': True},
                        'items': [
                            {
                                'type': 'kde2D', 'label': content.label ,
                                'datafile': content.datafilepath,
                                'variables': [self._variables[j], self._variables[i]],
                                'color': content.color,
                            }
                        for content in self.contents]
                    }))

        GridFigure(shape=[size, size], plots=plots).draw()

    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        if 'contents' in kwargs:
            _kwargs['contents'] = [DataFileContent.from_dict(**c) for c in kwargs['contents']]
        return Deserializable.make(cls, **_kwargs)


class FigureFactory:
    # Also build the documentation based on ordered registry
    # Initializer is well-defined for python version >= 3.6
    registry = OrderedDict({
        'single': SingleFigure, # default
        'grid':   GridFigure,
        'corner': CornerFigure,
    })

    @staticmethod
    def from_yaml(yaml_data:str):
        kwargs = _yaml.safe_load(yaml_data)
        return FigureFactory.from_dict(**kwargs)

    @staticmethod
    def from_dict(**kwargs):
        r"""Factory method to create a plot from a dictionary

        This method takes a dictionary description of a figure and uses matplotlib to render that
        figure conveniently and ata publication-grade level.

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
