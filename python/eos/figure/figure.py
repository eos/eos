# Copyright (c) 2026      Carolina Bolognani
# Copyright (c) 2023-2026 Danny van Dyk
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
from eos.analysis_file_context import AnalysisFileContext
from eos.deserializable import Deserializable
from eos.diagnostic import Diagnostic, Severity
import eos.data

from .plot import Plot, PlotFactory
from .common import Watermark
from .data import DataFile

import copy as _copy
import inspect
import matplotlib.pyplot as plt
import numpy as np
import yaml as _yaml
import re

@dataclass(kw_only=True)
class Figure(ABC, Deserializable):
    r"""Base class for figures to be drawn using matplotlib."""
    name:str=field(default=None)

    def validate_semantics(self, context):
        yield from ()

    def validate_structure(self):
        if self.name is None:
            return
        if '/' in self.name:
            yield Diagnostic(('name',), Severity.ERROR, f"Invalid character '/' in figure name '{self.name}'")
        if any(character.isspace() for character in self.name):
            yield Diagnostic(('name',), Severity.ERROR, f"Invalid whitespace in figure name '{self.name}'")

    @abstractmethod
    def draw(self, context:AnalysisFileContext=None):
        """Draw the figure.

        :param context: The analysis file context, which contains the paths to the data files and other relevant information.
        :type context: :class:`AnalysisFileContext <eos.analysis_file_context.AnalysisFileContext>`
        """

        raise NotImplementedError

    def save(self, output:str|list[str]):
        """Save the figure to one or more output files.

        The figure must have been drawn by a call to :meth:`draw <eos.figure.Figure.draw>` beforehand;
        saving an undrawn figure yields an empty figure. Drawing the same figure more than once stacks
        its elements on top of each other, which distorts the colors of elements with opacity < 1;
        draw once and pass all output files to a single call to this method instead.

        :param output: The path(s) to the output file(s) where the figure is saved. The file format of each
            output file is determined by its file name extension.
        :type output: str | list[str]
        """

        for out in [output] if isinstance(output, str) else output:
            self._figure.savefig(out, bbox_inches='tight')


@dataclass(kw_only=True)
class SingleFigure(Figure):
    """Produces a figure with a single plot.

    :param plot: The plot to be drawn in the figure. This can be any of the :class:`Plot <eos.figure.plot.Plot>` descendants.
    :type plot: :class:`Plot <eos.figure.plot.Plot>`
    :param size: The size of the figure in inches. Defaults to (6.4, 4.8).
    :type size: tuple[float, float]
    :param watermark: The optional specification where and how to draw the EOS watermark. For the default specification, see :class:`Watermark <eos.figure.Watermark>`.
    :type watermark: :class:`Watermark <eos.figure.Watermark>`
    """

    type:str=field(repr=False, init=False, default='single')

    plot:Plot
    size:tuple[float, float]=field(default=(6.4, 4.8))
    watermark:Watermark=field(default_factory=Watermark)

    _api_doc = inspect.cleandoc("""
    Producing a Figure with a single Plot
    -------------------------------------

    This figure's type is ``single``, which is the default type of figure. It displays a single plot.

    The following keys are mandatory:

        * ``plot`` (:class:`Plot <eos.figure.plot.Plot>`) -- The plot to be drawn in the figure. This can be any of the
            :class:`Plot <eos.figure.plot.Plot>` descendants.

    The following keys are optional:

        * ``size`` (*tuple[float, float]*) -- The size of the figure in inches. Defaults to (6.4, 4.8).


    """)

    def __post_init__(self):
        self._figure, self._ax = plt.subplots(figsize=self.size)

    def validate_semantics(self, context):
        yield from (
            diagnostic.prefixed('plot')
            for diagnostic in self.plot.validate_semantics(context)
        )

    def draw(self, context:AnalysisFileContext=None):
        """Draw the single-plot figure.

        Prepares and draws the figure's plot and adds the watermark. Use
        :meth:`save <eos.figure.Figure.save>` to store the drawn figure to one or more output files.

        :param context: The analysis file context, which contains the paths to the data files and other relevant information.
        :type context: :class:`AnalysisFileContext <eos.analysis_file_context.AnalysisFileContext>`
        """
        context = AnalysisFileContext() if context is None else context
        self.plot.prepare(context)
        self.plot.draw(self._ax)
        self.watermark.draw(self._ax)

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`SingleFigure` from its keyword description.

        Recursively deserializes the nested ``plot`` description via
        :class:`PlotFactory <eos.figure.plot.PlotFactory>` and the optional ``watermark`` description.

        :param kwargs: The figure description. Must contain a ``plot`` key.
        :returns: The instantiated figure.
        :rtype: SingleFigure
        """
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs['plot'] = PlotFactory.from_dict(**_kwargs['plot'])
        if 'watermark' in _kwargs:
            _kwargs['watermark'] = Watermark.from_dict(**_kwargs['watermark'])
        return Deserializable.make(cls, **_kwargs)


@dataclass(kw_only=True)
class Inset(Deserializable):
    """Represents the inset properties for an `InsetFigure`.

    :param plot: The inset plot to be drawn in the figure. This can be any of the :class:`Plot <eos.figure.plot.Plot>` descendants.
    :type plot: :class:`Plot <eos.figure.plot.Plot>`
    :param position: The position of the bottom left corner of the inset plot in the figure, specified as a tuple of (x, y) coordinates in normalized figure coordinates (0 to 1).
    :type position: tuple[float, float]
    :param size: The size of the inset plot, specified as a tuple of (width, height) in normalized figure coordinates (0 to 1).
    :type size: tuple[float, float]
    """

    plot:Plot
    position:tuple[float, float]
    size:tuple[float, float]

    def __post_init__(self):
        pass

    def validate_semantics(self, context):
        yield from (
            diagnostic.prefixed('plot')
            for diagnostic in self.plot.validate_semantics(context)
        )

    def prepare(self, context, ax):
        """Prepare the inset plot for drawing.

        Creates the inset axes at the configured position and size within the parent axes and prepares
        the inset's plot.

        :param context: The analysis file context forwarded to the inset's plot.
        :type context: AnalysisFileContext | None
        :param ax: The parent matplotlib axes into which the inset axes are inserted.
        :type ax: matplotlib.axes.Axes
        """
        self._inset_ax = ax.inset_axes([self.position[0], self.position[1], self.size[0], self.size[1]])
        self.plot.prepare(context)

    def draw(self, ax):
        """Draw the inset plot and indicate its zoom region on the parent axes.

        :param ax: The parent matplotlib axes on which the inset's zoom region is indicated.
        :type ax: matplotlib.axes.Axes
        """
        self.plot.draw(self._inset_ax)
        ax.indicate_inset_zoom(self._inset_ax, edgecolor="black")

    @classmethod
    def from_dict(cls, **kwargs):
        """Create an :class:`Inset` from its keyword description.

        Recursively deserializes the nested ``plot`` description via
        :class:`PlotFactory <eos.figure.plot.PlotFactory>`.

        :param kwargs: The inset description. Must contain a ``plot`` key.
        :returns: The instantiated inset.
        :rtype: Inset
        """
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs['plot'] = PlotFactory.from_dict(**_kwargs['plot'])
        return Deserializable.make(cls, **_kwargs)


@dataclass(kw_only=True)
class InsetFigure(Figure):
    """Produces an inset figure with a main plot and a smaller inset plot.

    :param plot: The main plot to be drawn in the figure. This can be any of the :class:`Plot <eos.figure.Plot>` descendants.
    :type plot: :class:`Plot <eos.figure.Plot>`
    :param inset: The inset plot to be drawn in the figure. This should be an instance of :class:`Inset <eos.figure.Inset>`.
    :type inset: :class:`Inset <eos.figure.Inset>`
    :param size: The size of the figure in inches. Defaults to (6.4, 4.8).
    :type size: tuple[float, float]
    :param watermark: The optional specification where and how to draw the EOS watermark.
    :type watermark: :class:`Watermark <eos.figure.Watermark>`
    """

    type:str=field(repr=False, init=False, default='inset')

    plot:Plot
    inset:Plot
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

    def validate_semantics(self, context):
        yield from (
            diagnostic.prefixed('plot')
            for diagnostic in self.plot.validate_semantics(context)
        )
        yield from (
            diagnostic.prefixed('inset')
            for diagnostic in self.inset.validate_semantics(context)
        )

    def draw(self, context:AnalysisFileContext=None):
        """Draw the inset figure.

        Use :meth:`save <eos.figure.Figure.save>` to store the drawn figure to one or more output files.

        :param context: The analysis file context, which contains the paths to the data files and other relevant information.
        :type context: :class:`AnalysisFileContext <eos.analysis_file_context.AnalysisFileContext>`
        """
        context = AnalysisFileContext() if context is None else context
        self.plot.prepare(context)
        self.plot.draw(self._ax)
        self.watermark.draw(self._ax)
        self.inset.prepare(context, self._ax)
        self.inset.draw(self._ax)

    @classmethod
    def from_dict(cls, **kwargs):
        """Create an :class:`InsetFigure` from its keyword description.

        Recursively deserializes the nested ``plot`` and ``inset`` descriptions, as well as the optional
        ``watermark`` description.

        :param kwargs: The figure description. Must contain ``plot`` and ``inset`` keys.
        :returns: The instantiated figure.
        :rtype: InsetFigure
        """
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs['plot'] = PlotFactory.from_dict(**_kwargs['plot'])
        _kwargs['inset'] = Inset.from_dict(**_kwargs['inset'])
        if 'watermark' in _kwargs:
            _kwargs['watermark'] = Watermark.from_dict(**_kwargs['watermark'])
        return Deserializable.make(cls, **_kwargs)


@dataclass(kw_only=True)
class GridFigure(Figure):
    """Produces a figure with a configurable number of plots, arranged in a grid.

    The list of plots are assigned to the grid positions in row-major order, i.e. the first plot is assigned to the first row
    and first column, the second plot to the first row and second column, and so on.

    :param padding: The tuple of horizontal and vertical padding between the plots in the grid, specified as a tuple of fractions of the average plot size. Defaults to (0.2, 0.2).
    :type padding: tuple[float, float]
    :param plots: The list of :class:`Plot <eos.figure.Plot>` objects to be drawn in the figure.
    :type plots: list[:class:`Plot <eos.figure.Plot>`]
    :param shape: The tuple specifying the shape of the figure's grid, specifying the number of rows and columns in that order.
    :type shape: tuple[int, int]
    :param size: The size of the figure in inches. Defaults to (3.0 * ncol, 3.0 * nrow).
    :type size: tuple[float, float]
    :param watermark_plot: The plot that carries the watermark, as a flattened (row-major) index or a 2D ``(row, col)`` address. If None, every plot is stamped.
    :type watermark_plot: int | tuple[int, int] | None
    :param tight_layout: If True (default), the grid spec is laid out with ``tight_layout``. Set to False to keep an explicit ``padding``, e.g. ``(0, 0)`` for abutting panels.
    :type tight_layout: bool
    :param shared_axes: Which axes the panels share, one of ``None``, ``'x'``, ``'y'``, or ``'both'``. ``'x'`` shares the x-axis per column (``sharex='col'``), ``'y'`` shares the y-axis per row (``sharey='row'``), and ``'both'`` shares both; a shared axis gives one common range per group and tick labels only on the outer edge. Defaults to None (independent axes).
    :type shared_axes: str | None
    :param height_ratios: The relative height of each row, e.g. ``[4, 1]`` for a 2-row grid whose second row is a quarter as tall as the first. Must have one entry per row. Defaults to None (equal heights).
    :type height_ratios: list[float] | None
    :param width_ratios: The relative width of each column, analogous to ``height_ratios``. Must have one entry per column. Defaults to None (equal widths).
    :type width_ratios: list[float] | None
    """

    type:str=field(repr=False, init=False, default='grid')

    plots:list[Plot]
    padding:tuple[float,float]=field(default=(0.2, 0.2))
    shape:tuple[int, int]
    size:tuple[float, float]|None=field(default=None)
    watermark:Watermark=field(default_factory=Watermark)
    watermark_plot:int|tuple[int, int]|None=field(default=None)
    tight_layout:bool=field(default=True)
    shared_axes:str|None=field(default=None)
    height_ratios:list[float]|None=field(default=None)
    width_ratios:list[float]|None=field(default=None)

    _api_doc = inspect.cleandoc("""
    Producing a Figure with a Grid of Plots
    ---------------------------------------

    This figure's type is ``grid``. It produces a grid of plots.

    The following keys are mandatory:

        * ``plots`` (*list* of :class:`Plot <eos.figure.plot.Plot>`) -- The list of plots to be drawn in the figure. Each plot can be any of the
            :class:`Plot <eos.figure.plot.Plot>` descendants. The list of plots are assigned to the grid positions in row-major order, i.e. the
            first plot is assigned to the first row and first column, the second plot to the first row and second column, etc.

        * ``shape`` (*tuple[int, int]*) -- The shape of the figure's grid, specifying the number of rows and columns in that order.

    The following keys are optional:

        * ``size`` (*tuple[float, float]*) -- The size of the figure in inches. Defaults to (3.0 * ncol, 3.0 * nrow).

        * ``watermark_plot`` (*int* or *tuple[int, int]*) -- The plot that carries the watermark, given either as a flattened
            (row-major) index or as a 2D ``(row, col)`` address. If omitted, every plot is stamped.

        * ``tight_layout`` (*bool*) -- Whether to lay out the grid with ``tight_layout``. Defaults to True. Set to False to keep an
            explicit ``padding`` (e.g. ``(0, 0)`` for abutting panels), which ``tight_layout`` would otherwise undo.

        * ``shared_axes`` (*str* or *None*) -- Which axes the panels share, one of ``None``, ``'x'``, ``'y'``, or ``'both'``. ``'x'`` shares
            the x-axis per column, ``'y'`` shares the y-axis per row, ``'both'`` shares both; a shared axis gives one common range per group
            and tick labels only on the outer edge. Defaults to None (independent axes).

        * ``height_ratios`` (*list[float]*) -- The relative height of each row, one entry per row, e.g. ``[4, 1]`` for a 2-row grid whose
            second row is a quarter as tall as the first. Defaults to None (equal heights).

        * ``width_ratios`` (*list[float]*) -- The relative width of each column, analogous to ``height_ratios``. Defaults to None (equal widths).


    """)
    def __post_init__(self):
        nrow, ncol = self.shape
        figsize = self.size if self.size is not None else (3.0 * ncol, 3.0 * nrow)
        self._figure = plt.figure(figsize=figsize)

        if self.height_ratios is not None and len(self.height_ratios) != nrow:
            raise ValueError(f"'height_ratios' must have {nrow} entries (one per row), got {len(self.height_ratios)}")
        if self.width_ratios is not None and len(self.width_ratios) != ncol:
            raise ValueError(f"'width_ratios' must have {ncol} entries (one per column), got {len(self.width_ratios)}")

        self._gridspec = self._figure.add_gridspec(nrow, ncol, hspace=self.padding[0], wspace=self.padding[1],
                                                    height_ratios=self.height_ratios, width_ratios=self.width_ratios)
        if self.shared_axes not in (None, 'x', 'y', 'both'):
            raise ValueError(f"'shared_axes' must be one of None, 'x', 'y', 'both', got {self.shared_axes!r}")
        sharex = 'col' if self.shared_axes in ('x', 'both') else False
        sharey = 'row' if self.shared_axes in ('y', 'both') else False
        # squeeze=False keeps axes a 2D array even for a 1x1 grid, so flatten() is always valid
        axes = self._gridspec.subplots(sharex=sharex, sharey=sharey, squeeze=False)
        self._axes = axes.flatten('C') # flatten to row-major style
        self._watermark_idx = self._resolve_watermark_plot(nrow, ncol)

    def validate_semantics(self, context):
        for index, plot in enumerate(self.plots):
            yield from (
                diagnostic.prefixed('plots', index)
                for diagnostic in plot.validate_semantics(context)
            )

    def _resolve_watermark_plot(self, nrow, ncol):
        "Resolve the watermark_plot field to a single flattened (row-major) index, or None for all plots."
        wp = self.watermark_plot
        if wp is None:
            return None

        nplots = len(self.plots)
        if isinstance(wp, bool):
            # bool is a subclass of int; reject it so that e.g. 'watermark_plot: true'
            # does not silently select plot index 1
            raise ValueError(f"'watermark_plot' must be an int or a (row, col) pair, got {wp!r}")
        elif isinstance(wp, int):
            idx = wp
        else: # 2D (row, col) address; a YAML list arrives here as well
            if len(wp) != 2:
                raise ValueError(f"'watermark_plot' must be an int or a (row, col) pair, got {wp}")
            row, col = wp
            if not (0 <= row < nrow and 0 <= col < ncol):
                raise ValueError(f"'watermark_plot' {tuple(wp)} is outside the {nrow}x{ncol} grid")
            idx = row * ncol + col

        if not (0 <= idx < nplots):
            raise ValueError(f"'watermark_plot' resolves to plot {idx}, but there are only {nplots} plots")

        return idx

    def draw(self, context:AnalysisFileContext=None):
        """Draw the grid figure.

        Use :meth:`save <eos.figure.Figure.save>` to store the drawn figure to one or more output files.

        :param context: The analysis file context, which contains the paths to the data files and other relevant information.
        :type context: :class:`AnalysisFileContext <eos.analysis_file_context.AnalysisFileContext>`
        """
        context = AnalysisFileContext() if context is None else context
        for idx, plot in enumerate(self.plots):
            plot.prepare(context)
            plot.draw(self._axes[idx])
            if self._watermark_idx is None or self._watermark_idx == idx:
                plot.draw_watermark(self._axes[idx], self.watermark)

        if self.tight_layout:
            self._gridspec.tight_layout(self._figure)

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`GridFigure` from its keyword description.

        Recursively deserializes each entry of the ``plots`` list via
        :class:`PlotFactory <eos.figure.plot.PlotFactory>`, as well as the optional ``watermark`` description.

        :param kwargs: The figure description. Must contain ``plots`` and ``shape`` keys.
        :returns: The instantiated figure.
        :rtype: GridFigure
        """
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs['plots'] = [PlotFactory.from_dict(**p) for p in _kwargs['plots']]
        if 'watermark' in _kwargs:
            _kwargs['watermark'] = Watermark.from_dict(**_kwargs['watermark'])
        return Deserializable.make(cls, **_kwargs)


@dataclass(kw_only=True)
class CornerFigure(Figure):
    r"""Produces a corner figure, i.e., a figure with a triangular arrangement of 1D and 2D marginal PDFs.

    Distributions of the variables are shown on the diagonal, while correlations are plotted on the lower-left corner.

    :param contents: The list of :class:`DataFile <eos.figure.DataFile>` objects to be used, each containing the path to the data file, its label, and optionally its color.
    :type contents: list[:class:`DataFile <eos.figure.DataFile>`]
    :param variables: The list of variable names to be shown. If not provided, all variables contained in the first data file are shown.
    :type variables: list[str] | None
    """

    type:str=field(repr=False, init=False, default='corner')

    contents:list[DataFile]
    variables:list[str]=None
    kde:bool=False

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

    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the corner figure for drawing.

        :param context: The analysis file context, which contains the paths to the data files and other relevant information.
        :type context: :class:`AnalysisFileContext <eos.analysis_file_context.AnalysisFileContext>`
        """
        context = AnalysisFileContext() if context is None else context
        for content in self.contents:
            content.prepare(context=context)

        if self.variables:
            self._variables = self.variables
        else:
            self._variables = self.contents[0].variables

        # determine useful ranges empirically
        absmin, absmax = np.array([+np.inf] * len(self._variables)), np.array([-np.inf] * len(self._variables))
        for content in self.contents:
            indices = [content._datafile.lookup_table[v] for v in self._variables]
            cmin, cmax = content.empirical_range
            for idx, cidx in enumerate(indices):
                absmin[idx] = cmin[cidx] if cmin[cidx] < absmin[idx] else absmin[idx]
                absmax[idx] = cmax[cidx] if cmax[cidx] > absmax[idx] else absmax[idx]

        # Check that the variables of the data files match
        for content in self.contents:
            unknown_variables = set(self._variables) - set(content.variables)
            if len(unknown_variables) > 0:
                raise ValueError(f"Unknown variables requested from data file '{content.path}': {list(unknown_variables)}")

        self._labels = self.contents[0].labels(self._variables)

        plots = []
        size = len(self._variables)

        for i in range(size):     # rows
            for j in range(size): # columns

                if i < j:
                    plots.append(PlotFactory.from_dict(**{
                        'type': 'empty'
                    }))

                elif i == j:
                    plots.append(PlotFactory.from_dict(**{
                        'xaxis': {
                            'ticks': { 'visible': True, 'position': 'both' },
                            'label': self._labels[j],
                            'range': [ absmin[j], absmax[j] ]
                        }
                        if (i == size - 1) else
                        {
                            'ticks': { 'visible': True, 'position': 'top' },
                            'range': [ absmin[j], absmax[j] ]
                        },
                        'yaxis': {
                            'ticks': { 'visible': False },
                            # 1D marginals, no label, no range
                        },
                        'grid': { 'visible': True, 'axis': 'x' },
                        'items': [
                            {
                                'type': 'kde1D', 'label': content.label,
                                'datafile': context.data_path(content.path),
                                'variable': self._variables[j],
                                'color': content.color,
                                'range': [absmin[j], absmax[j]]
                            } if content.kde else {
                                'type': 'histogram1D', 'label': content.label,
                                'datafile': context.data_path(content.path),
                                'variable': self._variables[j],
                                'color': content.color,
                                'range': [absmin[j], absmax[j]]
                            }
                        for content in self.contents]
                    }))

                else:
                    plots.append(PlotFactory.from_dict(**{
                        'xaxis': {
                            'ticks': { 'visible': True, 'position': 'bottom' },
                            'label': self._labels[j],
                            'range': [ absmin[j], absmax[j] ]
                        }
                        if (i == size - 1) else
                        {
                            'ticks': { 'visible': False, 'position': 'both' },
                            'range': [ absmin[j], absmax[j] ]
                        },
                        'yaxis': {
                            'ticks': { 'visible': True },
                            'label': self._labels[i],
                            'range': [ absmin[i], absmax[i] ]
                        }
                        if (j == 0) else
                        {
                            'ticks': { 'visible': False },
                            'range': [ absmin[i], absmax[i] ]
                        },
                        'grid': { 'visible': True},
                        'items': [
                            {
                                'type': 'kde2D', 'label': content.label,
                                'datafile': context.data_path(content.path),
                                'variables': [self._variables[j], self._variables[i]],
                                'color': content.color,
                                'contours': ['lines', 'areas'],
                                'xrange': [absmin[j], absmax[j]],
                                'yrange': [absmin[i], absmax[i]]
                            } if content.kde else
                            {
                                'type': 'histogram2D', 'label': content.label,
                                'datafile': context.data_path(content.path),
                                'variables': [self._variables[j], self._variables[i]],
                                'color': content.color,
                                'xrange': [absmin[j], absmax[j]],
                                'yrange': [absmin[i], absmax[i]]
                            }
                        for content in self.contents]
                    }))

        self._figure = GridFigure(shape=(size, size), plots=plots, padding=(0.0, 0.0))


    def draw(self, context:AnalysisFileContext=None):
        """Draw the corner figure.

        Use :meth:`save <eos.figure.CornerFigure.save>` to store the drawn figure to one or more output files.

        :param context: The analysis file context, which contains the paths to the data files and other relevant information.
        :type context: :class:`AnalysisFileContext <eos.analysis_file_context.AnalysisFileContext>`
        """
        context = AnalysisFileContext() if context is None else context
        if not hasattr(self, '_figure'):
            self.prepare(context=context)

        self._figure.draw(context=context)

    def save(self, output:str|list[str]):
        """Save the corner figure to one or more output files.

        This figure's contents are drawn by an underlying :class:`GridFigure <eos.figure.GridFigure>`,
        to which saving is delegated.

        :param output: The path(s) to the output file(s) where the figure is saved. The file format of each
            output file is determined by its file name extension.
        :type output: str | list[str]
        """

        self._figure.save(output)


    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`CornerFigure` from its keyword description.

        Recursively deserializes each entry of the ``contents`` list into a
        :class:`DataFile <eos.figure.data.DataFile>` instance.

        :param kwargs: The figure description. Must contain a ``contents`` key.
        :returns: The instantiated figure.
        :rtype: CornerFigure
        """
        _kwargs = _copy.deepcopy(kwargs)
        if 'contents' in _kwargs:
            _kwargs['contents'] = [DataFile.from_dict(**c) for c in _kwargs['contents']]
        return Deserializable.make(cls, **_kwargs)

@dataclass(kw_only=True)
class OverviewSource(Deserializable):
    type: str
    label: str
    color: str = None
    markerstyle: str = None
    linewidth: float = None
    linestyle: str = None
    alpha: float = None

    datafiles: list[str] = None
    names: list[str] = None

@dataclass(kw_only=True)
class OverviewFigure(Figure):
    """Produces a figure with a single plot, giving an overview of 1D predictions.
       Prints as eos.info the summary of overview data.

    :param legend: As in :class:`Plot <eos.figure.plot.Plot>`, position of the legend, always drawn
        outside of the main plot.
    :type legend: dict
    :param xaxis: Similar to :class:`Plot <eos.figure.plot.Plot>`, xaxis definition. May contain an
        optional ``normalize`` key (``None``, ``'constraint'``, or ``'posterior'``) and, if
        ``normalize == 'posterior'``, a ``normalize-source`` key naming the source label to normalize to.
    :type xaxis: dict
    :param yaxis: Optional yaxis spacing. Recognized keys are ``observable-offset`` (space between rows
        for different observables, default 1.0) and ``source-offset`` (space between sources within one
        row, default 0.05).
    :type yaxis: dict
    :param observables: List of named observables to plot; nesting into a list of lists causes a
        horizontal separator line to be drawn between the groups.
    :type observables: list[str] | list[list[str]]
    :param sources: Sources to compare. Each is either a prediction (``type: 'prediction'``, with a list
        of sample-file paths under ``datafiles``) or a set of named EOS database constraints
        (``type: 'constraint'``, with a list of constraint names under ``names``), plus optional
        ``label``, ``color``, ``markerstyle``, ``linewidth``, ``linestyle``, ``alpha``.
    :type sources: list[dict]
    :param size: The size of the figure in inches. Defaults to an automatic calculation based on the
        number of observables and sources.
    :type size: tuple[float, float]
    :param watermark: The optional specification of the EOS watermark. See
        :class:`Watermark <eos.figure.Watermark>`.
    :type watermark: :class:`Watermark <eos.figure.Watermark>`
    """

    type: str = field(repr=False, init=False, default='overview')

    legend: dict = field(default_factory=lambda: {'position': 'upper center'})
    xaxis: dict = None
    yaxis: dict = field(default_factory=dict)
    observables: list = None
    sources: list = None

    size: tuple[float, float] = None
    watermark: Watermark = field(default_factory=Watermark)

    _NORMALIZE_MODES = (None, 'constraint', 'posterior')
    _DEFAULT_MARKERS = ('o', 's', '^', 'D', 'v', 'P', 'X')

    _api_doc = inspect.cleandoc("""
    Producing a Figure with a single Plot
    -------------------------------------

    This figure's type is ``overview``. It produces a single plot containing an overview comparison of
    predictions and/or measurements for a set of observables.

    The following keys are mandatory:
       * ``xaxis`` (*dict*) -- The x-axis label and range, and optionally ``normalize``.
       * ``observables`` (*list[str]* or *list[list[str]]*) -- The observables entering the overview plot.
       * ``sources`` (*list* of :class:`DataFile <eos.figure.data.DataFile>`) -- The predictions/constraints
         to be drawn.

    The following keys are optional:
        * ``legend`` (*dict*) -- Position and number of columns of the legend. Defaults to ``{'position': 'upper center', 'ncol': 1}``.
        * ``yaxis`` (*dict*) -- ``observable-offset`` and ``source-offset`` spacing. Defaults to 1.0 and 0.05.
        * ``size`` (*tuple[float, float]*) -- Defaults to an automatic size.
    """)

    def __post_init__(self):
        if not self.xaxis:
            raise ValueError("xaxis must be defined.")
        if not self.observables:
            raise ValueError("observables must include at least one item to be plotted.")
        if not self.sources:
            raise ValueError("sources must include at least one item to be plotted.")

        self._normalize = self.xaxis.get('normalize', None)
        if self._normalize not in self._NORMALIZE_MODES:
            raise ValueError(f"xaxis['normalize'] must be one of {self._NORMALIZE_MODES}, got {self._normalize!r}.")

        self._normalize_source = self.xaxis.get('normalize-source', None)
        if self._normalize == 'posterior' and not self._normalize_source:
            raise ValueError("xaxis['normalize'] == 'posterior' requires xaxis['normalize-source'] "
                              "(the label of the source to normalize to).")

    @staticmethod
    def _group(entries:list):
        #if not any(isinstance(entry, list) for entry in entries):
        #    groups = [[entry] for entry in entries]
        #else:
        #    groups = [entry if isinstance(entry, list) else [entry] for entry in entries]
        groups = [group if isinstance(group, list) else [group] for group in entries]
        flat, group_ids, labels = [], [], []
        for group_id, group in enumerate(groups):
            for entry in group:
                if not isinstance(entry, dict):
                    raise TypeError("Each observable entry must be a dictionary of the form {observable name: latex label}.")
                if len(entry) != 1:
                    raise ValueError("Each observable entry must contain exactly one key-value pair.")
                observable, label = next(iter(entry.items()))
                flat.append(observable)
                group_ids.append(group_id)
                labels.append(label)

        return flat, group_ids, labels

    #@staticmethod
    #def _observable_label(observable):
    #    name, _, options = observable.partition(';')
    #    base_latex = Observables()[name].latex()
    #    if not (base_latex.startswith('$') and base_latex.endswith('$')):
    #        base_latex = f'${base_latex}$'
    #    return f'{base_latex} | {options}' if options else base_latex

    @staticmethod
    def _prediction_quantiles(samples, weights, level=68.27e-2):
        """Weighted median and asymmetric 1-sigma uncertainty widths as for `uncertainty` items."""
        half = level / 2.0
        interval = [0.5 - half, 0.5, 0.5 + half]
        lower, central, higher = np.quantile(samples, q = interval, weights = weights, method='inverted_cdf', axis=0)
        err_low = central - lower
        err_high = higher - central
        return central, err_low, err_high

    @staticmethod
    def _observable_matches(query, candidate):
        query_name, _, query_options = query.partition(';')

        candidate_no_kin = candidate.split('[', 1)[0]
        candidate_name, _, candidate_options = candidate_no_kin.partition(';')

        if query_name != candidate_name:
            return False

        query_options = set(query_options.split(',')) if query_options else set()
        candidate_options = set(candidate_options.split(',')) if candidate_options else set()

        return query_options.issubset(candidate_options)

    @staticmethod
    def _prediction_observable_index(prediction, observable):
        for key, index in prediction.lookup_table.items():
            if OverviewFigure._observable_matches(observable, key):
                return index
        return None

    def _prediction_entries_for_source(self, source, flat_obs):
        results = {}
        datafiles = source.datafiles or []
        if isinstance(datafiles, str):
            datafiles = [datafiles]

        for obs in flat_obs:
            found = False
            for path in datafiles:
                prediction = eos.data.Prediction(path)
                idx = self._prediction_observable_index(prediction, obs)
                if idx is None:
                    continue

                if obs in results:
                    raise ValueError(
                        f"source '{source.label}' provides observable '{obs}' "
                        f"from more than one data file."
                    )

                central, err_low, err_high = self._prediction_quantiles(
                    prediction.samples[:, idx],
                    prediction.weights
                )

                results[obs] = (central, err_low, err_high)
                found = True
                break
        if not found:
            raise KeyError(f"No matching observable found for {obs} in source '{source.label}'")
        return results

    @staticmethod
    def _constraint_entry(name):
        constraints = eos.Constraints()
        if name not in constraints:
            raise KeyError(f"Constraint '{name}' not found in EOS database.")

        entry = constraints[name]
        kind = entry.type() if callable(getattr(entry, 'type', None)) else getattr(entry, 'type', 'Gaussian')
        if kind != 'Gaussian':
            raise NotImplementedError(
                f"Constraint '{name}' has type '{kind}'; only Gaussian constraints are "
                "currently supported in OverviewFigure."
            )

        return entry

    def _constraint_entries_for_source(self, source, flat_obs):
        results = {}
        matched_constraints = {}

        for name in source.names:
            entry = self._constraint_entry(name)
            const = _yaml.safe_load(entry.serialize())
            raw_observables = const.get('observable', [])
            nameobs, _, optionsobs = raw_observables.partition(';')
            if optionsobs == '':
                optionsobs = const.get('options', {})
                option_suffix =  ';' + ','.join(f'{k}={v}' for k, v in optionsobs.items()) if optionsobs else ''
            else:
                option_suffix = ';' + optionsobs
            mean = const['mean']
            sigma_stat = const['sigma-stat']
            sigma_sys = (const['sigma-sys'] if 'sigma-sys' in const else {'hi': 0.0, 'lo': 0.0})
            err_high = np.sqrt(sigma_stat['hi'] ** 2 + sigma_sys['hi'] ** 2)
            err_low = np.sqrt(sigma_stat['lo'] ** 2 + sigma_sys['lo'] ** 2)

            obsfull = nameobs + option_suffix

            for obs in flat_obs:
                if not self._observable_matches(obs, obsfull):
                    continue
                if obs in results:
                    raise ValueError(
                        f"observable '{obs}' is constrained by both '{matched_constraints[obs]}' and '{name}'; "
                        "only one named constraint per observable is allowed."
                    )
                results[obs] = (mean, err_low, err_high)
                matched_constraints[obs] = name

        return results

    def info_summary(self, context: AnalysisFileContext = None):
        """Print each observable's raw (un-normalized) central value and uncertainty for every
        source, as a sanity check before drawing (e.g. to catch an inverted or mismatched entry)."""
        context = AnalysisFileContext() if context is None else context
        flat_obs, _, _ = self._group(self.observables)

        per_source_entries = {}
        for source in self.sources:
            if source.type == 'prediction':
                per_source_entries[source.label] = self._prediction_entries_for_source(source, flat_obs)
            elif source.type == 'constraint':
                per_source_entries[source.label] = self._constraint_entries_for_source(source, flat_obs)
            else:
                raise ValueError(f"unknown source type '{source.type}' for source '{source.label}'.")

        headers = ['observable'] + [source.label for source in self.sources]
        rows = []
        for obs in flat_obs:
            row = [obs]
            for source in self.sources:
                entry = per_source_entries[source.label].get(obs)
                if entry is None:
                    row.append('--')
                else:
                    central, err_low, err_high = entry
                    row.append(f'{central:.4g} (+{err_high:.2g}/-{err_low:.2g})')
            rows.append(row)

        widths = [max(len(str(row[i])) for row in ([headers] + rows)) for i in range(len(headers))]
        def fmt_row(row):
            return '  '.join(str(cell).ljust(w) for cell, w in zip(row, widths))

        eos.info(fmt_row(headers))
        eos.info('  '.join('-' * w for w in widths))
        for row in rows:
            eos.info(fmt_row(row))

    def prepare(self, context: AnalysisFileContext = None):
        """Prepare the overview figure for drawing.

        :param context: The analysis file context, which contains the paths to the data files and other
            relevant information.
        :type context: :class:`AnalysisFileContext <eos.analysis_file_context.AnalysisFileContext>`
        """
        context = AnalysisFileContext() if context is None else context
        flat_obs, obs_group_ids, obs_labels = self._group(self.observables)
        n_obs = len(flat_obs)
        n_sources = len(self.sources)

        self._obs_offset = self.yaxis.get('observable-offset', 1.0)
        if self._obs_offset < 1.0:
            raise ValueError("observable-offset must be at least 1.0")
        self._source_offset = self.yaxis.get('source-offset', 0.05)


        self._yvals = [1.0 + i * self._obs_offset for i in range(n_obs)]
        self._yvals.reverse()
        self._ylabels = obs_labels #[self._observable_label(obs) for obs in flat_obs]

        self._sourced_entries = self._build_sourced_entries(flat_obs, self._yvals, n_sources)

        # horizontal separators between observable groups, at the midpoint
        # between the last row of one group and the first row of the next
        separator_ys = [
            (self._yvals[i - 1] + self._yvals[i]) / 2.0
            for i in range(1, n_obs)
            if obs_group_ids[i] != obs_group_ids[i - 1]
        ]
        items = []
        if self._normalize:
            items.append({'type': 'vertical', 'x': 1.0, 'color': 'lightgray'})

        items += [
            {
                'type': 'errorbars',
                'positions': entry['positions'],
                'xerrors': entry['xerrors'],
                'yerrors': [0.0] * len(entry['positions']),
                'color': entry['color'],
                'alpha': entry['alpha'],
                'marker': entry['marker'],
                'linestyle': entry['linestyle'],
                'linewidth': entry['linewidth'],
                'label': entry['label'],
            }
            for entry in self._sourced_entries
        ]
        items += [{'type': 'expression', 'expression': f'{y}', 'color': 'lightgray', 'range': self.xaxis['range']} for y in separator_ys]
        plot = PlotFactory.from_dict(**{
            'legend': {'position': self.legend.get('position', 'upper center')},
            'xaxis': {
                'ticks': {'visible': True, 'position': 'bottom'},
                'label': self.xaxis['label'],
                'range': self.xaxis['range'],
            },
            'yaxis': {
                'ticks': {'visible': True},
                'range': [1.5 - self._obs_offset, n_obs * self._obs_offset + self._obs_offset],
            },
            'items': items,
        })

        size = self.size if self.size is not None else (6.4, 0.6 * (n_obs * self._obs_offset) + 1.2)
        self._figure = SingleFigure(size=size, plot=plot)
        self._ax = self._figure._ax
        self._ax.yaxis.set_ticks(self._yvals)
        self._ax.yaxis.set_ticklabels(self._ylabels)

    def _build_sourced_entries(self, flat_obs, yvals, n_sources):
        """Build one drawable entry per source, aggregating across every observable it covers."""
        default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

        # each source's per-observable (central, err_low, err_high), computed once
        per_source_entries = {}
        for source in self.sources:
            if source.type == 'prediction':
                per_source_entries[id(source)] = self._prediction_entries_for_source(source, flat_obs)
            elif source.type == 'constraint':
                per_source_entries[id(source)] = self._constraint_entries_for_source(source, flat_obs)
            else:
                raise ValueError(f"unknown source type '{source.type}' for source '{source.label}'.")

        # resolve the normalization reference, if any
        reference = {}
        if self._normalize == 'constraint':
            claimed_by = {}
            for source in self.sources:
                if source.type != 'constraint':
                    continue
                for obs, value in per_source_entries[id(source)].items():
                    if obs in reference:
                        raise ValueError(
                            f"observable '{obs}' is constrained by more than one constraint source "
                            f"('{claimed_by[obs]}' and '{source.label}'); cannot normalize unambiguously."
                        )
                    reference[obs] = value
                    claimed_by[obs] = source.label
        elif self._normalize == 'posterior':
            ref_source = next((s for s in self.sources if s.label == self._normalize_source), None)
            if ref_source is None:
                raise ValueError(f"normalize-source '{self._normalize_source}' matches no source label.")
            reference = per_source_entries[id(ref_source)]

        # symmetric stagger of sources around each observable's row
        stagger = [(-(n_sources - 1) / 2.0 + i) * self._source_offset for i in range(n_sources)]

        entries = []
        for source_index, source in enumerate(self.sources):
            source_data = per_source_entries[id(source)]
            positions, xerrors = [], []

            for obs_index, obs in enumerate(flat_obs):
                if obs not in source_data:
                    continue

                central, err_low, err_high = source_data[obs]

                if self._normalize and obs in reference:
                    ref_central = reference[obs][0]
                    err_low, err_high = err_low / ref_central, err_high / ref_central
                    central = central / ref_central

                y = yvals[obs_index] + stagger[source_index]
                positions.append([central, y])
                xerrors.append([err_low, err_high])

            if not positions:
                continue

            entries.append({
                'positions': positions,
                'xerrors': xerrors,
                'color': getattr(source, 'color', None) or default_colors[source_index % len(default_colors)],
                'marker': getattr(source, 'markerstyle', None) or self._DEFAULT_MARKERS[source_index % len(self._DEFAULT_MARKERS)],
                'linestyle': getattr(source, 'linestyle', None) or 'none',
                'alpha': getattr(source, 'alpha', None) or 1.0,
                'linewidth': getattr(source, 'linewidth', None) or 1.5,
                'label': source.label,
            })

        return entries

    _LEGEND_OUTSIDE = {
        'upper center': ('lower center', (0.5, 1.02)),
        'lower center': ('upper center', (0.5, -0.02)),
        'upper left':   ('lower left',   (0.0, 1.02)),
        'upper right':  ('lower right',  (1.0, 1.02)),
        'lower left':   ('upper left',   (0.0, -0.02)),
        'lower right':  ('upper right',  (1.0, -0.02)),
        'center left':  ('center right', (-0.02, 0.5)),
        'center right': ('center left',  (1.02, 0.5)),
    }

    def _place_legend_outside(self, ax):
        handles = [
            ax.errorbar(
                [], [], xerr=[],
                fmt=entry['marker'], color=entry['color'],
                linestyle=entry['linestyle'], linewidth=entry['linewidth'],
                label=entry['label'],
            )
            for entry in self._sourced_entries
        ]
        position = self.legend.get('position', 'upper center')
        loc, anchor = self._LEGEND_OUTSIDE.get(position, ('lower center', (0.5, 1.02)))
        ax.legend(handles=handles, loc=loc, bbox_to_anchor=anchor,
                  ncol=self.legend.get('ncol', 1),
                  frameon=self.legend.get('frameon', True))

    def draw(self, context: AnalysisFileContext = None):
        """Draw the single-plot overview figure.

        :param context: The analysis file context, which contains the paths to the data files and other
            relevant information.
        :type context: :class:`AnalysisFileContext <eos.analysis_file_context.AnalysisFileContext>`
        """
        context = AnalysisFileContext() if context is None else context
        if not hasattr(self, '_figure'):
            self.prepare(context)
        self._figure.draw(context)
        self._ax.yaxis.set_ticks(self._yvals)
        self._ax.yaxis.set_ticklabels(self._ylabels)
        self._place_legend_outside(self._ax)
        self.watermark.draw(self._ax)
        self.info_summary(context)

    def save(self, output: str | list[str]):
        """Save the overview figure to one or more output files.

        :param output: The path(s) to the output file(s) where the figure is saved.
        :type output: str | list[str]
        """
        self._figure.save(output)

    @classmethod
    def from_dict(cls, **kwargs):
        _kwargs = _copy.deepcopy(kwargs)
        _kwargs['sources'] = [
            OverviewSource.from_dict(**source)
            for source in _kwargs['sources']
        ]
        if 'watermark' in _kwargs:
            _kwargs['watermark'] = Watermark.from_dict(**_kwargs['watermark'])
        return Deserializable.make(cls, **_kwargs)

class FigureFactory:
    r"""Factory class to create figures from a dictionary description.

    This class provides a convenient way to create figures from a dictionary description, which can be used to deserialize figures.
    The factory method :py:meth:`from_dict` should be used to create a figure from a dictionary description.
    The description contains the arguments used for initializing objects of classes within the :py:mod:`eos.figure` module.

    A figure can contain one or more individual plots. Each plot contains one or more plottable items, such as a data series from the evaluation of an observable or a kernel density estimate to visualize a set of statistical samples.

    The factory method :py:meth:`from_yaml` can be used to create a figure from a YAML description.
    """

    # Also build the documentation based on ordered registry
    # Initializer is well-defined for python version >= 3.6
    registry = {
        'single':   SingleFigure, # default
        'inset':    InsetFigure,
        'grid':     GridFigure,
        'corner':   CornerFigure,
        'overview': OverviewFigure
    }

    @staticmethod
    def from_yaml(yaml_data:str):
        """Factory method to create a figure from a YAML description.

        This method takes a YAML description of a figure, turns it into a `dict`,
        and forwards it to the :py:meth:`from_dict` method.

        :param yaml_data: A YAML description of the figure.
        :type yaml_data: str
        :returns: A descendant of the :py:class:`eos.figure.Figure` class
        :rtype: :py:class:`eos.figure.Figure`
        """
        kwargs = _yaml.safe_load(yaml_data)
        return FigureFactory.from_dict(**kwargs)

    @staticmethod
    def from_dict(**kwargs):
        r"""Factory method to create a plot from a dictionary

        This method takes a dictionary description of a figure and uses matplotlib to render that
        figure conveniently and at a publication-grade level.

        The description contains the arguments used for initializing objects of classes within the
        :py:mod:`eos.figure` module.

        A figure can contain one or more individual plots. Each plot contains one or more plottable
        items, such as a data series from the evaluation of an observable or a kernel density estimate to
        visualize a set of statistical samples.

        :param description: A figure description as required by any of the descendant classes of :py:class:`eos.figure.Figure`.
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
