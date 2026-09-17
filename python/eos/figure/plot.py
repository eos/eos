# Copyright (c) 2026      Carolina Bolognani
# Copyright (c) 2023-2026 Danny van Dyk
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
from eos.analysis_file_context import AnalysisFileContext
from eos.deserializable import Deserializable
from eos.diagnostic import Diagnostic, Severity, _check_qualified

from .item import Item, ItemFactory, ItemColorCycler

import copy as _copy
import eos
import eos.data
import inspect
import matplotlib
import matplotlib.pyplot as plt
import numpy as _np
import yaml as _yaml


@dataclass
class Legend(Deserializable):
    r"""Represents the legend properties in a plot.

    :param position: The position of the legend in the plot. Defaults to 'best', which automatically places the legend in the best location.
    :type position: str
    :param ncol: The number of columns of the legend. Defaults to 1.
    :type ncol: int
    :param outside: Whether to draw the legend outside of the axes, adjacent to the edge named by
        ``position``. Defaults to False, which draws the legend inside the axes. Requires a
        ``position`` other than 'best', which has no edge to attach to.
    :type outside: bool
    """

    position:str=field(default='best')
    ncol:int=field(default=1)
    outside:bool=field(default=False)

    # for each position, the anchor point of the legend's own box and the point of the axes it is
    # placed against, so that the legend sits just outside the named edge
    _OUTSIDE = {
        'upper center': ('lower center', (0.5, 1.02)),
        'lower center': ('upper center', (0.5, -0.02)),
        'upper left':   ('lower left',   (0.0, 1.02)),
        'upper right':  ('lower right',  (1.0, 1.02)),
        'lower left':   ('upper left',   (0.0, -0.02)),
        'lower right':  ('upper right',  (1.0, -0.02)),
        'center left':  ('center right', (-0.02, 0.5)),
        'center right': ('center left',  (1.02, 0.5)),
    }

    def __post_init__(self):
        if self.ncol < 1:
            raise ValueError(f"'ncol' must be at least 1, got {self.ncol}")

        if self.outside and self.position not in self._OUTSIDE:
            raise ValueError(
                f"Legend position '{self.position}' cannot be drawn outside of the axes; "
                f"must be one of {', '.join(sorted(self._OUTSIDE))}"
            )

    def draw(self, ax, entries=None):
        """Draw the legend on the axes.

        :param ax: The axes to draw the legend on.
        :type ax: matplotlib.axes.Axes
        :param entries: A list of tuples (handle, label) for the legend entries. If None, uses the entries from the items.
        :type entries: list[tuple[matplotlib.artist.Artist, str]] | None
        """
        if self.outside:
            loc, anchor = self._OUTSIDE[self.position]
            kwargs = { 'loc': loc, 'bbox_to_anchor': anchor, 'ncol': self.ncol }
        else:
            kwargs = { 'loc': self.position, 'ncol': self.ncol }

        if entries is not None:
            # Nothing labelled to show: drawing an empty legend is pointless and, as of
            # matplotlib 3.10, passing empty handles/labels raises instead of warning.
            if not entries:
                return
            handles = [entry[0] for entry in entries]
            labels = [entry[1] for entry in entries]
            ax.legend(handles=handles, labels=labels, **kwargs)
        else:
            ax.legend(**kwargs)


@dataclass
class Grid(Deserializable):
    r"""Represents the grid properties in a plot.

    :param visible: Whether the grid is visible. Defaults to False.
    :type visible: bool
    :param axis: The axis on which the grid is drawn. Can be 'both', 'x', or 'y'. Defaults to 'both'.
    :type axis: str
    """

    visible:bool=field(default=False)
    axis:str=field(default='both')

    def __post_init__(self):
        if self.axis not in ['both', 'x', 'y']:
            raise ValueError(f'Unknown axis: {self.axis}')

    def draw(self, ax):
        """Draw the grid on the provided axes if it is visible.

        :param ax: The matplotlib axes onto which the grid is drawn.
        :type ax: matplotlib.axes.Axes
        """
        if self.visible:
            ax.grid(visible=self.visible, axis=self.axis, alpha=0.3)


@dataclass
class XTicks(Deserializable):
    r"""Represents the x axis ticks properties in a plot.

    :param minor: Whether to show minor ticks. Defaults to True.
    :type minor: bool
    :param position: The position of the ticks. Can be 'bottom', 'top', or 'both'. Defaults to 'bottom'.
    :type position: str
    :param visible: Whether the ticks are visible. Defaults to True.
    :type visible: bool
    :param format: A printf-style format string for the major tick labels, e.g. '%.2f'. Defaults to None (matplotlib's default formatter).
    :type format: str | None
    :param scaling_factor: A non-zero factor by which each major tick value is divided before being
        displayed, e.g. ``1e-3`` to label an axis in units of :math:`10^{-3}`. The (rescaled) value is
        rendered with ``format`` if given, or with ``'%g'`` otherwise. The magnitude removed in this
        way should be communicated through the axis ``label`` or ``unit``. Intended for linear axes.
        Defaults to None (no rescaling).
    :type scaling_factor: float | None
    :param locations: Explicit positions for the major ticks, e.g. ``[-2, 0, 2]``. Overrides the
        automatic tick placement, which otherwise chooses a step (and so which values end up
        labelled) based on the axes' rendered size rather than on the data range alone. Defaults to
        None (automatic placement).
    :type locations: list[float] | None
    """

    minor:bool=field(default=True)
    position:str=field(default='bottom')
    visible:bool=field(default=True)
    format:str|None=field(default=None)
    scaling_factor:float|None=field(default=None)
    locations:list[float]|None=field(default=None)

    def __post_init__(self):
        POSITIONS = ['bottom', 'top', 'both']
        if self.position not in POSITIONS:
            raise ValueError(f'Unknown position: {self.position}. Must be one of {",".join(POSITIONS)}')

        if self.scaling_factor is not None and self.scaling_factor == 0.0:
            raise ValueError("'scaling_factor' must be non-zero")

    def draw(self, ax):
        """Apply the x axis tick settings to the provided axes.

        Configures the visibility, position (bottom/top/both), and major/minor locators of the
        x axis ticks, accounting for linear and logarithmic scales.

        :param ax: The matplotlib axes whose x axis ticks are configured.
        :type ax: matplotlib.axes.Axes
        """
        if not self.visible:
            ax.xaxis.set_major_formatter(plt.NullFormatter())
            ax.xaxis.set_minor_formatter(plt.NullFormatter())
            ax.xaxis.set_tick_params(bottom=False, top=False)
        else:
            log_scale = ax.xaxis.get_scale() == 'log'
            if self.locations is not None:
                ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(self.locations))
            elif not log_scale:
                ax.xaxis.set_major_locator(matplotlib.ticker.AutoLocator())
            else:
                ax.xaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10.0))
            if self.minor:
                if not log_scale:
                    ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
                else:
                    ax.xaxis.set_minor_locator(matplotlib.ticker.LogLocator(base=10.0, subs='auto'))
            if self.scaling_factor is not None:
                if log_scale:
                    eos.warn("Tick 'scaling_factor' is intended for linear axes; applying it to a logarithmic x axis may yield misleading tick labels")
                fmt = self.format if self.format is not None else '%g'
                ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
                    lambda x, pos, s=self.scaling_factor, f=fmt: f % (x / s)))
            elif self.format is not None:
                ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter(self.format))
            ax.xaxis.set_tick_params(
                which = 'both' if self.minor else 'major',
                bottom=(self.position == 'bottom' or self.position == 'both'),
                labelbottom=(self.position == 'bottom' or self.position == 'both'),
                top=(self.position == 'top' or self.position == 'both'),
                labeltop=(self.position == 'top' or self.position == 'both')
            )


@dataclass
class XAxis(Deserializable):
    r"""Represents the x axis properties in a plot.

    :param label: The label for the x axis.
    :type label: str
    :param range: The range of the x axis as a tuple (min, max). Defaults to None, which means the range is determined empirically.
    :type range: tuple[float, float] | None
    :param ticks: The x axis ticks properties. Defaults to an instance of :class:`XTicks <eos.figure.XTicks>`.
    :type ticks: :class:`XTicks <eos.figure.XTicks>`
    :param unit: The unit of the x axis, e.g., 'GeV'. Defaults to None.
    :type unit: str | None
    :param scale: The scale of the x axis, e.g., 'linear' or 'log'. Defaults to 'linear'.
    :type scale: str
    """

    label:str=field(default=None)
    range:tuple[float, float]=field(default=None)
    ticks:XTicks=field(default_factory=XTicks)
    unit:str=None
    scale:str='linear'

    def __post_init__(self):
        if self.range is not None and len(self.range) != 2:
            raise ValueError("Range must be a tuple of two values (min, max)")

        if self.range is not None:
            self.range = tuple(float(x) for x in self.range)

    def draw(self, ax):
        """Apply the x axis settings to the provided axes.

        Sets the axis label (including the unit, if given), the axis range and scale, and draws the ticks.

        :param ax: The matplotlib axes whose x axis is configured.
        :type ax: matplotlib.axes.Axes
        """
        if self.label is not None and self.unit is not None:
            ax.set_xlabel(f'{self.label} [{self.unit}]')
        elif self.label is not None:
            ax.set_xlabel(self.label)

        if self.range is not None:
            ax.set_xlim(self.range)
        ax.set_xscale(self.scale)
        self.ticks.draw(ax)

    @classmethod
    def from_dict(cls, **kwargs):
        """Create an :class:`XAxis` from its keyword description.

        Recursively deserializes the nested ``ticks`` description into an :class:`XTicks` instance.

        :param kwargs: The x axis description.
        :returns: The instantiated x axis.
        :rtype: XAxis
        """
        _kwargs = _copy.deepcopy(kwargs)
        if 'ticks' in _kwargs:
            _kwargs['ticks'] = XTicks.from_dict(**_kwargs['ticks'])
        return Deserializable.make(cls, **_kwargs)


@dataclass
class YTicks(Deserializable):
    r"""Represents the y axis ticks properties in a plot.

    :param minor: Whether to show minor ticks. Defaults to True.
    :type minor: bool
    :param position: The position of the ticks. Can be 'left', 'right', or 'both'. Defaults to 'left'.
    :type position: str
    :param visible: Whether the ticks are visible. Defaults to True.
    :type visible: bool
    :param format: A printf-style format string for the major tick labels, e.g. '%.2f'. Defaults to None (matplotlib's default formatter).
    :type format: str | None
    :param scaling_factor: A non-zero factor by which each major tick value is divided before being
        displayed, e.g. ``1e-3`` to label an axis in units of :math:`10^{-3}`. The (rescaled) value is
        rendered with ``format`` if given, or with ``'%g'`` otherwise. The magnitude removed in this
        way should be communicated through the axis ``label`` or ``unit``. Intended for linear axes.
        Defaults to None (no rescaling).
    :type scaling_factor: float | None
    :param locations: Explicit positions for the major ticks, e.g. ``[-2, 0, 2]``. Overrides the
        automatic tick placement, which otherwise chooses a step (and so which values end up
        labelled) based on the axes' rendered size rather than on the data range alone. Defaults to
        None (automatic placement).
    :type locations: list[float] | None
    :param labels: Explicit labels for the major ticks, one per entry of ``locations``, e.g. the
        LaTeX names of the quantities shown along a categorical axis. Requires ``locations`` and
        replaces any formatting that ``format`` or ``scaling_factor`` would apply. Defaults to None
        (the tick values themselves are labelled).
    :type labels: list[str] | None
    """

    minor:bool=field(default=True)
    position:str=field(default='left')
    visible:bool=field(default=True)
    format:str|None=field(default=None)
    scaling_factor:float|None=field(default=None)
    locations:list[float]|None=field(default=None)
    labels:list[str]|None=field(default=None)

    def __post_init__(self):
        POSITIONS = ['left', 'right', 'both']
        if self.position not in POSITIONS:
            raise ValueError(f'Unknown position: {self.position}. Must be one of {",".join(POSITIONS)}')

        if self.scaling_factor is not None and self.scaling_factor == 0.0:
            raise ValueError("'scaling_factor' must be non-zero")

        if self.labels is not None:
            if self.locations is None:
                raise ValueError("'labels' requires 'locations'")
            if len(self.labels) != len(self.locations):
                raise ValueError(
                    f"'labels' has {len(self.labels)} entries, but 'locations' has {len(self.locations)}"
                )

    def draw(self, ax):
        """Apply the y axis tick settings to the provided axes.

        Configures the visibility, position (left/right/both), and major/minor locators of the
        y axis ticks, accounting for linear and logarithmic scales.

        :param ax: The matplotlib axes whose y axis ticks are configured.
        :type ax: matplotlib.axes.Axes
        """
        if not self.visible:
            ax.yaxis.set_major_formatter(plt.NullFormatter())
            ax.yaxis.set_minor_formatter(plt.NullFormatter())
            ax.yaxis.set_tick_params(left=False, right=False)
        else:
            log_scale = ax.yaxis.get_scale() == 'log'
            if self.locations is not None:
                ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(self.locations))
            elif not log_scale:
                ax.yaxis.set_major_locator(matplotlib.ticker.AutoLocator())
            else:
                ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10.0))
            if self.minor:
                if not log_scale:
                    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
                else:
                    ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(base=10.0, subs='auto'))
            if self.labels is not None:
                ax.yaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(self.labels))
            elif self.scaling_factor is not None:
                if log_scale:
                    eos.warn("Tick 'scaling_factor' is intended for linear axes; applying it to a logarithmic y axis may yield misleading tick labels")
                fmt = self.format if self.format is not None else '%g'
                ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
                    lambda y, pos, s=self.scaling_factor, f=fmt: f % (y / s)))
            elif self.format is not None:
                ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter(self.format))
            ax.yaxis.set_tick_params(
                left=(self.position == 'left' or self.position == 'both'),
                right=(self.position == 'right' or self.position == 'both')
            )

@dataclass
class YAxis(Deserializable):
    r"""Represents the y axis properties in a plot.

    :param label: The label for the y axis.
    :type label: str
    :param range: The range of the y axis as a tuple (min, max). Defaults to None, which means the range is determined empirically.
    :type range: tuple[float, float] | None
    :param ticks: The y axis ticks properties. Defaults to an instance of :class:`YTicks <eos.figure.YTicks>`.
    :type ticks: :class:`YTicks <eos.figure.YTicks>`
    :param unit: The unit of the y axis, e.g., 'GeV'. Defaults to None.
    :type unit: str | None
    :param scale: The scale of the y axis, e.g., 'linear' or 'log'. Defaults to 'linear'.
    :type scale: str
    """

    label:str=None
    range:tuple[float, float]=None
    ticks:YTicks=field(default_factory=YTicks)
    unit:str=None
    scale:str='linear'

    def __post_init__(self):
        if self.range is not None and len(self.range) != 2:
            raise ValueError("Range must be a tuple of two values (min, max)")

        if self.range is not None:
            self.range = tuple(float(y) for y in self.range)

    def draw(self, ax):
        """Apply the y axis settings to the provided axes.

        Sets the axis label (including the unit, if given), the axis range and scale, and draws the ticks.

        :param ax: The matplotlib axes whose y axis is configured.
        :type ax: matplotlib.axes.Axes
        """
        if self.label is not None and self.unit is not None:
            ax.set_ylabel(f'{self.label} [{self.unit}]')
        elif self.label is not None:
            ax.set_ylabel(self.label)

        if self.range is not None:
            ax.set_ylim(self.range)
        ax.set_yscale(self.scale)
        self.ticks.draw(ax)

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`YAxis` from its keyword description.

        Recursively deserializes the nested ``ticks`` description into a :class:`YTicks` instance.

        :param kwargs: The y axis description.
        :returns: The instantiated y axis.
        :rtype: YAxis
        """
        _kwargs = _copy.deepcopy(kwargs)
        if 'ticks' in _kwargs:
            _kwargs['ticks'] = YTicks.from_dict(**_kwargs['ticks'])
        return Deserializable.make(cls, **_kwargs)


class Plot(ABC, Deserializable):
    r"""Base class for plots to be drawn into a figure."""

    def validate_semantics(self, context):
        yield from ()

    @abstractmethod
    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the plot for drawing.

        Subclasses override this method to prepare all of their items before drawing.

        :param context: The analysis file context used to resolve relative paths to data files.
            If ``None``, a default context rooted at the current working directory is used.
        :type context: AnalysisFileContext | None
        """
        raise NotImplementedError

    @abstractmethod
    def draw(self, ax):
        """Draw the plot on the provided axes.

        Subclasses override this method to render their items, axes, grid, and legend.

        :param ax: The matplotlib axes onto which the plot is drawn.
        :type ax: matplotlib.axes.Axes
        """
        raise NotImplementedError

    @abstractmethod
    def draw_watermark(self, ax, watermark):
        """Draw the watermark on the provided axes.

        :param ax: The matplotlib axes onto which the watermark is drawn.
        :type ax: matplotlib.axes.Axes
        :param watermark: The watermark to draw.
        :type watermark: eos.figure.common.Watermark
        """
        watermark.draw(ax)

@dataclass(kw_only=True)
class TwoDimensionalPlot(Plot):
    """Draws a 2D plot along a single set of axes.

    :param aspect: The aspect ratio of the plot. If None, the aspect ratio defaults to matplotlib's default.
    :type aspect: float | None
    :param grid: The grid properties of the plot. For the default properties, see :class:`Grid <eos.figure.Grid>`.
    :type grid: :class:`Grid <eos.figure.Grid>`
    :param items: A list of items to be drawn in the plot.
    :type items: list[:class:`Item <eos.figure.Item>`]
    :param legend: The legend properties of the plot. For the default properties, see :class:`Legend <eos.figure.Legend>`.
    :type legend: :class:`Legend <eos.figure.Legend>`
    :param title: The title of the plot. Defaults to None, which means no title is displayed.
    :type title: str | None
    :param xaxis: The x axis properties of the plot. For the default properties, see :class:`XAxis <eos.figure.XAxis>`.
    :type xaxis: :class:`XAxis <eos.figure.XAxis>`
    :param yaxis: The y axis properties of the plot. For the default properties, see :class:`YAxis <eos.figure.YAxis>`.
    :type yaxis: :class:`YAxis <eos.figure.YAxis>`
    """

    aspect:float|None=field(default=None)
    grid:Grid=field(default_factory=Grid)
    items:list[Item]
    legend:Legend=field(default=None)
    title:str=None
    xaxis:XAxis=field(default_factory=XAxis)
    yaxis:YAxis=field(default_factory=YAxis)

    def validate_semantics(self, context):
        for index, item in enumerate(self.items):
            if hasattr(item, 'validate_semantics'):
                yield from (
                    diagnostic.prefixed('items', index)
                    for diagnostic in item.validate_semantics(context)
                )

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
        """Prepare the plot by preparing each of its items.

        :param context: The analysis file context forwarded to each item's ``prepare`` method.
            If ``None``, a default context rooted at the current working directory is used.
        :type context: AnalysisFileContext | None
        """
        context = AnalysisFileContext() if context is None else context
        for item in self.items:
            item.prepare(context=context)

    def draw(self, ax):
        """Draw the plot and all of its items on the provided axes.

        Sets the title, configures the grid, aspect ratio, and axes, draws each item, and collects
        the items' legend entries into the legend.

        :param ax: The matplotlib axes onto which the plot is drawn.
        :type ax: matplotlib.axes.Axes
        """
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
        legend_entries = []
        for item in self.items:
            item.draw(ax)
            legend_entries.extend(item.legend())

        # Draw legend
        if self.legend is not None:
            self.legend.draw(ax=ax, entries=legend_entries)

    def draw_watermark(self, ax, watermark):
        """Draw the watermark on the provided axes.

        :param ax: The matplotlib axes onto which the watermark is drawn.
        :type ax: matplotlib.axes.Axes
        :param watermark: The watermark to draw.
        :type watermark: eos.figure.common.Watermark
        """
        watermark.draw(ax)

    @classmethod
    def from_dict(cls, **kwargs):
        """Create a :class:`TwoDimensionalPlot` from its keyword description.

        Recursively deserializes the nested ``grid``, ``legend``, ``xaxis``, and ``yaxis`` descriptions,
        as well as each entry of ``items`` via :class:`ItemFactory <eos.figure.item.ItemFactory>`.

        :param kwargs: The plot description. Must contain an ``items`` key.
        :returns: The instantiated plot.
        :rtype: TwoDimensionalPlot
        """
        _kwargs = _copy.deepcopy(kwargs)
        if 'grid' in _kwargs:
            _kwargs['grid'] = Grid.from_dict(**_kwargs['grid'])
        if 'legend' in _kwargs:
            _kwargs['legend'] = Legend.from_dict(**_kwargs['legend'])
        if 'xaxis' in _kwargs:
            _kwargs['xaxis'] = XAxis.from_dict(**_kwargs['xaxis'])
        if 'yaxis' in _kwargs:
            _kwargs['yaxis'] = YAxis.from_dict(**_kwargs['yaxis'])
        _kwargs['items'] = [ItemFactory.from_dict(**i) for i in _kwargs['items']]
        return Deserializable.make(cls, **_kwargs)


@dataclass
class EmptyPlot(Plot):
    """Draws an empty plot.

    Can be used in a grid as empty space instead of a plot. It does not accept any parameters.
    """

    _api_doc = inspect.cleandoc("""
    Drawing an Empty Plot
    ---------------------

    This plot's type is ``empty``. It produces an empty plot, which can be used in a grid as empty space.

    """)

    def __post_init__(self):
        pass

    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the empty plot for drawing.

        This plot requires no preparation.

        :param context: The analysis file context. Accepted for interface consistency and unused.
        :type context: AnalysisFileContext | None
        """
        pass

    def draw(self, ax):
        """Draw the empty plot by turning off the provided axes.

        :param ax: The matplotlib axes to turn off.
        :type ax: matplotlib.axes.Axes
        """
        ax.set_axis_off()

    def draw_watermark(self, ax, watermark):
        """Draw the watermark on the provided axes.

        Empty plots carry no watermark, so this is a no-op.

        :param ax: The matplotlib axes onto which the watermark would be drawn.
        :type ax: matplotlib.axes.Axes
        :param watermark: The watermark to draw.
        :type watermark: eos.figure.common.Watermark
        """
        pass

@dataclass(kw_only=True)
class OverviewSource(Deserializable):
    """One set of values entering an overview plot, drawn with a common style.

    :param type: The kind of source: ``'prediction'`` for values read from data files produced by
        the ``predict-observables`` task, or ``'constraint'`` for values taken from the EOS database
        of constraints.
    :type type: str
    :param label: The label of the source, used in the legend.
    :type label: str
    :param datafiles: The paths to the data files holding the predictions. Mandatory for, and
        permitted only with, ``type: 'prediction'``.
    :type datafiles: list[str]
    :param names: The names of the EOS constraints. Mandatory for, and permitted only with,
        ``type: 'constraint'``.
    :type names: list[str]
    :param color: The color used to draw the source. Defaults to the next color in the plot's cycle.
    :type color: str
    :param markerstyle: The style of the markers. Defaults to the next style in the plot's cycle.
    :type markerstyle: str
    :param linewidth: The width of the error bars. Defaults to 1.5.
    :type linewidth: float
    :param linestyle: The style of the line joining the markers. Defaults to 'none'.
    :type linestyle: str
    :param alpha: The opacity of the source. Defaults to 1.0.
    :type alpha: float
    """

    type: str
    label: str
    color: str = None
    markerstyle: str = None
    linewidth: float = None
    linestyle: str = None
    alpha: float = None

    datafiles: list[str] = None
    names: list[str] = None

    _TYPES = ('prediction', 'constraint')

    def __post_init__(self):
        if self.type not in self._TYPES:
            raise ValueError(f"source '{self.label}' has unknown type '{self.type}'; must be one of {self._TYPES}.")

        # a prediction is read from sample files, a constraint from the EOS database; providing the
        # other type's key means the source would be read from somewhere it does not describe
        required, forbidden = ('datafiles', 'names') if self.type == 'prediction' else ('names', 'datafiles')
        if not getattr(self, required):
            raise ValueError(f"source '{self.label}' of type '{self.type}' requires a non-empty '{required}'.")
        if getattr(self, forbidden):
            raise ValueError(f"source '{self.label}' of type '{self.type}' must not provide '{forbidden}'.")


@dataclass(kw_only=True)
class OverviewPlot(Plot):
    """Draws an overview of 1D predictions and measurements, one row per observable.

    Each observable occupies one row, and each source is drawn as a marker with error bars, the
    sources within one row staggered symmetrically about it. The plot's height should be chosen to
    suit the number of rows; as a rule of thumb, allow 0.6 inches per observable in the enclosing
    figure's ``size``, plus 1.2 inches for the axes and the legend.

    :param legend: The legend properties of the plot. Defaults to a legend drawn outside of the
        axes, above the plot. For the available properties, see :class:`Legend <eos.figure.Legend>`.
    :type legend: :class:`Legend <eos.figure.Legend>`
    :param xaxis: The x axis properties of the plot. Its ``range`` is mandatory, since the rows'
        separators span it. For the remaining properties, see :class:`XAxis <eos.figure.XAxis>`.
    :type xaxis: :class:`XAxis <eos.figure.XAxis>`
    :param yaxis: The y axis properties of the plot. Its range and ticks are determined by the
        observables and cannot be set. For the remaining properties, see
        :class:`YAxis <eos.figure.YAxis>`.
    :type yaxis: :class:`YAxis <eos.figure.YAxis>`
    :param observables: The observables entering the plot, each a single-entry mapping from the
        observable's name to its label. Nesting them into a list of lists draws a horizontal
        separator between the groups.
    :type observables: list[dict] | list[list[dict]]
    :param sources: The predictions and measurements to compare.
    :type sources: list[:class:`OverviewSource <eos.figure.OverviewSource>`]
    :param normalize: Whether to divide each observable's values by a reference value, and where to
        take that reference from: ``'constraint'`` uses the constraint source constraining that
        observable, ``'posterior'`` uses the source named by ``normalize_source``. Defaults to None,
        which draws the values themselves.
    :type normalize: str | None
    :param normalize_source: The label of the source to normalize to. Mandatory if
        ``normalize == 'posterior'``.
    :type normalize_source: str | None
    :param observable_offset: The spacing between the rows of two observables. Defaults to 1.0.
    :type observable_offset: float
    :param source_offset: The spacing between two sources within one row. Defaults to 0.05.
    :type source_offset: float
    """

    legend:Legend=field(default_factory=lambda: Legend(position='upper center', outside=True))
    xaxis:XAxis=field(default_factory=XAxis)
    yaxis:YAxis=field(default_factory=YAxis)
    observables:list=field(default=None)
    sources:list=field(default=None)

    normalize:str=field(default=None)
    normalize_source:str=field(default=None)
    observable_offset:float=field(default=1.0)
    source_offset:float=field(default=0.05)

    _NORMALIZE_MODES = (None, 'constraint', 'posterior')
    _DEFAULT_MARKERS = ('o', 's', '^', 'D', 'v', 'P', 'X')

    _api_doc = inspect.cleandoc("""
    Drawing an Overview of 1D Predictions and Measurements
    ------------------------------------------------------

    This plot's type is ``overview``. It compares predictions and measurements for a set of
    observables, one row per observable.

    The following keys are mandatory:

        * ``observables``: The observables entering the plot, each a single-entry mapping from the
          observable's name to its label. Nesting them into a list of lists draws a horizontal
          separator between the groups.
        * ``sources``: A list of objects of type :class:`eos.figure.OverviewSource`, which contains
          the predictions and measurements to compare.
        * ``xaxis``: An object of type :class:`eos.figure.XAxis`, whose ``range`` is mandatory here.

    The following keys are optional:

        * ``legend``: An object of type :class:`eos.figure.Legend`, which contains the legend
          properties. Defaults to a legend drawn outside of the axes, above the plot.
        * ``yaxis``: An object of type :class:`eos.figure.YAxis`, which contains the y axis
          properties. Its range and ticks are determined by the observables and cannot be set.
        * ``normalize``: One of ``'constraint'`` or ``'posterior'``, which normalizes each
          observable's values to a reference value.
        * ``normalize_source``: The label of the source to normalize to. Mandatory if
          ``normalize == 'posterior'``.
        * ``observable_offset``: The spacing between the rows of two observables. Defaults to 1.0.
        * ``source_offset``: The spacing between two sources within one row. Defaults to 0.05.

    """)

    def __post_init__(self):
        if not self.observables:
            raise ValueError("observables must include at least one item to be plotted.")
        if not self.sources:
            raise ValueError("sources must include at least one item to be plotted.")
        if self.xaxis.range is None:
            raise ValueError("xaxis must provide a range.")

        if self.normalize not in self._NORMALIZE_MODES:
            raise ValueError(f"normalize must be one of {self._NORMALIZE_MODES}, got {self.normalize!r}.")
        if self.normalize == 'posterior' and not self.normalize_source:
            raise ValueError("normalize == 'posterior' requires normalize_source "
                             "(the label of the source to normalize to).")

        if self.observable_offset < 1.0:
            raise ValueError("observable_offset must be at least 1.0")

        self._per_source_entries = None

    def validate_semantics(self, context):
        try:
            flat_obs, _, _ = self._group(self.observables)
        except (TypeError, ValueError) as error:
            yield Diagnostic(('observables',), Severity.ERROR, str(error))
            return

        for index, observable in enumerate(flat_obs):
            yield from _check_qualified(context, observable, 'observable', ('observables', index))

        for index, source in enumerate(self.sources):
            for name in source.names or []:
                yield from _check_qualified(context, name, 'constraint', ('sources', index, 'names'))

        if self.normalize == 'posterior' and not any(source.label == self.normalize_source for source in self.sources):
            yield Diagnostic(
                ('normalize_source',), Severity.ERROR,
                f"normalize_source '{self.normalize_source}' matches no source label"
            )

    @staticmethod
    def _group(entries:list):
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

    @staticmethod
    def _prediction_quantiles(samples, weights, level=68.27e-2):
        """Weighted median and asymmetric 1-sigma uncertainty widths as for `uncertainty` items."""
        half = level / 2.0
        interval = [0.5 - half, 0.5, 0.5 + half]
        lower, central, higher = _np.quantile(samples, q = interval, weights = weights, method='inverted_cdf', axis=0)
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
            if OverviewPlot._observable_matches(observable, key):
                return index
        return None

    def _prediction_entries_for_source(self, source, flat_obs, context):
        results = {}
        datafiles = source.datafiles or []
        if isinstance(datafiles, str):
            datafiles = [datafiles]

        predictions = [eos.data.Prediction(context.data_path(path)) for path in datafiles]

        for obs in flat_obs:
            for prediction in predictions:
                idx = self._prediction_observable_index(prediction, obs)
                if idx is None:
                    continue

                if obs in results:
                    raise ValueError(
                        f"source '{source.label}' provides observable '{obs}' "
                        f"from more than one data file."
                    )

                results[obs] = self._prediction_quantiles(
                    prediction.samples[:, idx],
                    prediction.weights
                )

            if obs not in results:
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
                "currently supported in OverviewPlot."
            )

        return entry

    def _constraint_entries_for_source(self, source, flat_obs):
        results = {}
        matched_constraints = {}

        for name in source.names:
            entry = self._constraint_entry(name)
            const = _yaml.safe_load(entry.serialize())
            raw_observable = const.get('observable', None)
            if not isinstance(raw_observable, str):
                raise ValueError(f"Constraint '{name}' does not provide a single observable.")
            nameobs, _, optionsobs = raw_observable.partition(';')
            if optionsobs == '':
                optionsobs = const.get('options', {})
                option_suffix =  ';' + ','.join(f'{k}={v}' for k, v in optionsobs.items()) if optionsobs else ''
            else:
                option_suffix = ';' + optionsobs
            # a serialized value such as '8e-05' carries no decimal point and so is a str, not a
            # float, under the YAML 1.1 resolver that safe_load applies
            mean = float(const['mean'])
            sigma_stat = const['sigma-stat']
            sigma_sys = const.get('sigma-sys', {'hi': 0.0, 'lo': 0.0})
            err_high = _np.sqrt(float(sigma_stat['hi']) ** 2 + float(sigma_sys['hi']) ** 2)
            err_low = _np.sqrt(float(sigma_stat['lo']) ** 2 + float(sigma_sys['lo']) ** 2)

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

    def _entries_per_source(self, flat_obs, context):
        """Each source's per-observable (central, err_low, err_high), read from disk only once."""
        if self._per_source_entries is None:
            entries = {}
            for source in self.sources:
                if source.type == 'prediction':
                    entries[id(source)] = self._prediction_entries_for_source(source, flat_obs, context)
                else:
                    entries[id(source)] = self._constraint_entries_for_source(source, flat_obs)
            self._per_source_entries = entries

        return self._per_source_entries

    def info_summary(self, context:AnalysisFileContext=None):
        """Print each observable's raw (un-normalized) central value and uncertainty for every
        source, as a sanity check before drawing (e.g. to catch an inverted or mismatched entry)."""
        context = AnalysisFileContext() if context is None else context
        flat_obs, _, _ = self._group(self.observables)

        per_source_entries = self._entries_per_source(flat_obs, context)

        headers = ['observable'] + [source.label for source in self.sources]
        rows = []
        for obs in flat_obs:
            row = [obs]
            for source in self.sources:
                entry = per_source_entries[id(source)].get(obs)
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

    def prepare(self, context:AnalysisFileContext=None):
        """Prepare the overview plot by reading every source and building its items.

        :param context: The analysis file context used to resolve relative paths to data files.
            If ``None``, a default context rooted at the current working directory is used.
        :type context: AnalysisFileContext | None
        """
        context = AnalysisFileContext() if context is None else context
        flat_obs, obs_group_ids, obs_labels = self._group(self.observables)
        n_obs = len(flat_obs)
        n_sources = len(self.sources)

        yvals = [1.0 + i * self.observable_offset for i in range(n_obs)]
        yvals.reverse()

        self._sourced_entries = self._build_sourced_entries(flat_obs, yvals, n_sources, context)

        # the rows are placed by this plot, so the y axis is labelled by the observables rather
        # than by the values of a quantity
        self.yaxis.range = (1.5 - self.observable_offset, n_obs * self.observable_offset + self.observable_offset)
        self.yaxis.ticks.minor = False
        self.yaxis.ticks.locations = yvals
        self.yaxis.ticks.labels = obs_labels

        # horizontal separators between observable groups, at the midpoint
        # between the last row of one group and the first row of the next
        separator_ys = [
            (yvals[i - 1] + yvals[i]) / 2.0
            for i in range(1, n_obs)
            if obs_group_ids[i] != obs_group_ids[i - 1]
        ]

        items = []
        if self.normalize:
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
        items += [
            {'type': 'expression', 'expression': f'{y}', 'color': 'lightgray', 'range': self.xaxis.range}
            for y in separator_ys
        ]

        self._items = [ItemFactory.from_dict(**item) for item in items]
        for item in self._items:
            item.prepare(context=context)

        self.info_summary(context)

    def draw(self, ax):
        """Draw the overview plot and all of its items on the provided axes.

        :param ax: The matplotlib axes onto which the plot is drawn.
        :type ax: matplotlib.axes.Axes
        """
        # Remove default margin used by matplotlib
        ax.margins(0.0)

        # Handle axes
        self.xaxis.draw(ax)
        self.yaxis.draw(ax)

        # Draw all items
        legend_entries = []
        for item in self._items:
            item.draw(ax)
            legend_entries.extend(item.legend())

        # Draw legend
        if self.legend is not None:
            self.legend.draw(ax=ax, entries=legend_entries)

    def draw_watermark(self, ax, watermark):
        """Draw the watermark on the provided axes.

        :param ax: The matplotlib axes onto which the watermark is drawn.
        :type ax: matplotlib.axes.Axes
        :param watermark: The watermark to draw.
        :type watermark: eos.figure.common.Watermark
        """
        watermark.draw(ax)

    def _build_sourced_entries(self, flat_obs, yvals, n_sources, context):
        """Build one drawable entry per source, aggregating across every observable it covers."""
        default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

        per_source_entries = self._entries_per_source(flat_obs, context)

        # resolve the normalization reference, if any
        reference = {}
        if self.normalize == 'constraint':
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
        elif self.normalize == 'posterior':
            ref_source = next((s for s in self.sources if s.label == self.normalize_source), None)
            if ref_source is None:
                raise ValueError(f"normalize_source '{self.normalize_source}' matches no source label.")
            reference = per_source_entries[id(ref_source)]

        # symmetric stagger of sources around each observable's row
        stagger = [(-(n_sources - 1) / 2.0 + i) * self.source_offset for i in range(n_sources)]

        entries = []
        for source_index, source in enumerate(self.sources):
            source_data = per_source_entries[id(source)]
            positions, xerrors = [], []

            for obs_index, obs in enumerate(flat_obs):
                if obs not in source_data:
                    continue

                central, err_low, err_high = source_data[obs]

                if self.normalize and obs in reference:
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
                'color': source.color or default_colors[source_index % len(default_colors)],
                'marker': source.markerstyle or self._DEFAULT_MARKERS[source_index % len(self._DEFAULT_MARKERS)],
                'linestyle': source.linestyle or 'none',
                'alpha': source.alpha or 1.0,
                'linewidth': source.linewidth or 1.5,
                'label': source.label,
            })

        return entries

    @classmethod
    def from_dict(cls, **kwargs):
        """Create an :class:`OverviewPlot` from its keyword description.

        Recursively deserializes the nested ``legend``, ``xaxis``, and ``yaxis`` descriptions, as
        well as each entry of ``sources``.

        :param kwargs: The plot description. Must contain ``observables``, ``sources``, and ``xaxis`` keys.
        :returns: The instantiated plot.
        :rtype: OverviewPlot
        """
        _kwargs = _copy.deepcopy(kwargs)
        if 'legend' in _kwargs:
            _kwargs['legend'] = Legend.from_dict(**_kwargs['legend'])
        if 'xaxis' in _kwargs:
            _kwargs['xaxis'] = XAxis.from_dict(**_kwargs['xaxis'])
        if 'yaxis' in _kwargs:
            _kwargs['yaxis'] = YAxis.from_dict(**_kwargs['yaxis'])
        _kwargs['sources'] = [OverviewSource.from_dict(**source) for source in _kwargs.get('sources') or []]
        return Deserializable.make(cls, **_kwargs)


class PlotFactory:
    """Factory that creates :class:`Plot` instances from their YAML or dictionary description.

    The concrete plot class is selected from the optional ``type`` key using the :attr:`registry`,
    which maps each supported type string to its corresponding :class:`Plot` subclass. If no ``type``
    is given, the default ``'2D'`` plot is created.
    """

    # Also build the documentation based on ordered registry
    # Initializer is well-defined for python version >= 3.6
    registry = {
        '2D':       TwoDimensionalPlot, # default
        'empty':    EmptyPlot,
        'overview': OverviewPlot
    }

    @staticmethod
    def from_yaml(yaml_data:str):
        """Create a plot from a YAML description.

        :param yaml_data: A YAML string describing a single plot, optionally including its ``type`` key.
        :type yaml_data: str
        :returns: The instantiated plot.
        :rtype: Plot
        """
        kwargs = _yaml.safe_load(yaml_data)
        return PlotFactory.from_dict(**kwargs)

    @staticmethod
    def from_dict(**kwargs):
        """Create a plot from its keyword description.

        The optional ``type`` key selects the concrete :class:`Plot` subclass from :attr:`registry`,
        defaulting to ``'2D'``; the remaining keyword arguments are forwarded to that subclass. The
        item color cycle is reset before the plot is created.

        :param kwargs: The plot description. May contain a ``type`` key identifying a registered plot type.
        :returns: The instantiated plot.
        :rtype: Plot
        :raises ValueError: If the ``type`` key names an unknown plot type.
        """

        if 'type' not in kwargs:
            plot_type = '2D'
        else:
            plot_type = kwargs.pop('type')

        if plot_type not in PlotFactory.registry:
            raise ValueError(f'Unknown plot type: {plot_type}')

        ItemColorCycler.reset()

        return PlotFactory.registry[plot_type].from_dict(**kwargs)
