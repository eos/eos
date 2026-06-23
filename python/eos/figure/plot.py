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

from .item import Item, ItemFactory, ItemColorCycler

import copy as _copy
import eos
import inspect
import matplotlib
import matplotlib.pyplot as plt
import yaml as _yaml


@dataclass
class Legend(Deserializable):
    r"""Represents the legend properties in a plot.

    :param position: The position of the legend in the plot. Defaults to 'best', which automatically places the legend in the best location.
    :type position: str
    """

    position:str=field(default='best')

    def __post_init__(self):
        pass

    def draw(self, ax, entries=None):
        """Draw the legend on the axes.

        :param ax: The axes to draw the legend on.
        :type ax: matplotlib.axes.Axes
        :param entries: A list of tuples (handle, label) for the legend entries. If None, uses the entries from the items.
        :type entries: list[tuple[matplotlib.artist.Artist, str]] | None
        """
        if entries is not None:
            handles = [entry[0] for entry in entries]
            labels = [entry[1] for entry in entries]
            ax.legend(loc=self.position, handles=handles, labels=labels)
        else:
            ax.legend(loc=self.position)


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
    """

    minor:bool=field(default=True)
    position:str=field(default='bottom')
    visible:bool=field(default=True)
    format:str|None=field(default=None)

    def __post_init__(self):
        POSITIONS = ['bottom', 'top', 'both']
        if self.position not in POSITIONS:
            raise ValueError(f'Unknown position: {self.position}. Must be one of {",".join(POSITIONS)}')

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
            if not log_scale:
                ax.xaxis.set_major_locator(matplotlib.ticker.AutoLocator())
            else:
                ax.xaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10.0))
            if self.minor:
                if not log_scale:
                    ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
                else:
                    ax.xaxis.set_minor_locator(matplotlib.ticker.LogLocator(base=10.0, subs='auto'))
            if self.format is not None:
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
    """

    minor:bool=field(default=True)
    position:str=field(default='left')
    visible:bool=field(default=True)
    format:str|None=field(default=None)

    def __post_init__(self):
        POSITIONS = ['left', 'right', 'both']
        if self.position not in POSITIONS:
            raise ValueError(f'Unknown position: {self.position}. Must be one of {",".join(POSITIONS)}')

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
            if not log_scale:
                ax.yaxis.set_major_locator(matplotlib.ticker.AutoLocator())
            else:
                ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10.0))
            if self.minor:
                if not log_scale:
                    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
                else:
                    ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(base=10.0, subs='auto'))
            if self.format is not None:
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

class PlotFactory:
    """Factory that creates :class:`Plot` instances from their YAML or dictionary description.

    The concrete plot class is selected from the optional ``type`` key using the :attr:`registry`,
    which maps each supported type string to its corresponding :class:`Plot` subclass. If no ``type``
    is given, the default ``'2D'`` plot is created.
    """

    # Also build the documentation based on ordered registry
    # Initializer is well-defined for python version >= 3.6
    registry = {
        '2D':    TwoDimensionalPlot, # default
        'empty': EmptyPlot
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
