# Copyright (c) 2025 Danny van Dyk
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

from dataclasses import dataclass, field
from eos.deserializable import Deserializable
from _eos import __version__ as eos_version

import numpy as np


@dataclass(kw_only=True)
class Watermark(Deserializable):
    r"""Inserts an EOS watermark into a figure.

    :param position: The position of the watermark as ``'<vpos> <hpos>'``.
        The vertical position can be any of ``'right'``, ``'left'`` and ``'center'``.
        The horizontal position can be any of ``'upper'``, ``'lower'``, and ``'center'``.
    :type position: str
    :param preliminary: If true, marks the figure as a preliminary version.
    :type preliminary: bool
    """
    position:str=field(default='upper right')
    preliminary:bool=field(default=False)

    def __post_init__(self):
        xdelta, ydelta = (0.04, 0.04)
        vpos, hpos = self.position.split(' ')

        if hpos == 'right':
            self._x = 1 - xdelta
        elif hpos == 'left':
            self._x = xdelta
        elif hpos == 'center':
            self._x = 0.5
        else:
            raise ValueError(f'invalid horizontal position \'{hpos}\'')
        self._halign = hpos

        if vpos == 'lower':
            self._y = 0 + ydelta
            self._valign = 'bottom'
        elif vpos == 'upper':
            self._y = 1 - ydelta
            self._valign = 'top'
        else:
            raise ValueError(f'invalid vertical position \'{vpos}\'')

        if self.preliminary:
            self._color = 'red'
            self._version = 'Preliminary'
        else:
            self._color = 'OrangeRed'
            self._version = f'v{eos_version}'

    def draw(self, ax):
        ax.text(self._x, self._y, fr'\textsf{{\textbf{{EOS {self._version}}}}}',
                transform=ax.transAxes,
                color=self._color, alpha=0.5, bbox=dict(facecolor='white', alpha=0.5, lw=0),
                horizontalalignment=self._halign, verticalalignment=self._valign, zorder=+5)
