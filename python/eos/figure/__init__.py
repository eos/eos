# Copyright (c) 2018 Frederik Beaujean
# Copyright (c) 2018 Danny van Dyk
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

"""
The `eos.figure` submodule provides means to produce publication-quality figures.

A figure contains one or more plots, each of which is a set of axis that can be
thought of as a 'canvas' to draw on.
A plot further contains items for display, for example a constraint in the form of
data points with error bars or an observable as a function of a kinematic variable
in the form or a line graph.

"""

#__all__ = ['make']

## minimal versions checks
# matplotlib
import matplotlib._version
try:
    if matplotlib._version.get_versions()['version'] < '2.0':
        raise ImportError('eos.figure requires matplotlib in version 2.0 or higher')
except AttributeError:
    # matplotlib._version.get_versions() was removed in version 3.5 or higher
    pass

from . import config
from .figure import *
from .plot import *
from .item import *
from .common import Watermark
from .data import DataFile

## Inject EOS color names
import matplotlib.colors as mcolors
for i, cspec in enumerate(ItemColorCycler._colors):
    mcolors._colors_full_map[f"eos:C{i}"] = cspec
