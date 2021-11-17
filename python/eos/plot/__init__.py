# Copyright (c) 2018 Frederik Beaujean
# Copyright (c) 2018 Danny van Dyk
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

## minimal versions checks
# matplotlib
import matplotlib._version
try:
    if matplotlib._version.get_versions()['version'] < '2.0':
        raise ImportError('eos.plot requires matplotlib in version 2.0 or higher')
except AttributeError:
    # matplotlib._version.get_versions() was removed in version 3.5 or higher
    pass

from . import config

from .plotter import *
