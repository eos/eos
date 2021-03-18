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

import matplotlib
from matplotlib import rcParams

try:
    if __IPYTHON__:
        pass
except NameError as e:
    matplotlib.use('pgf')

# set some default values for plotting
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = 'Computer Modern Sans serif'
matplotlib.rcParams['font.size'] = 14
matplotlib.rcParams['font.weight'] = 'normal'

matplotlib.rcParams['axes.labelsize'] = 16
matplotlib.rcParams['axes.linewidth'] = 1
matplotlib.rcParams['axes.titlepad'] = 12
matplotlib.rcParams['axes.formatter.use_mathtext'] = True

matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['savefig.pad_inches'] = 0.1

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

matplotlib.rcParams['text.usetex'] = True

matplotlib.rcParams['pgf.preamble'] = r'''
\usepackage{amsmath}
\usepackage{xcolor}
'''

matplotlib.rcParams['errorbar.capsize'] = 5
