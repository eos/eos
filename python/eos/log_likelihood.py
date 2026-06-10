# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2026 Danny van Dyk
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

from _eos import LogLikelihoodBlock
import numpy as np


def _unbinned_1d(cache, pdf_name, kinematics, options, resolution, observations):
    """
    Create a new unbinned log-likelihood block with a resolution function.

    The named SignalPDF is evaluated on the supplied grid and convolved with the resolution
    function via a discrete Fourier transform, yielding the resolution-smeared PDF sampled on
    the grid. The log-likelihood is then the sum over the observed events of the logarithm of
    the smeared PDF, evaluated at each observation by linear interpolation from the surrounding
    grid points.

    :param cache: The observable cache used by the total log-likelihood.
    :type cache: eos.ObservableCache
    :param pdf_name: The name of the SignalPDF to evaluate at each grid point.
    :type pdf_name: eos.QualifiedName
    :param kinematics: One Kinematics object per grid point.
    :type kinematics: list of eos.Kinematics
    :param options: Options forwarded to the SignalPDF constructor.
    :type options: eos.Options
    :param resolution: Discretised resolution function on the same grid, in natural (centred) order,
        i.e. the zero offset (the peak of a symmetric kernel) sits at the centre of the array.
        The kernel is converted to the wrap-around order expected by the underlying circular
        convolution via :func:`numpy.fft.ifftshift` before the block is constructed.
    :type resolution: list of float
    :param observations: The observed events, expressed as kinematic variables.
    :type observations: list of eos.Kinematics

    :returns: The new block.
    :rtype: eos.LogLikelihoodBlock

    .. note::
        The convolution underlying the unbinned likelihood is circular. The grid must therefore be
        padded with a sufficiently large region in which both the PDF and the resolution function
        are negligible, so that no appreciable density wraps across the boundary.
    """
    # Convert the resolution from natural (centred) order to the wrap-around ("FFT-native") order
    # expected by the underlying binding: ifftshift moves the central (zero-offset) sample to index 0.
    resolution = np.fft.ifftshift(np.asarray(resolution, dtype=float))

    return LogLikelihoodBlock._Unbinned1D(cache, pdf_name, kinematics, options, resolution.tolist(), observations)


# Expose the wrapper as the public factory method on the native LogLikelihoodBlock class.
LogLikelihoodBlock.Unbinned1D = staticmethod(_unbinned_1d)
