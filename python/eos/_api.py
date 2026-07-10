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

"""
Single source of truth for the public, C++ (boost::python) binding classes of EOS that are
part of the documented Python API.

This module intentionally has no dependencies (in particular it does not import ``eos``), so that
it can be imported cheaply by both:

  * ``eos/_api_TEST.py``, which asserts that every method bound on these classes names its
    arguments (no ``argN`` placeholder) and carries a docstring; and
  * ``doc/reference/python.rst.py``, which generates the ``.. autoclass::`` directives for the
    Python API reference.

The classes are grouped exactly as in the Python API reference: ``API_BASIC_CLASSES`` are the
"Basic Classes" and ``API_COMMON_CLASSES`` are the "Common Classes". Each dict maps the attribute
name of a class in the ``eos`` module to whether its ``.. autoclass::`` directive uses
``:inherited-members:``. The flag is ``True`` for the thin Python wrapper classes (e.g.
``eos.Parameters``) whose bound methods are inherited from an underscore-prefixed C++ base (e.g.
``_eos._Parameters``) and must therefore be documented through inheritance; it is ``False`` for
classes that are documented with ``:members:`` only. The insertion order is the order in which the
classes are documented.
"""

API_BASIC_CLASSES = {
    'Constraint':         False,
    'ConstraintEntry':    False,
    'Constraints':        True,
    'Kinematics':         True,
    'LogLikelihood':      False,
    'LogLikelihoodBlock': False,
    'LogPrior':           False,
    'LogPosterior':       False,
    'Model':              False,
    'Observable':         False,
    'ObservableEntry':    False,
    'Observables':        True,
    'Options':            False,
    'Parameter':          False,
    'Parameters':         True,
    'QualifiedName':      False,
    'Reference':          False,
    'References':         True,
    'SignalPDF':          True,
    'SignalPDFEntry':     False,
    'SignalPDFs':         True,
    'Unit':               False,
}

API_COMMON_CLASSES = {
    'Analysis':        False,
    'AnalysisFile':    False,
    'BestFitPoint':    False,
    'GoodnessOfFit':   False,
    'ObservableCache': False,
}
