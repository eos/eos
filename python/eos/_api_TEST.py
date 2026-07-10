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

import re
import unittest

import eos
from eos._api import API_BASIC_CLASSES, API_COMMON_CLASSES

# The set of documented classes, mapping each class' attribute name in the ``eos`` module to whether
# its ``.. autoclass::`` directive uses ``:inherited-members:`` (see eos._api).
_DOCUMENTED_CLASSES = {**API_BASIC_CLASSES, **API_COMMON_CLASSES}

# For each documented class, every method that boost::python binds must:
#
#   (a) name all of its arguments (including the receiver ``self``) via boost::python::args,
#       so that no ``argN`` placeholder leaks into the rendered signature; and
#   (b) carry a user-provided docstring.
#
# This test introspects the live bindings rather than the rendered HTML, so it is
# independent of the documentation theme or Sphinx version. The set of classes and the
# ``:inherited-members:`` flag are shared with the documentation via eos._api.

# Matches an unnamed boost::python argument placeholder in a signature, e.g. ``(Parameter)arg1``.
_ARGN = re.compile(r'\)arg\d')

# Matches the signature of a boost::python iterator method bound via ``range(...)``, e.g.
# ``sections( (_Parameters)arg1) -> object``. Such methods take only the receiver and return a
# Python iterator; boost::python does not accept ``args(...)`` for them, so their receiver
# unavoidably renders as ``arg1``. They expose no user-facing arguments and are therefore exempt
# from the unnamed-argument check. (The trailing `` :`` appears when the method has a docstring.)
_ITERATOR_SIG = re.compile(r'\w+\(\s*\(\w+\)arg1\s*\)\s*->\s*object\b')


def _is_iterator_method(doc):
    """Return whether ``doc``'s signature is that of a ``range(...)``-bound iterator method."""
    first = doc.strip().split('\n', 1)[0]
    return bool(_ITERATOR_SIG.match(first))


def _is_boost_function(obj):
    """Return whether ``obj`` is a live boost::python binding (as opposed to a Python wrapper)."""
    t = type(obj)
    return t.__module__ == 'Boost.Python' and t.__name__ == 'function'


def _bound_methods(cls, include_inherited):
    """
    Yield ``(name, doc)`` for every boost::python method rendered for ``cls``.

    When ``include_inherited`` is set, the class' full MRO is walked so that C++ methods
    inherited by a Python wrapper class (e.g. ``eos.Parameters`` deriving from
    ``_eos._Parameters``) are covered; this mirrors the ``:inherited-members:`` autodoc option.
    Otherwise only the class' own methods are considered.

    Methods that a Python wrapper overrides (e.g. ``eos.SignalPDF.make``) resolve to the Python
    function and are therefore skipped: they are documented on the Python side, not in the
    bindings. Dunder methods are excluded.
    """
    bases = cls.__mro__ if include_inherited else [cls]
    seen = set()
    for base in bases:
        if base is object:
            continue
        for name in vars(base):
            if name.startswith('__') or name in seen:
                continue
            seen.add(name)
            try:
                attr = getattr(cls, name)
            except Exception:
                continue
            if not _is_boost_function(attr):
                continue
            yield name, (attr.__doc__ or '')


def _signature_only(name, doc):
    """Return whether ``doc`` contains only boost's auto-generated signature line(s), i.e. no
    user-provided description follows."""
    for line in doc.split('\n'):
        stripped = line.strip()
        if not stripped:
            continue
        # boost renders one signature line per overload, each starting with ``name(``.
        if stripped.startswith(name + '('):
            continue
        return False
    return True


def _is_zero_argument_static(name, doc):
    """Return whether ``doc`` describes a receiver-less, argument-less binding, e.g. ``GeV() -> Unit``."""
    return bool(re.match(r'\s*' + re.escape(name) + r'\(\s*\)\s*->', doc.strip()))


class BindingDocstringTests(unittest.TestCase):

    def test_no_unnamed_arguments(self):
        "every documented binding names all of its arguments (no ``argN`` placeholder leaks)"

        offenders = []
        for clsname, include_inherited in _DOCUMENTED_CLASSES.items():
            cls = getattr(eos, clsname)
            for name, doc in _bound_methods(cls, include_inherited):
                # iterator methods bound via boost::python::range cannot carry args(...); their
                # receiver necessarily renders as ``arg1`` and they take no user arguments.
                if _is_iterator_method(doc):
                    continue
                if _ARGN.search(doc):
                    offenders.append(f'{clsname}.{name}')

        self.assertEqual(
            offenders, [],
            "the following bindings leak an unnamed 'argN' placeholder into their signature "
            "(add boost::python::args(\"self\", ...) to their .def in python/_eos.cc):\n  "
            + '\n  '.join(offenders)
        )

    def test_methods_have_docstrings(self):
        "every documented binding carries a user-provided docstring"

        offenders = []
        for clsname, include_inherited in _DOCUMENTED_CLASSES.items():
            cls = getattr(eos, clsname)
            for name, doc in _bound_methods(cls, include_inherited):
                # eos.Unit's zero-argument static factory methods (GeV(), GeV2(), ...) are
                # self-explanatory and intentionally left un-docstringed; the eos.Unit class
                # docstring enumerates them. See the Group B documentation decision.
                if clsname == 'Unit' and _is_zero_argument_static(name, doc):
                    continue
                if _signature_only(name, doc):
                    offenders.append(f'{clsname}.{name}')

        self.assertEqual(
            offenders, [],
            "the following bindings lack a docstring (add one to their .def in python/_eos.cc):\n  "
            + '\n  '.join(offenders)
        )


if __name__ == '__main__':
    unittest.main(verbosity=5)
