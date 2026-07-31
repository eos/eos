#!/usr/bin/python
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

import unittest

from eos.diagnostic import Diagnostic, Severity


class DiagnosticTests(unittest.TestCase):

    def test_location(self):
        cases = (
            ((), ''),
            ((0, 'min'), '[0]/min'),
            (('items', 0, 1, 'min'), 'items[0][1]/min'),
            (('likelihoods', 'LLH2', 'constraints', 0), 'likelihoods/LLH2/constraints[0]'),
            (('priors', 'FF', 'descriptions', 2, 'min'), 'priors/FF/descriptions[2]/min'),
            (('observable', 'B_c->J/psi::FormFactors@HPQCD:2025A'), 'observable/B_c->J/psi::FormFactors@HPQCD:2025A'),
        )

        for path, expected in cases:
            with self.subTest(path=path):
                diagnostic = Diagnostic(path, Severity.ERROR, 'message')
                self.assertEqual(diagnostic.location(), expected)

    def test_prefixed(self):
        original = Diagnostic(('min',), Severity.WARNING, 'message')
        prefixed = original.prefixed('priors', 'FF', 'descriptions', 2)

        self.assertEqual(original.path, ('min',))
        self.assertEqual(prefixed.path, ('priors', 'FF', 'descriptions', 2, 'min'))
        self.assertIsNot(original, prefixed)

    def test_str(self):
        diagnostic = Diagnostic(('priors', 0), Severity.ERROR, 'invalid prior')
        self.assertEqual(str(diagnostic), 'priors[0]: invalid prior')

    def test_sorting_with_mixed_path_segments(self):
        diagnostics = [
            Diagnostic((0,), Severity.WARNING, 'index'),
            Diagnostic(('priors',), Severity.WARNING, 'warning'),
            Diagnostic(('priors',), Severity.ERROR, 'error'),
        ]

        self.assertEqual(
            sorted(diagnostics),
            [diagnostics[2], diagnostics[1], diagnostics[0]],
        )


if __name__ == '__main__':
    unittest.main()
