#!/usr/bin/python
# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2026 Danny van Dyk
# Copyright (c) 2026 Lorenz Gärtner
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

import tempfile
import textwrap
import unittest

from eos.diagnostic import Diagnostic, Severity, attach_line_numbers


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

    def test_str_with_line(self):
        diagnostic = Diagnostic(('priors', 0), Severity.ERROR, 'invalid prior', line=3)
        self.assertEqual(str(diagnostic), 'line 3: priors[0]: invalid prior')

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


class AttachLineNumbersTests(unittest.TestCase):

    def setUp(self):
        # line 1 is blank; line numbers below refer to this file as written.
        yaml_text = textwrap.dedent('''
            steps:
              - id: step-a
                title: A
              - title: B
            priors:
              - name: prior-a
                type: uniform
              - name: bad/name
                id: prior-b
                type: uniform
            observables:
              clean-name:
                latex: '$x$'
              bad/name:
                latex: '$y$'
        ''')
        self.file = tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False)
        self.file.write(yaml_text)
        self.file.close()

    def _line(self, path):
        [diagnostic] = attach_line_numbers([Diagnostic(path, Severity.ERROR, 'message')], self.file.name)
        return diagnostic.line

    def test_top_level_key(self):
        self.assertEqual(self._line(('steps',)), 2)

    def test_id_identified_sequence_element(self):
        self.assertEqual(self._line(('steps', 'step-a', 'title')), 4)

    def test_unnamed_sequence_element_falls_back_to_index(self):
        self.assertEqual(self._line(('steps', 1, 'title')), 5)

    def test_named_sequence_element(self):
        self.assertEqual(self._line(('priors', 'prior-a', 'type')), 8)

    def test_missing_key_falls_back_to_enclosing_mapping(self):
        self.assertEqual(self._line(('priors', 'prior-a', 'min')), 7)

    def test_unresolvable_path_falls_back_to_document_start(self):
        self.assertEqual(self._line(('nonexistent',)), 2)

    def test_sequence_identifier_does_not_leak_across_sections(self):
        # 'priors' always uses 'name' (the 'id' identifier is steps-only, analysis_file_description.py:1651):
        # a dirty 'name' must fall back to the index even though this entry also has a clean 'id' field.
        self.assertEqual(self._line(('priors', 1, 'type')), 11)

    def test_name_keyed_mapping_section_uses_clean_key_as_segment(self):
        self.assertEqual(self._line(('observables', 'clean-name', 'latex')), 14)

    def test_name_keyed_mapping_section_falls_back_to_index_for_dirty_key(self):
        # 'bad/name' is not usable as a path segment (analysis_file_description.py's _segments()
        # falls back to the entry's position), so the index key is the integer 1, not 'bad/name'.
        self.assertEqual(self._line(('observables', 1, 'latex')), 16)

    def test_empty_diagnostics_skips_file_access(self):
        self.assertEqual(attach_line_numbers([], '/nonexistent/file.yaml'), [])


if __name__ == '__main__':
    unittest.main()
