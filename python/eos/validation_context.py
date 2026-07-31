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

from collections import Counter
from types import MappingProxyType

from .analysis_file_description import MaskExpressionComponent


class ValidationContext:
    """Read-only entity lookup context for one analysis file.

    The context snapshots the names declared by an
    :class:`~eos.analysis_file_description.AnalysisFileDescription`. Lookups first
    consult this per-file shadow registry and then, for EOS entities, the
    corresponding real registry. Constructing and querying a context never inserts
    into or otherwise changes an EOS registry.

    :param description: The deserialized analysis file whose declarations should
        be visible during semantic validation.
    :type description: eos.analysis_file_description.AnalysisFileDescription
    """

    _LOCAL_KINDS = frozenset(('prior', 'likelihood', 'mask', 'posterior'))

    def __init__(self, description):
        import eos

        observables = {observable.name for observable in description.observables}
        for mask in description.masks:
            observables.update(
                item.name for item in mask.description
                if isinstance(item, MaskExpressionComponent)
            )

        parameters = {parameter.name for parameter in description.parameters}
        parameter_aliases = {}
        for parameter in description.parameters:
            parameters.update(parameter.alias_of)
            for alias in parameter.alias_of:
                parameter_aliases.setdefault(alias, set()).add(parameter.name)

        constraints = set()
        likelihood_constraints = {}
        for likelihood in description.likelihoods:
            names = {constraint.name for constraint in likelihood.manual_constraints}
            constraints.update(names)
            if names:
                likelihood_constraints[likelihood.name] = names

        self._shadow = MappingProxyType({
            'observable': frozenset(observables),
            'parameter': frozenset(parameters),
            'constraint': frozenset(constraints),
            'prior': frozenset(prior.name for prior in description.priors),
            'likelihood': frozenset(likelihood.name for likelihood in description.likelihoods),
            'mask': frozenset(mask.name for mask in description.masks),
            'posterior': frozenset(posterior.name for posterior in description.posteriors),
        })
        self._counts = {
            kind: Counter({name: 0 for name in names})
            for kind, names in self._shadow.items()
        }
        self._parameter_aliases = MappingProxyType({
            alias: frozenset(declarations)
            for alias, declarations in parameter_aliases.items()
        })
        self._likelihood_constraints = MappingProxyType({
            likelihood: frozenset(names)
            for likelihood, names in likelihood_constraints.items()
        })
        self._real_registries = MappingProxyType({
            'observable': eos.Observables(),
            'parameter': eos.Parameters(),
            'constraint': eos.Constraints(),
        })

    def lookup(self, kind, qn):
        """Return whether *qn* exists in the file shadow or the real EOS registry.

        ``kind`` is one of ``observable``, ``parameter``, ``constraint``,
        ``prior``, ``likelihood``, ``mask``, or ``posterior``. The final four
        namespaces are local to an analysis file and therefore have no real-registry
        fallback.
        """
        if kind not in self._shadow:
            raise ValueError(f'Unknown validation entity kind: {kind!r}')

        name = str(qn)
        if name in self._shadow[kind]:
            self._counts[kind][name] += 1
            if kind == 'parameter':
                for declaration in self._parameter_aliases.get(name, ()):
                    self._counts[kind][declaration] += 1
            elif kind == 'likelihood':
                # A manual constraint is never referenced by name: it is used by virtue of being
                # declared inside a likelihood, whose manual_constraints AnalysisFile.analysis()
                # passes straight to eos.Analysis. Credit a likelihood's manual constraints along
                # with the likelihood, so that only those in an unused likelihood are reported --
                # and that likelihood is itself reported alongside them.
                for constraint in self._likelihood_constraints.get(name, ()):
                    self._counts['constraint'][constraint] += 1
            return True

        if kind in self._LOCAL_KINDS:
            return False

        import eos

        qn = eos.QualifiedName(name)
        if kind == 'parameter':
            found = self._real_registries[kind].has(qn)
        elif kind == 'observable':
            found = qn in self._real_registries[kind]
        elif kind == 'constraint':
            try:
                self._real_registries[kind][qn]
                found = True
            except RuntimeError:
                found = False
        else:
            raise AssertionError(f'Unhandled validation entity kind: {kind!r}')

        if found:
            self._counts[kind][name] += 1
            return True

        # EOS resolves an observable name through Observable::make, which synthesises a parameter
        # observable whenever the name names a parameter rather than a registered observable. A
        # manual constraint may therefore constrain 'decay-constant::K_u' directly, and an
        # expression may reference it as <<decay-constant::K_u>>. Such a name is a valid observable
        # name, so fall back to the parameter namespace before reporting it as unknown; the
        # reference is then credited to the parameter, which is what it actually uses.
        if kind == 'observable':
            return self.lookup('parameter', qn)

        return False

    def unused(self, kind):
        """Return the names of shadow entities of *kind* that have not been looked up."""
        if kind not in self._shadow:
            raise ValueError(f'Unknown validation entity kind: {kind!r}')

        return frozenset(
            name for name in self._shadow[kind]
            if self._counts[kind][name] == 0
        )
