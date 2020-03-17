# Copyright (c) 2017-2020 Danny van Dyk
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

from _eos import _Parameters, QualifiedName
import re

class Parameters(_Parameters):
    """
    Represents the set of parameters known to EOS.

    An independent instance of the default parameters can be obtained from the :meth:`Defaults <eos.Parameters.Defaults>` method.
    See also `the complete list of parameters <../parameters.html>`_.
    """
    def __new__(cls, *args, **kwargs):
        instance = _Parameters.Defaults()
        instance.__class__ = Parameters
        return(instance)

    def __init__(self, name=None):
        self.name = name
        self._replacements = [
            (re.compile(r'(\\GeV)'),      r'\\text{GeV}'),
            (re.compile(r'\$([^\$]*)\$'), r'$$\1$$'),        # ensure display mode MathJax presentation
        ]

    def filter_entry(self, name):
        if self.name and not self.name in name:
            return False

        return True

    @staticmethod
    def _key(p):
        n = p.name()
        try:
            n = QualifiedName(n)
        except RuntimeError:
            pass

        if type(n) is str:
            return(('', '', n))
        else:
            return((str(n.prefix_part()), str(n.suffix_part()), str(n.name_part())))

    def _latex_refine(self, s):
        result = s
        for regexp, repl in self._replacements:
            result = regexp.sub(repl, result)

        return(result)

    def _repr_html_(self):
        result = '<table>\n'
        for section in self.sections():
            section_entries = 0
            section_result = '  <tr><th style="text-align:left" colspan=3><big>{section}</big></th></tr>\n'.format(section=section.name())
            for group in section:
                group_entries = 0
                group_result = '    <tr><th style="text-align:left" colspan=3>{group}</th></tr>\n'.format(group=group.name())

                group_parameters = [p for p in group]
                #group_parameters.sort(key = lambda p: p.name())
                group_parameters.sort(key = Parameters._key)
                for parameter in group_parameters:
                    name  = parameter.name()
                    if not self.filter_entry(name):
                        continue

                    latex = self._latex_refine(parameter.latex())
                    if 0 == len(latex):
                        latex = '---'

                    value = parameter.evaluate()

                    group_result += r'      <tr><th><tt style="color:grey">{name}</tt></th><td style="text-align:left">{latex}</td><td style="text-align:right">{value}</td></tr>'.format(name=name,latex=latex,value=value)
                    group_entries += 1

                group_result += '    <tr><td style="text-align:left" colspan=3>{desc}</td></tr>\n'.format(desc=group.description())
                section_entries += group_entries
                if group_entries > 0:
                    section_result += group_result
            if section_entries > 0:
                result += section_result
        result += r'</table>'

        return(result)


    @staticmethod
    def FromWCxf(wc):

        # the following coefficients are treated as real-valued in EOS
        real_coeffs = [
            'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c8', 'c8\''
        ]

        parameters = _Parameters.Defaults()
        for name in wc.dict:
            qn     = eos.QualifiedName(name)
            prefix = qn.prefix_part()
            coeff  = qn.name_part()
            value  = wc.dict[name]
            if (not value.imag == 0) and (coeff in real_coeffs):
                raise ValueError('WC {0} does not support non-zero imaginary part'.format(name))

            # Add values provided by WCxf to EOS central (SM) values
            if str(coeff) in real_coeffs:
                p = parameters[str(prefix) + '::' + str(coeff)]
                p.set(p.central() + value.real)
            else:
                pr = parameters[str(prefix) + '::Re{' + str(coeff) + '}']
                pi = parameters[str(prefix) + '::Im{' + str(coeff) + '}']
                pr.set(pr.central() + value.real)
                pi.set(pi.central() + value.imag)

        return parameters


