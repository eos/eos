# Copyright (c) 2017-2021 Danny van Dyk
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

from _eos import _Parameters, QualifiedName, Unit
import re

class Parameters(_Parameters):
    """
    Represents the set of parameters known to EOS.

    Objects of this class are visualized as tables in Jupyter notebooks for easy
    overview. Filters can be applied through keyword arguments to the initialization.
    An independent instance of the default parameters is obtained by each call to the constructor.

    :param prefix: Only show parameters whose qualified names contain the provided ``prefix`` in their prefix part.
    :type prefix: str
    :param name: Only show parameters whose qualified names contain the provided ``name`` in their name part.
    :type name: str
    :param suffix: Only show parameters whose qualified names contain the provided ``suffix`` in their suffix part.
    :type suffix: str

    See also `the complete list of parameters <../reference/parameters.html>`_.
    """
    def __new__(cls, *args, **kwargs):
        instance = _Parameters.Defaults()
        instance.__class__ = Parameters
        return(instance)

    def __init__(self, prefix=None, name=None, suffix=None):
        self.prefix=prefix
        self.name=name
        self.suffix=suffix
        self._replacements = [
            (re.compile(r'(\\GeV)'),      r'\\text{GeV}'),
            (re.compile(r'\$([^\$]*)\$'), r'$$\1$$'),        # ensure display mode MathJax presentation
        ]

    def filter_entry(self, qn):
        qn = QualifiedName(str(qn))
        if self.prefix and not self.prefix in str(qn.prefix_part()):
            return False

        if self.name and not self.name in str(qn.name_part()):
            return False

        if self.suffix and not self.suffix in str(qn.suffix_part()):
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
        result = r'''
        <table>
            <colgroup>
                <col width="25%" id="qn"      style="min-width: 200px">
                <col width="20%" id="symbol"  style="min-width: 200px">
                <col width="5%"  id="unit"    style="min-width:  50px">
                <col width="10%" id="value"   style="min-width: 100px">
            </colgroup>
            <thead>
                <tr>
                    <th>qualified name</th>
                    <th>symbol</th>
                    <th>unit</th>
                    <th>value</th>
                </tr>
            </thead>
        '''
        for section in self.sections():
            section_entries = 0
            section_result = fr'''
              <tr><th style="text-align:left" colspan=4><big>{section.name()}</big></th></tr>'''
            for group in section:
                group_entries = 0
                group_result = fr'''
                    <tr><th style="text-align:left" colspan=4>{group.name()}</th></tr>
                    <tr><td style="text-align:left" colspan=4>{group.description()}</td></tr>'''
                group_parameters = [p for p in group]
                group_parameters.sort(key = Parameters._key)
                for parameter in group_parameters:
                    qn = parameter.name()
                    if not self.filter_entry(qn):
                        continue

                    latex = self._latex_refine(parameter.latex())
                    if 0 == len(latex):
                        latex = '---'

                    unit = parameter.unit()
                    if unit == Unit.Undefined() or unit.latex() == '1':
                        unit = '&mdash;'
                    else:
                        unit = fr'$$\left[ {unit.latex()} \right]$$'

                    value = parameter.evaluate()

                    group_result += fr'''
                        <tbody>
                            <tr>
                                <th><tt style="color:grey">{qn}</tt></th>
                                <td style="text-align:left">{latex}</td>
                                <td style="text-align:left">{unit}</td>
                                <td style="text-align:right">{value}</td>
                            </tr>
                        </tbody>'''
                    group_entries += 1

                section_entries += group_entries
                if group_entries > 0:
                    section_result += group_result
            if section_entries > 0:
                result += section_result
        result += r'''
            </table>
        '''

        return(result)


    def to_yaml(self, names=None):
        """
        Converts an eos.Parameters object to a YAML representation.

        :param names: Names of the parameters that shall be converted.
        :type name: iterable of str (optional)
        """

        import yaml

        parameters = []
        if names is None:
            parameters = [p for p in self]
        else:
            for p in self:
                if p.name() not in names:
                    continue

                parameters.append(p)

        contents = {}
        for p in parameters:
            contents.update({
                p.name(): {
                    'central': p.evaluate(),
                    'min':     p.min(),
                    'max':     p.max(),
                    'latex':   p.latex()
                }
            })

        return yaml.dump(contents)


    def dump(self, file, **kwargs):
        """
        Dumps an eos.Parameters object to a YAML file.

        :param file: Name of the file to which the parameters shall be written.
        :type file: str
        :param names: Names of the parameters that shall be converted.
        :type names: iterable of str (optional)
        """
        with open(file, 'w') as f:
            f.write(self.to_yaml(**kwargs))


    @staticmethod
    def FromWCxf(w):
        """
        Imports Wilson coefficients into an eos.Parameters object from a wilson.Wilson object.

        :param w: WET Wilson coefficients in the EOS basis.
        :type w: wilson.Wilson
        """

        # extract wc data
        wc = w.wc

        if wc.basis != 'EOS':
            raise ValueError('Wilson coefficients must be passed to this function in the EOS basis')

        # Needs to be revisited, once EOS support WET-4 or WET-3 bases
        if wc.scale != 4.2:
            raise ValueError('Wilson coefficients must be passed to this function at the reference scale 4.2 GeV')

        # the following coefficients are treated as real-valued in EOS
        real_coeffs = [
            'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c8', 'c8\''
        ]

        parameters = _Parameters.Defaults()
        for name, value in wc.dict.items():
            qn     = QualifiedName(name)
            prefix = qn.prefix_part()
            coeff  = qn.name_part()

            # Add values provided by WCxf to EOS central (SM) values
            if (prefix == 'b->s') and (coeff in real_coeffs):
                if (abs(value.imag) > 1.0e-10):
                    raise ValueError(f'imaginary part of WC {name} larger than 10^-10 threshold, which is not supported at present')

                p = parameters[str(prefix) + '::' + str(coeff)]
                p.set(p.central() + value.real)
            else:
                pr = parameters[str(prefix) + '::Re{' + str(coeff) + '}']
                pi = parameters[str(prefix) + '::Im{' + str(coeff) + '}']
                pr.set(pr.central() + value.real)
                pi.set(pi.central() + value.imag)

        return parameters
