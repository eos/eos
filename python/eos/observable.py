# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2019 Danny van Dyk
# Copyright (c) 2021 Philip Lüghausen
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

from _eos import _Observables

class Observables(_Observables):
    """
    Represents the complete list of observables known to EOS.

    Objects of this class are visualized as tables in Jupyter notebooks for easy
    overview. Filters can be applied through keyword arguments to the initialization.

    :param prefix: Only show observables whose qualified names contain the provided ``prefix`` in their prefix part.
    :type prefix: str
    :param name: Only show observables whose qualified names contain the provided ``name`` in their name part.
    :type name: str
    :param suffix: Only show observables whose qualified names contain the provided ``suffix`` in their suffix part.
    :type suffix: str

    See also `the complete list of observables <../reference/observables.html>`_ in this documentation.
    """
    def __init__(self, prefix=None, name=None, suffix=None, showall=False):
        super().__init__()
        self.prefix=prefix
        self.name=name
        self.suffix=suffix
        self.showall=showall

    def filter_entry(self, qn):
        if self.prefix and not self.prefix in str(qn.prefix_part()):
            return False

        if self.name and not self.name in str(qn.name_part()):
            return False

        if self.suffix and not self.suffix in str(qn.suffix_part()):
            return False

        return True

    def _repr_html_(self):
        result = r'''
        <script>
            function toggle_group(group_title, id) {
                var table = group_title.parentNode.parentNode.parentNode.parentNode
                var query = 'tbody[id="' + id + '"]'
                var group = table.querySelector(query)
                if (group.style.visibility == "collapse") {
                    group.style.visibility = "visible"
                } else {
                    group.style.visibility = "collapse"
                }
            }
            function toggle_av(opt_anchor, id) {
                var query_dots   = 'span.dots[id="' + id + '"]'
                var query_values = 'span.values[id="' + id + '"]'
                var dots   = opt_anchor.querySelector(query_dots)
                var values = opt_anchor.querySelector(query_values)
                if (dots.style.display == "none") {
                    dots.style.display   = "inline"
                    values.style.display = "none"
                } else {
                    dots.style.display   = "none"
                    values.style.display = "inline"
                }
            }
        </script>
        <style>
            td.qn     { text-align: left;   }
            td.sym    { text-align: center; }
            td.unit   { text-align: right;  }
            td.optkey { text-align: left;   }
            td.optav  { text-align: left;   }
            td.optdef { text-align: left;   }
        </style>
        <table>
            <colgroup>
                <col width="25%" id="qn"          style="min-width: 200px; text-align: left">
                <col width="20%" id="symbol"      style="min-width: 200px">
                <col width="5%"  id="unit"        style="min-width:  50px">
                <col width="20%" id="kv"          style="min-width: 200px">
                <col width="10%" id="opt-key"     style="min-width:  75px">
                <col width="10%" id="opt-allowed" style="min-width:  75px">
                <col width="10%" id="opt-default" style="min-width:  75px">
            </colgroup>
            <thead>
                <tr>
                    <th rowspan="2">qualified name</th>
                    <th rowspan="2">symbol</th>
                    <th rowspan="2">unit</th>
                    <th rowspan="2">kinematic<br> variables</th>
                    <th colspan=3>options</th>
                </tr>
                <tr>
                    <th>key</th>
                    <th>values</th>
                    <th>default</th>
                </tr>
            </thead>
        '''
        group_id = 0
        observable_id = 0
        for section in self.sections():
            section_entries = 0
            section_result = fr'''
                <tr>
                    <th style="text-align:left" colspan=8><big>{section.name()}</big></th>
                </tr>'''
            for group in section:
                group_entries = 0
                group_result  = fr'''
                    <tbody>
                        <tr>
                            <th style="text-align:left" colspan=8>
                                <a style="text-decoration: none" onclick="toggle_group(this, 'grp{group_id}')">{group.name()}</a>
                            </th>
                        </tr>
                    </tbody>
                '''#.format(group=group.name(), id=group_id)
                group_result += fr'''
                    <tbody style="visibility:collapse" id="grp{group_id}">
                    <tr>
                        <td style="text-align:left" colspan=8>{group.description()}</td>
                    </tr>
                '''
                for qn, entry in group:
                    if not self.filter_entry(qn):
                        continue

                    latex = entry.latex()
                    if (0 == len(latex)) and not self.showall:
                        continue

                    unit = entry.unit().latex()
                    if unit == '1':
                        unit = '&mdash;'
                    else:
                        unit = fr'$$\left[ {unit} \right]$$'
                    if not self.filter_entry(qn):
                        continue
                    if (0 == len(latex)) and not self.showall:
                        continue

                    kinematic_variables = r'<br>'.join(['<tt>' + str(kv) + '</tt>' for kv in entry.kinematic_variables()])
                    if len(kinematic_variables) == 0:
                        kinematic_variables = '&mdash;'

                    keys           = []
                    allowed_values = []
                    default_values = []
                    for option_idx, option in enumerate(entry.options()):
                        id  = fr'grp{group_id}-obs{observable_id}-opt{option_idx}'
                        av  = fr'''<a onclick="toggle_av(this, '{id}')">
                            <span class="dots"   id="{id}" style="display: inline; text-align: left">...</span>
                            <span class="values" id="{id}" style="display: none;   text-align: left">
                        '''
                        av += '   ' + r'<br/>'.join([fr'''<tt>{av}</tt>''' for av in option.allowed_values])
                        av += r'''
                            </span>
                        </a>'''

                        keys           += [fr'''<tt>{option.key}</tt>''']
                        default_values += [fr'''<tt>{option.default_value}</tt>''']
                        allowed_values += [av]
                    if len(keys) == 0:
                        keys           += ['&mdash;']
                        allowed_values += ['&mdash;']
                        default_values += ['&mdash;']

                    rows = len(keys)
                    #print(f'qn = {qn}')
                    #print(f'len(keys) = {len(keys)}')
                    #print(f'len(allowed_values) = {len(allowed_values)}')
                    #print(f'len(default_values) = {len(default_values)}')
                    #print(f'len(spans) = {len(spans)}')

                    group_result += fr'''
                        <tr>
                            <th class="qn"     rowspan="{rows}"><tt>{qn}</tt></th>
                            <td class="sym"    rowspan="{rows}">$${latex}$$</td>
                            <td class="unit"   rowspan="{rows}">{unit}</td>
                            <td class="kv"     rowspan="{rows}">{kinematic_variables}</td>
                            <td class="optkey" rowspan="1">{keys[0]}</td>
                            <td class="optav"  rowspan="1">{allowed_values[0]}</td>
                            <td class="optdef" rowspan="1">{default_values[0]}</td>
                        </tr>
                    '''
                    options = list(zip(keys, allowed_values, default_values))
                    for key, allowed_value, default_value in options[1:]:
                        group_result += fr'''
                            <tr>
                                <td class="optkey" rowspan="1">{key}</td>
                                <td class="optav"  rowspan="1">{allowed_value}</td>
                                <td class="optdef" rowspan="1">{default_value}</td>
                            </tr>
                        '''

                    group_entries += 1
                    observable_id += 1

                group_result += '    </tbody>'
                group_id += 1
                section_entries += group_entries
                if group_entries > 0:
                    section_result += group_result
            if section_entries > 0:
                result += section_result
        result += r'</table>'

        return(result)
