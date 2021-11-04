# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2019 Danny van Dyk
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

from _eos import _Constraints, _References

class Constraints(_Constraints):
    def __init__(self, prefix=None, name=None, suffix=None):
        super().__init__()
        self.prefix=prefix
        self.name=name
        self.suffix=suffix

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
            function toggle_obs(obs_anchor, id) {
                var query_dots   = 'span.dots[id="' + id + '"]'
                var query_values = 'span.values[id="' + id + '"]'
                var dots   = obs_anchor.querySelector(query_dots)
                var values = obs_anchor.querySelector(query_values)
                if (dots.style.display == "none") {
                    dots.style.display   = "inline"
                    values.style.display = "none"
                } else {
                    dots.style.display   = "none"
                    values.style.display = "inline"
                }
            }
        </script>
        <table>
            <colgroup>
                <col width="50%" id="qn"     style="min-width: 200px">
                <col width="25%" id="type"   style="min-width: 200px">
                <col width="15%" id="type"   style="min-width: 100px">
                <col width="10%" id="ref"    style="min-width: 100px">
            </colgroup>
            <thead>
                <tr>
                    <th>qualified name</th>
                    <th>observables</th>
                    <th>type</th>
                    <th>reference</th>
                </tr>
            </thead>'''

        constraint_id = 0
        references = _References()
        for qn, entry in self:
            if not self.filter_entry(qn):
                continue

            id           = fr'con{constraint_id}-obs'
            observables  = fr'''<a onclick="toggle_obs(this, '{id}')">
                <span class="dots"   id="{id}" style="display: inline; text-align: left">...</span>
                <span class="values" id="{id}" style="display: none;   text-align: left">
            '''
            observables += '   ' + r'<br/>'.join({fr'''<tt>{o}</tt>''' for o in entry.observables()})
            observables += r'''
                </span>
            </a>'''
            refname     = str(qn.suffix_part())
            reflink     = ''
            try:
                ref         = references[refname]
                if 'arXiv' == ref.eprint_archive():
                    reflink = fr' href="https://arxiv.org/abs/{ref.eprint_id().split(":")[-1]}"'
            except:
                pass

            result += fr'''
                <tr>
                    <td><tt>{qn}</tt></td>
                    <td>{observables}</td>
                    <td>{entry.type()}</td>
                    <td><a "{reflink}">{refname}</a></td>
                </tr>'''
        result += '''
            </table>
        '''

        return(result)
