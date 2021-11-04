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

from _eos import _References

class References(_References):
    def __init__(self, year=None, index=None):
        super().__init__()
        self.year  = year
        self.index = index

    def filter_entry(self, rn):
        if self.year and not self.year in str(rn.year_part()):
            return False

        if self.index and not self.index in str(rn.index_part()):
            return False

        return True

    def _repr_html_(self):
        result = r'''<table>
            <colgroup>
                <col width="10%" id="rn"      style="min-width: 50px">
                <col width="40%" id="title"   style="min-width: 200px">
                <col width="15%" id="authors" style="min-width: 100px">
            </colgroup>
            <thead>
                <tr>
                    <th>name</th>
                    <th>title</th>
                    <th>authors</th>
                </tr>
            </thead>'''
        for rn, entry in self:
            if not self.filter_entry(rn):
                continue

            title          = entry.title()
            authors        = '<br/>'.join([fr'{a.strip()}' for a in entry.authors().split(' and ')])
            eprint_archive = entry.eprint_archive()
            eprint_id      = entry.eprint_id()

            link = ''
            if 'arXiv' == eprint_archive:
                link = fr' href="https://arxiv.org/abs/{eprint_id.split(":")[-1]}"'

            result += fr'''
                <tr>
                    <td><a {link}><tt>{rn}</tt></a></td>
                    <td>{title}</td>
                    <td>{authors}</td>
                </tr>
            '''

        result += r'</table>'

        return(result)
