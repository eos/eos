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
        result  = r'<table>'
        result += r'<tr><th style="text-align:left">name</th><th>details</th></tr>'
        for rn, entry in self:
            if not self.filter_entry(rn):
                continue

            authors        = entry.authors()
            title          = entry.title()
            eprint_archive = entry.eprint_archive()
            eprint_id      = entry.eprint_id()

            details  = r'<table>'
            details += r'<tr><th>authors</th><td>{}</td></tr>'.format(authors)
            details += r'<tr><th>title</th><td>{}</td></tr>'.format(title)
            if eprint_archive == 'arXiv' and eprint_id:
                pos = eprint_id[::-1].find(':')
                details += '<tr><th>eprint</th><td><a href="https://arxiv.org/abs/{0}">{1}</a></td></tr>'.format(eprint_id[pos:], eprint_id)
            details += r'</table>'

            result += r'<tr><td>[{}]</td><td>{}</td></tr>'.format(rn, details)

        result += r'</table>'

        return(result)
