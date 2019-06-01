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

from _eos import _Constraints

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
        result = '<table>\n'
        result += '<tr><th style="text-align:left">Name</th><th style="text-align:left">Type</th></tr>'
        for qn, entry in self:
            if not self.filter_entry(qn):
                continue
            result += '<tr><td><tt style="color:grey">{qn}</tt></td><td style="text-align:left">{t}</td></tr>'.format(qn=qn,t=entry.type())
        result += '</table>'

        return(result)
