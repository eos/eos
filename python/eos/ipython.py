# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2018 Danny van Dyk
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

def __format_Parameter(p):
    name = ''
    latex = p.latex()
    if latex:
        name = r'$' + latex + r'$'
    else:
        name = p.name()

    return(
    """
    <table>
        <tr>
            <th>{name}</th>
            <td>(eos.Parameter)</td>
        </tr>
        <tr>
            <th>current value</th>
            <td><tt>{value}</tt></td>
        </tr>
        <tr>
            <th>default value</th>
            <td><tt>{central}</tt></td>
        </tr>
    </table>""".format(
        name=name,
        value=p.evaluate(),
        central=p.central()
    ))

def __format_KinematicVariable(kv):
    name = kv.name()

    return("""
        <table>
            <tr>
                <th>{name}</th>
                <td>(eos.KinematicVariable)</td>
            </tr>
            <tr>
                <th>current value</th>
                <td><tt>{value}</tt></td>
            </tr>
        </table>""".format(
            name=name,
            value=kv.evaluate()
    ))

def __format_Observable(obs):
    name = obs.name()

    return("""
        <table>
            <tr>
                <th>{name}</th>
                <td>(eos.Observable)</td>
            </tr>
            <tr>
                <th>current value</th>
                <td><tt>{value:1.4g}</tt></td>
            </tr>
        </table>""".format(
            name=name,
            value=obs.evaluate()
    ))
