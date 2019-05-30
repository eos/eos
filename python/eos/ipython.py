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

import scipy

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

def __format_ObservableEntry(e):
    result = '<table>\n'
    result += '<tr><th>QualifedName</th><td><tt style="color:grey">{qn}</tt></td></tr>'.format(qn=e.name())
    result += '<tr><th>Description</th><td>$${latex}$$</td></tr>'.format(latex=e.latex())
    kvs = [kv for kv in e.kinematic_variables()]
    if len(kvs) > 0:
        result += '<tr><th rowspan={rows}>Kinematic Variables</th><td>{kv}</td></tr>'.format(rows=len(kvs),kv=kvs[0])
    for i in range(1, len(kvs)):
        result += '<tr><td>{kv}</td></tr>'.format(kv=kvs[i])
    result += '</table>'
    return(result)

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

def __format_GoodnessOfFit(gof):
    result = '<table>\n'
    result += '<tr><th>constraint</th><th>&chi;<sup>2</sup></th><th>d.o.f.</th></tr>\n'
    for entry in gof:
        result += '<tr><td><tt>{name}</tt></td><td>{chi2:6.4f}</td><td>{dof}</td></tr>\n'.format(
            name=entry[0], chi2=entry[1].chi2, dof=entry[1].dof)
    result += '</table><br/>\n'
    chi2 = gof.total_chi_square()
    dof  = gof.total_degrees_of_freedom()
    pvalue = 1.0 - scipy.stats.chi2(dof).cdf(chi2)
    result += '<table>\n'
    result += '<tr><th>total &chi;<sup>2</sup></th><td>{chi2:6.4f}</td></tr>\n'.format(chi2=chi2)
    result += '<tr><th>total degrees of freedom</th><td>{dof}</td></tr>\n'.format(dof=dof)
    result += '<tr><th>p-value</th><td>{p:6.4f}%</td></tr>\n'.format(p=pvalue * 100)
    result += '</table>\n'
    return(result)
