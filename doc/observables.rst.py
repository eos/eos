import eos
import re
from jinja_util import print_template, qn_to_link_map

def latex_to_rst(s):
    return(re.sub(r'\$([^\$]*)\$', r':math:`\1`', s))

# Mirror the EOS observables hierarchy as structured string data: arrays of
# dicts. Observables are organised in groups, which are organised in sections.

def make_doc_sections():
    return [{
        'name'  : latex_to_rst(section.name()),
        'groups': make_doc_groups(section)
        } for section in eos.Observables().sections()]

def make_doc_groups(section):
    return [{
        'name'       : latex_to_rst(group.name()),
        'description': latex_to_rst(group.description()),
        'observables': make_doc_observables(group),
        } for group in section]

def make_doc_observables(group):
    observables = []
    for qn, entry in group:

        unit_string = None
        if entry.unit() == eos.Unit.Unity():
            unit_string = ''
        else:
            unit_string = fr'\, \left[ {entry.unit().latex()} \right]'

        latex = entry.latex()
        if len(latex) > 0:
            description = latex_to_rst(latex + unit_string)
        else:
            description = r'' # signal latex string n/a

        qualified_name = str(qn)
        kinematics     = ', '.join([('``' + kv + '``') for kv in entry.kinematic_variables()])

        observables.append({
            'qualified_name': qualified_name,
            'description'   : description,
            'kinematics'    : kinematics,
            'link_key'      : qualified_name.translate(qn_to_link_map).lower(),
        })

    return observables


if __name__ == '__main__':

    print_template(__file__,
        version = eos.__version__,
        sections = make_doc_sections(),
        len = len,
    )
