import eos
import re
from jinja_util import print_template

def latex_to_rst(s):
    return(re.sub(r'\$([^\$]*)\$', r':math:`\1`', s))

qn_translation = {
    ord(':'): 'co', ord('@'): 'at', ord('/'): 'sl', ord('_'): 'un',
    ord('('): 'po', ord(')'): 'pc', ord('+'): 'pp', ord('-'): 'mm',
    ord('>'): 'to'
}

def make_doc_observables(group):
    "Make structured string data from eos group"
    observables = []
    for qn, entry in group:
        latex = entry.latex()

        unit_string = None
        if entry.unit() == eos.Unit.Unity():
            unit_string = ''
        else:
            unit_string = r'\, \left[ {unit_string} \right]'.format(unit_string=entry.unit().latex())

        if len(latex) > 0:
            description = latex_to_rst(latex + unit_string)
        else:
            description = r'' # signal latex string n/a

        kinematics = ', '.join([('``' + kv + '``') for kv in entry.kinematic_variables()])

        observables.append({
            'qn': str(qn),
            'description': description,
            'kinematics': kinematics,
            'link_key': str(qn).translate(qn_translation).lower(),
        })
    return observables


# mirror the eos object structure with strings
sections = []
for section in eos.Observables().sections():
    groups = []
    for group in section:
        groups.append({
            'name': latex_to_rst(group.name()),
            'observables': make_doc_observables(group),
            'description': latex_to_rst(group.description()),
        })

    sections.append({
        'name': latex_to_rst(section.name()),
        'groups': groups
    })


print_template(__file__,
    len = len,
    version = eos.version(),
    sections = sections
)

