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

tex_translation = {
    'charmonium': r'\Psi', 'gamma': r'\gamma'
}

eos_sections = eos.Observables().sections()

def make_doc_observables(group):
    "Make structured string data from eos group"
    doc_obss = []
    for qn, entry in group:

        latex = entry.latex()

        unit_string = None
        if entry.unit() == eos.Unit.Unity():
            unit_string = ''
        else:
            unit_string = r'\, \left[ {unit_string} \right]'.format(unit_string=entry.unit().latex())

        if len(latex) > 0:
            entry_desc = latex_to_rst(latex + unit_string)
        else:
            entry_desc = r'' # signal latex string n/a

        entry_name = str(qn)
        entry_kv   = ', '.join([('``' + kv + '``') for kv in entry.kinematic_variables()])

        doc_obss.append({
            'qn': entry_name,
            'description': entry_desc,
            'kinematics': entry_kv,
            'link_key': entry_name.translate(qn_translation).lower(),
        })
    links = []
    return doc_obss


# mirror the eos object structure with strings
doc_sections = []
for section in eos_sections:

    doc_groups = []
    for group in section:

        doc_obss = make_doc_observables(group)
        group_title = group.name()
        # TODO: this substitution is not safe if the group title contains
        # the substituted names as TeX commands. Example: we substitute gamma
        # to \gamma, which gives \\gamma if the command is already used somewhere.
        for k, v in tex_translation.items():
            group_title = group_title.replace(k, v)
        doc_group = {
            'title': latex_to_rst(group_title),
            'obss': doc_obss,
            'description': latex_to_rst(group.description()),
        }
        doc_groups.append(doc_group)

    doc_section = {
        'title': latex_to_rst(section.name()),
        'groups': doc_groups
    }
    doc_sections.append(doc_section)


print_template(__file__,
    len = len,
    eos_version = eos.version(),
    sections = doc_sections
)

