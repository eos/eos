import eos
import re
from jinja_util import print_template, qn_to_link_map

def latex_to_rst(s):
    s = re.sub(r'\$([^\$]*)\$', r':math:`\1`', s) # inline math
    s = re.sub(r'(\\begin{equation})([^\$]*?)(\\end{equation})', r'\n.. math::\n\2', s) # equation
    s = re.sub(r'(\\begin{align})([^\$]*?)(\\end{align})', r'\n.. math::\n\2', s) # align

    return(s)

# Mirror the EOS parameters hierarchy as structured string data: arrays of
# dicts. Parameters are organised in groups, which are organised in sections.

def make_doc_sections():
    return [{
        'name'  : latex_to_rst(section.name()),
        'groups': make_doc_groups(section)
        } for section in eos.Parameters().sections()]

def make_doc_groups(section):
    return [{
        'name'       : latex_to_rst(group.name()),
        'description': latex_to_rst(group.description()),
        'parameters' : make_doc_parameters(group),
        } for group in section]

def make_doc_parameters(group):
    parameters = []
    for param in group:
        qn = str(param.name())

        parameters.append({
            'qualified_name': qn,
            'latex'         : latex_to_rst(f'{param.latex()}'),
            'unit'          : latex_to_rst(f'${param.unit().latex()}$'),
            'value'         : param.evaluate(),
            'link_key'      : qn.translate(qn_to_link_map).lower(),
        })

    return parameters


if __name__ == '__main__':

    print_template(__file__,
        version = eos.__version__,
        sections = make_doc_sections(),
        len = len,
    )
