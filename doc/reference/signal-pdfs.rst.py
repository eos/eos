import eos
import re
from jinja_util import print_template, qn_to_link_map

def latex_to_rst(s):
    return(re.sub(r'\$([^\$]*)\$', r':math:`\1`', s))

# Mirror the EOS signal PDFs hierarchy as structured string data: arrays of
# dicts. Observables are organised in groups, which are organised in sections.

def make_doc_sections():
    return [{
        'name'  : latex_to_rst(section.name()),
        'groups': make_doc_groups(section)
        } for section in eos.SignalPDFs().sections()]

def make_doc_groups(section):
    return [{
        'name'       : latex_to_rst(group.name()),
        'description': latex_to_rst(group.description()),
        'signal_pdfs': make_doc_signal_pdfs(group),
        } for group in section]

def make_doc_signal_pdfs(group):
    signal_pdfs = []
    for qn, _ in group:
        qualified_name = str(qn)

        signal_pdfs.append({
            'qualified_name': qualified_name,
            'link_key': qualified_name.translate(qn_to_link_map).lower(),
        })

    return signal_pdfs


if __name__ == '__main__':

    print_template(__file__,
        version = eos.__version__,
        sections = make_doc_sections(),
        len = len,
    )
