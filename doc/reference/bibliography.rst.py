#!/usr/bin/python3
# vim: set sw=4 sts=4 et tw=120 :

import eos
import re
from jinja_util import print_template

replacements = [
    (re.compile(r'(\\GeV)'),      r'\\text{GeV}'),
    (re.compile(r'\$([^\$]*)\$'), r':math:`\1`'),
]

def latex_to_rst(s):
    result = s
    for regexp, repl in replacements:
        result = regexp.sub(repl, result)

    return(result)


def ref_to_eprint(reference):
    eprint_id = reference.eprint_id()
    eprint_archive = reference.eprint_archive()
    if eprint_id.startswith('oai:arXiv.org:'):
        id = eprint_id.lstrip('oai:arXiv.org:')
        result = {
                'id': id,
                'url': f'https://arxiv.org/abs/{id}',
                'badge': f'https://img.shields.io/badge/arXiv-{ id.replace("-", "--") }-red.svg',
                'alt': f'arXiv:{id}'
            }
        return result
    elif eprint_archive == 'CDS':
        id = eprint_id
        result = {
                'id': id,
                'url': f'https://cds.cern.ch/record/{id}?ln=en',
                'badge': f'https://img.shields.io/badge/CDS-{ id.replace("-", "--") }-blue.svg',
                'alt': f'CDS:{id}'
            }
        return result

    return None


def make_references():
    result = []
    for handle, reference in eos.References():
        title = latex_to_rst(reference.title())
        eprint = ref_to_eprint(reference)
        data = {
            'authors': reference.authors(),
            'title':   title,
            'eprint':  eprint
        }
        result.append((handle, data))
    return result


if __name__ == '__main__':

    print_template(__file__,
        version = eos.__version__,
        references = make_references(),
        len = len,
    )
