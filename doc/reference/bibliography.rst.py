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


def id_to_eprint_url(eprint_id):
    if eprint_id.startswith('oai:arXiv.org:'):
        return 'arXiv:' + eprint_id.lstrip('oai:arXiv.org:'), 'https://arxiv.org/abs/' + eprint_id.lstrip('oai:arXiv.org:')

    return None, None


def make_references():
    result = []
    for handle, reference in eos.References():
        title = latex_to_rst(reference.title())
        eprint, url = id_to_eprint_url(reference.eprint_id())
        data = {
            'authors': reference.authors(),
            'title':   title,
            'eprint':  eprint,
            'url':     url
        }
        result.append((handle, data))
    return result


if __name__ == '__main__':

    print_template(__file__,
        version = eos.__version__,
        references = make_references(),
        len = len,
    )
