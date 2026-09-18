#!/usr/bin/python3
# vim: set sw=4 sts=4 et tw=120 :

import eos
import re
from jinja_util import print_template

replacements = [
    (re.compile(r'(\\GeV)'),      r'\\text{GeV}'),
    # docutils rejects inline markup that is padded with whitespace, so strip the delimited text
    (re.compile(r'\$\s*([^\$]*?)\s*\$'), r':math:`\1`'),
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
        id = eprint_id.removeprefix('oai:arXiv.org:')
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


def ref_to_inspire(reference):
    inspire_id = reference.inspire_id()
    if not inspire_id:
        return None

    return {
            'id':    inspire_id,
            'url':   f'https://inspirehep.net/literature?q=texkey:{inspire_id}',
            'badge': f'https://img.shields.io/badge/INSPIRE-{ inspire_id.replace("-", "--") }-yellow.svg',
            'alt':   f'INSPIRE:{inspire_id}'
        }


def ref_to_url(reference, handle):
    url = reference.url()
    if not url:
        return None

    return {
            'id':    handle,
            'url':   url,
            'badge': f'https://img.shields.io/badge/LINK-blue.svg',
            'alt':   f'URL:{url}'
        }


def make_references():
    result = []
    for handle, reference in eos.References():
        title = latex_to_rst(reference.title())
        eprint = ref_to_eprint(reference)
        inspire = ref_to_inspire(reference)
        url = ref_to_url(reference, handle)
        data = {
            'authors': reference.authors(),
            'title':   title,
            'eprint':  eprint,
            'inspire': inspire,
            'url': url,
        }
        result.append((handle, data))
    return result


def make_links(references):
    # several references may cite the same link, while each substitution may only be defined once
    result = {}
    for _, reference in references:
        for link in (reference['eprint'], reference['inspire'], reference['url']):
            if link:
                result.setdefault(link['id'], link)
    return list(result.values())


if __name__ == '__main__':

    references = make_references()

    print_template(__file__,
        version = eos.__version__,
        references = references,
        links = make_links(references),
        len = len,
    )
