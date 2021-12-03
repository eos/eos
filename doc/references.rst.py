#!/usr/bin/python3
# vim: set sw=4 sts=4 et tw=120 :

import eos
import functools
import re

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


page_title = 'List of References'
print('#' * len(page_title))
print(page_title)
print('#' * len(page_title))
print('\n')
print('''
The following is the full list of references that were used to implement theory codes and provide likelihood constraints
or parameter values within EOS v{}.\n\n'''.format(eos.__version__))
for handle, reference in eos.References():
    title = latex_to_rst(reference.title())
    eprint, url = id_to_eprint_url(reference.eprint_id())
    print('.. [{handle}] {authors} {delim}'.format(handle=handle, authors=reference.authors(), delim='--' if title else ''))
    print('   {title} {delim}'.format(title=title, delim='--' if eprint else ''))
    if eprint:
        print('   `{eprint} <{url}>`_'.format(eprint=eprint,url=url))
    print()
