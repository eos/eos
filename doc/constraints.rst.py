#!/usr/bin/python3
# vim: set sw=4 sts=4 et tw=120 :

import eos
import re
import yaml

constraints = eos.Constraints()

def latex_to_rst(s):
    return(re.sub(r'\$([^\$]*)\$', r':math:`\1`', s))

page_title = 'List of Constraints'
print('#' * len(page_title))
print(page_title)
print('#' * len(page_title))
print('\n')
print('The following is the full list of likelihood constraints (both experimental and theoretical) included in EOS v{}.\n\n'.format(eos.__version__))
print('\n')
print('.. list-table::')
print('   :widths: 50,50')
print('   :header-rows: 1')
print('')
print('   * - Qualified Name')
print('     - Observables')
for qn, entry in constraints:
    print('   * - ``{qn}``'.format(qn=qn))
    data = yaml.load(entry.serialize(), Loader=yaml.SafeLoader)
    translation = {
        ord(':'): 'co', ord('@'): 'at', ord('/'): 'sl', ord('_'): 'un',
        ord('('): 'po', ord(')'): 'pc', ord('+'): 'pp', ord('-'): 'mm',
        ord('>'): 'to'
    }
    unique_observables = ''
    if 'observable' in data:
        unique_observables = set([str(data['observable'])])
    elif 'observables' in data:
        unique_observables = set([str(o) for o in data['observables']])

    entries = ['`{qn} <observables.html#{anchor}>`_'.format(qn=qn, anchor=qn.translate(translation).lower()) for qn in unique_observables]
    observables = ', '.join(entries)
    print('     - {}'.format(observables))
print('\n\n')
