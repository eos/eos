#!/usr/bin/python3
# vim: set sw=4 sts=4 et tw=120 :

import eos
import re

constraints = eos.Constraints()

def latex_to_rst(s):
    return(re.sub(r'\$([^\$]*)\$', r':math:`\1`', s))

page_title = 'List of Constraints'
print('#' * len(page_title))
print(page_title)
print('#' * len(page_title))
print('\n')
print('The following is the full list of constraints (both experimental and theoretical as of EOS v{}.\n\n'.format(eos.version()))
print('\n')
print('.. list-table::')
print('   :widths: 25')
print('   :header-rows: 1')
print('')
print('   * - Qualified Name')
for qn, entry in constraints:
    print('   * - ``{qn}``'.format(qn=qn))
print('\n\n')
