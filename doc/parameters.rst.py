#!/usr/bin/python3
# vim: set sw=4 sts=4 et tw=120 :

import eos
import re

parameters = eos.Parameters.Defaults()

def latex_to_rst(s):
    return(re.sub(r'\$([^\$]*)\$', r':math:`\1`', s))

page_title = 'List of Parameters'
print('#' * len(page_title))
print(page_title)
print('#' * len(page_title))
print('\n')
print('The following the full list of parameters and their values as of EOS v{}.\n\n'.format(eos.version()))
print('\n')
print('.. list-table::')
print('   :widths: 25, 25, 25')
print('   :header-rows: 1')
print('')
print('   * - Qualified Name')
print('     - Representation')
print('     - Default Value')
for parameter in parameters:
    qn = parameter.name()
    latex = parameter.latex()
    value = parameter.evaluate()
    print('   * - ``{qn}``'.format(qn=qn))
    if len(latex) > 0:
        print('     - :math:`{latex}`'.format(latex=latex))
    else:
        print('     - ---')
    print('     - {value}'.format(value=value))
print('\n\n')
