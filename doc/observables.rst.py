#!/usr/bin/python3
# vim: set sw=4 sts=4 et tw=120 :

import eos
import re

observables = eos.Observables()

def latex_to_rst(s):
    return(re.sub(r'\$([^\$]*)\$', r':math:`\1`', s))

all_data = []
page_title = 'List of Observables'
print('#' * len(page_title))
print(page_title)
print('#' * len(page_title))
print('\n')
print('The following is a (incomplete) list observables that EOS can predict.\n\n')
for section in observables.sections():
    section_title = latex_to_rst(section.name())
    print('*' * len(section_title))
    print(section_title)
    print('*' * len(section_title))

    for group in section:
        group_title = latex_to_rst(group.name())
        print(group_title)
        print('=' * len(group_title))

        print('\n')
        print('.. list-table::')
        print('   :widths: 25, 25')
        print('   :header-rows: 1')
        print('')
        print('   * - Qualified Name')
        print('     - Description')
        for qn, entry in group:
            latex = entry.latex()
            if 0 == len(latex):
                continue

            entry_name = str(qn)
            entry_desc = latex_to_rst(latex)
            print('   * - ``{qn}``'.format(qn=entry_name))
            print('     - :math:`{desc}`'.format(desc=entry_desc))


        print('\n\n')
        group_desc  = latex_to_rst(group.description())
        print(group_desc)
        print('\n\n')

    print('\n\n')
print('\n\n')
