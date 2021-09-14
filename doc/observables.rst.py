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
print('The following is a list of observables that is available in EOS v{}.\n\n'.format(eos.version()))
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
        translation = {
            ord(':'): 'co', ord('@'): 'at', ord('/'): 'sl', ord('_'): 'un',
            ord('('): 'po', ord(')'): 'pc', ord('+'): 'pp', ord('-'): 'mm',
            ord('>'): 'to'
        }
        for qn, _ in group:
            print('.. _{qn}:'.format(qn=str(qn).translate(translation).lower()))
        print('.. list-table::')
        print('   :widths: 50, 25, 25')
        print('   :header-rows: 1')
        print('')
        print('   * - Qualified Name')
        print('     - Description')
        print('     - Kinematic Variables')
        for qn, entry in group:
            latex = entry.latex()
            unit = entry.unit()
            if unit == eos.Unit.Unity():
                unit_string = ''
            else:
                print(entry.name())
                unit_string = r'\, \left[ {unit_string} \right]'.format(unit_string=unit.latex())

            if 0 == len(latex):
                continue

            entry_name = str(qn)
            entry_desc = latex_to_rst(latex + unit_string)
            entry_kv   = ', '.join([('``' + kv + '``') for kv in entry.kinematic_variables()])
            print('   * - ``{qn}``'.format(qn=entry_name))
            print('     - :math:`{desc}`'.format(desc=entry_desc))
            print('     - {kv}'.format(kv=entry_kv))


        print('\n\n')
        group_desc  = latex_to_rst(group.description())
        print(group_desc)
        print('\n\n')

    print('\n\n')
print('\n\n')
