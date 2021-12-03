#!/usr/bin/python3
# vim: set sw=4 sts=4 et tw=120 :

import eos
import functools
import re

parameters = eos.Parameters.Defaults()

replacements = [
    (re.compile(r'(\\GeV)'),      r'\\text{GeV}'),
    (re.compile(r'\$([^\$]*)\$'), r':math:`\1`'),
]

def latex_to_rst(s):
    result = s
    for regexp, repl in replacements:
        result = regexp.sub(repl, result)

    return(result)

def key(x):
    if type(x[0]) is eos.QualifiedName:
        return (str(x[0].prefix_part()), str(x[0].suffix_part()), str(x[0].name_part()))
    elif type(x[0]) is str:
        return ('', '', x[0])
    else:
        raise ValueError('x[0] is of unexpected type {}'.format(str(type(x[0]))))

page_title = 'List of Parameters'
print('#' * len(page_title))
print(page_title)
print('#' * len(page_title))
print('\n')
print('The following is the full list of parameters and their values used in EOS v{}.\n\n'.format(eos.__version__))
for section in parameters.sections():
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
        print('   :widths: 25, 25, 25')
        print('   :header-rows: 1')
        print('')
        print('   * - Qualified Name')
        print('     - Representation')
        print('     - Default Value')

        group_parameters = []
        for p in group:
            name = None
            try:
                name = eos.QualifiedName(p.name())
            except RuntimeError:
                name = p.name()

            latex = latex_to_rst(p.latex())
            value = p.evaluate()

            group_parameters.append((name, latex, value))

        group_parameters.sort(key=key)

        for qn, latex, value in group_parameters:
            print('   * - ``{qn}``'.format(qn=qn))
            if len(latex) > 0:
                print('     - {latex}'.format(latex=latex))
            else:
                print('     - ---')
            print('     - {value}'.format(value=value))
        print('\n\n')
