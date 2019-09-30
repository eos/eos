#!/usr/bin/env python3
# vim: set sw=4 sts=4 et tw=120 :

import eos

observables = eos.Observables()

for section in observables.sections():

    section_title = section.name()
    print('\section{' + section_title + '}\n')

    for group in section:
        group_title = group.name()
        print('\subsection{' + group_title + '}\n')

        print(r'{%')
        print(r'\captionsetup{type=table}')
        print(r'\renewcommand{\arraystretch}{1.1}')
        print(r'\newcolumntype{L}[1]{>{\raggedright\arraybackslash}p{#1}}')
        print(r'\begin{center}')
        print(r'\begin{longtable}{@{} p{.5\textwidth} p{.5\textwidth} @{}}')
        print(r'\toprule')
        print(r'    \textbf{Qualified Name} & \textbf{Description} \\')
        print(r'\midrule')
        for qn, entry in group:
            latex = entry.latex()
            if 0 == len(latex):
                continue

            print(r'    \verb|{0}| & {1} \\'.format(str(qn), r'$' + latex + r'$'))


        print(r'\bottomrule')
        print(r'\end{longtable}')
        print(r'\end{center}')
        print(r'\captionof{{table}}{{ {0} }}'.format(group.description()))
        print(r'}')
