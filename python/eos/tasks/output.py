# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2021-2026 Danny van Dyk
# Copyright (c) 2024      Carolina Bolognani
#
# This file is part of the EOS project. EOS is free software;
# you can redistribute it and/or modify it under the terms of the GNU General
# Public License version 2, as published by the Free Software Foundation.
#
# EOS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 59 Temple
# Place, Suite 330, Boston, MA  02111-1307  USA

import eos
import eos.analysis_file_context
import numpy as _np
import os
import warnings

from dataclasses import asdict
from ._framework import task, _backed_up, _staged_file


def _get_references(analysis_file):
    all_references = eos.References()
    all_constraints = eos.Constraints()
    result = []
    for likelihood in analysis_file.likelihoods.values():
        constraints         = likelihood.constraints
        for constraint in constraints:
            constraint_name = constraint.constraint.split(';')[0]
            for reference in all_constraints[constraint_name].references():

                try:
                    result.append(all_references[reference].inspire_id())
                except KeyError:
                    warnings.warn(f'No reference found for {reference}')

    result = _np.unique(result)
    return result

# Create a report
@task('report', 'reports', staged=False)
def report(analysis_file:str, template_file:str, base_directory:str='./', generate_pdf:bool=True):
    """
    Generates a report from an analysis file and a Jinja2 template file.

    The processed file will be stored in the same directory as the template file, with the the extension '.jinja' or '.jinja2' stripped from the file name.
    If `generate_pdf` is set to True, the intermediate file will be converted to PDF using Pandoc and XeLaTeX. This requires a working installation of both programs.

    :param analysis_file: The name of the analysis file that shall be used for the report.
    :type analysis_file: str or `eos.AnalysisFile`
    :param template_file: The name of the Jinja2 template file that shall be used for the report.
    :type template_file: str
    :param base_directory: The base directory for the storage of data files. Can also be set via the EOS_BASE_DIRECTORY environment variable.
    :type base_directory: str, optional
    :param generate_pdf: The flag that enables the conversion of the intermediate file to PDF. Defaults to `True`.
    :type generate_pdf: bool, optional
    """
    import jinja2

    eos.inprogress(f'Generating report from template \'{template_file}\' for analysis file \'{analysis_file}\'')

    resource_directory = os.path.dirname(template_file)

    intermediate_file, template_ext = os.path.splitext(template_file)
    if template_ext not in ['.jinja', '.jinja2']:
        raise ValueError(f'Template file is expected to have an extension of \'.jinja\' or \'.jinja2\', but got \'{template_ext}\'')

    eos.info(f'Processing template \'{template_file}\' to intermediate \'{intermediate_file}\'')

    basename, intermediate_ext = os.path.splitext(intermediate_file)
    if intermediate_ext not in ['.md']:
        raise ValueError(f'Intermediate file is expected to have an extension of \'.md\', but got \'{intermediate_ext}\'')

    from eos.reporting import AnalysisData

    # Lazy, disk-driven access to all recorded data below EOS_BASE_DIRECTORY/data/. Nothing is loaded
    # until the template actually references a posterior's samples, modes, predictions, etc.
    analysis = AnalysisData(base_directory=base_directory, analysis_file=analysis_file)

    # Corner figures are drawn into a 'figures/' subdirectory next to the report and embedded through a
    # path relative to the report, so that both pandoc and a reader of the intermediate file resolve it.
    figures_directory = os.path.join(resource_directory, 'figures') if resource_directory else 'figures'

    def corner_figure(posterior, **kwargs):
        _backed_up(os.path.join(figures_directory, f'corner-{posterior.name}.{kwargs.get("format", "pdf")}'))
        filename = posterior.corner_figure(output_directory=figures_directory, **kwargs)
        return os.path.join('figures', filename)

    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(searchpath='./'),
        undefined=jinja2.StrictUndefined, # error for undefined variables in
                                          # template
        trim_blocks=True,
        lstrip_blocks=True)
    template = env.get_template(template_file)
    result = template.render(
        eos={ 'version': eos.__version__ },
        analysis_file=analysis_file,
        metadata=asdict(analysis_file._metadata),
        analysis=analysis,
        base_directory=base_directory,
        corner_figure=corner_figure,
        references=_get_references(analysis_file),
    )

    staged_intermediate_file = _staged_file(intermediate_file)
    with open(staged_intermediate_file, 'w') as f:
        f.write(result)

    if not generate_pdf:
        eos.completed('... report generation finished')
        return

    output_file = basename + '.pdf'

    eos.info(f'Processing intermediate file \'{intermediate_file}\' to final output \'{output_file}\'')

    import pypandoc

    pypandoc.convert_file(
        source_file=staged_intermediate_file,
        outputfile=_staged_file(output_file),
        to='pdf',
        extra_args=[
            f'--resource-path=./:{resource_directory}/:{eos._pkg_data_dir}/',
            f'--variable=eoslogo:{eos._pkg_data_dir}/report-logo.pdf',
            f'--variable=eosversion:{eos.__version__}',
            '--template=report-template.tex',
            '--pdf-engine=xelatex'
        ]
    )

    eos.completed('... report generation finished')


# Draw figures
@task('draw-figure', 'figures', mode=lambda **kwargs: 'a', staged=False)
def draw_figure(analysis_file:str, figure_name:str, base_directory:str='./', format:str|list[str]='pdf'):
    """
    Draws figures from the analysis file.

    The output files will be stored in EOS_BASE_DIRECTORY/figures/. If more than one format is
    requested, the figure is drawn once and stored once per format, since drawing the same figure
    repeatedly would stack its elements on top of each other.

    :param analysis_file: The name of the analysis file that describes the named posterior, or an object of class `eos.AnalysisFile`.
    :type analysis_file: str or `eos.AnalysisFile`
    :param figure_name: The name of the figure to draw.
    :type figure_name: str
    :param base_directory: The base directory for the storage of data files. Can also be set via the EOS_BASE_DIRECTORY environment variable.
    :type base_directory: str, optional
    :param format: The file extension of the figure files. Can also be a list of file extensions.
    :type format: str or list of str, optional
    """
    if figure_name not in analysis_file._figures.keys():
        raise ValueError(f'Figure with name \'{figure_name}\' not found in analysis file')

    if isinstance(format, str):
        format = [format]

    SUPPORTED_FORMATS = ['pdf', 'png', 'svg']
    for fmt in format:
        if fmt not in SUPPORTED_FORMATS:
            raise ValueError(f'Figure format \'{fmt}\' not supported. Supported formats are: {",".join(SUPPORTED_FORMATS)}')

    eos.inprogress(f'Drawing figure {figure_name}')
    context = eos.analysis_file_context.AnalysisFileContext(base_directory=base_directory)
    figure = analysis_file._figures[figure_name]
    # the figure must be drawn only once: drawing it once per format would stack its elements on
    # top of each other, thereby distorting the colors of elements with opacity < 1
    figure.draw(context)
    for fmt in format:
        eos.info(f'Saving figure {figure_name} in \'{fmt}\' format')
        figure.save(_staged_file(os.path.join(base_directory, 'figures', f'{figure_name}.{fmt}')))
    eos.completed('... drawing finished')
