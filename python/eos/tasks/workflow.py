# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2021-2026 Danny van Dyk
# Copyright (c) 2024-2025 Matthew Kirk
# Copyright (c) 2026      Lorenz Gärtner
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

from ..analysis_file_description import AnalysisFileDescription
from ..diagnostic import Severity, attach_line_numbers
from ..validation_context import ValidationContext
from ._framework import task, _tasks


# Run one analysis step
@task('run', '')
def run(analysis_file:str, id:str, base_directory:str='./', dry_run:bool=False):
    """
    Runs one out of a list of predefined steps recorded in the analysis file.

    Each step corresponds to a call to one or more of the common tasks, e.g.,
     - sample-mcmc
     - find-cluster
     - sample-pmc
     - predict-observables

    :param analysis_file: The name of the analysis file that describes the named posterior, or an object of class `eos.AnalysisFile`.
    :type analysis_file: str or `eos.AnalysisFile`
    :param id: The id of the step to run.
    :type id: str
    :param base_directory: The base directory for the storage of data files. Can also be set via the EOS_BASE_DIRECTORY environment variable.
    :type base_directory: str, optional
    :param dry_run: The flag that disables execution and instead prints the full information on the tasks that would be run to standard output. Defaults to `False`.
    :type dry_run: bool, optional
    """
    try:
        step = analysis_file._steps[id]
    except KeyError:
        raise ValueError(f'Step with id \'{id}\' not found in analysis file')

    for task in step.tasks:
        if dry_run:
            arguments = step.default_arguments[task.task] | task.arguments | { 'analysis_file': analysis_file.analysis_file, 'base_directory': base_directory }
            print(f"{task.task}({','.join([f'{k}={v}' for k, v in arguments.items()])})")
        else:
            arguments = step.default_arguments[task.task] | task.arguments | { 'analysis_file': analysis_file, 'base_directory': base_directory }
            task_function = _tasks[task.task]

            task_function(**arguments)


@task('validate', '', logfile=False, load_analysis_file=False)
def validate(analysis_file:str, deep:bool=True):
    """
    Validates the analysis file structurally and semantically, then by default checks that all
    posteriors and prediction sets can be created. The deep phase inserts declarations into the
    process-wide EOS registries; once it has run, those registries are contaminated for the
    remainder of the process. Use ``deep=False`` for idempotent, side-effect-free validation.

    :param analysis_file: The name of the analysis file to validate.
    :type analysis_file: str
    :param deep: If True, additionally check that all posteriors and prediction sets can be created. Defaults to True.
    :type deep: bool
    """
    description = AnalysisFileDescription.from_yaml_file(analysis_file)
    structural_diagnostics = list(description.validate_structure())
    if any(diagnostic.severity is Severity.ERROR for diagnostic in structural_diagnostics):
        structural_diagnostics = attach_line_numbers(structural_diagnostics, analysis_file)
        for diagnostic in structural_diagnostics:
            print(diagnostic)
        return structural_diagnostics

    semantic_diagnostics = list(description.validate_semantics(ValidationContext(description)))
    diagnostics = structural_diagnostics + semantic_diagnostics

    # Abort before the deep phase if the file is already known to be invalid: constructing an
    # eos.AnalysisFile declares this file's custom parameters and inserts its custom observables
    # into the process-wide EOS registries, which cannot be undone. Running it for a file whose
    # semantics do not hold contaminates the process and yields only cascading errors.
    if not deep or any(diagnostic.severity is Severity.ERROR for diagnostic in semantic_diagnostics):
        diagnostics = attach_line_numbers(diagnostics, analysis_file)
        for diagnostic in diagnostics:
            print(diagnostic)
        return diagnostics

    # Only warnings remain at this point. Print them here, since eos.AnalysisFile.validate() reports
    # the semantic and deep phases but knows nothing about the structural one.
    structural_diagnostics = attach_line_numbers(structural_diagnostics, analysis_file)
    for diagnostic in structural_diagnostics:
        print(diagnostic)

    return structural_diagnostics + eos.AnalysisFile(analysis_file).validate(deep=True)


@task('list-figures', '', logfile=False)
def list_figures(analysis_file:str):
    """
    Lists all figures that can be created from the analysis file.

    :param analysis_file: The name of the analysis file, or an object of class `eos.AnalysisFile`.
    :type analysis_file: str or :class:`eos.AnalysisFile`
    """
    return [name for name in analysis_file._figures.keys()]


@task('list-steps', '', logfile=False)
def list_steps(analysis_file:str):
    """
    Lists the steps in the analysis file.

    :param analysis_file: The name of the analysis file that shall be inspected`.
    :type analysis_file: str or :class:`eos.AnalysisFile`
    """
    return analysis_file._steps.keys()


@task('list-step-dependencies', '', logfile=False)
def list_step_dependencies(analysis_file:str, id:str):
    """
    Lists all steps required to be completed before the given step can be executed.

    :param analysis_file: The name of the analysis file that shall be inspected`.
    :type analysis_file: str or :class:`eos.AnalysisFile`
    :param id: The id of the step to inspect.
    :type id: str
    """
    steps = analysis_file._steps

    if id not in steps:
        raise ValueError(f'Step with id \'{id}\' not found in analysis file')

    return steps[id].depends_on
