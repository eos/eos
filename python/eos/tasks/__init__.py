# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2020-2026 Danny van Dyk
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


from ._framework import LogfileHandler, task, task_output_templates, task_output_fields, expand_task_output
from ._framework import _tasks, _task_outputs, _task_output_fields, _staged, _staged_file, _backed_up, _mask_output_fields
from .inference import find_mode, sample_prior, DynestyResultLogger, sample_nested
from .pmc import sample_mcmc, find_clusters, mixture_product, sample_pmc
from .predictions import predict_observables, create_mask, GoodnessOfFitResult, goodness_of_fit_tests, ks_test, create_constraint
from .comparison import model_comparison, _bayes_factor_strength
from .output import report, draw_figure
from .workflow import run, validate, list_figures, list_steps, list_step_dependencies

__all__ = [
    'LogfileHandler', 'task', 'task_output_templates', 'task_output_fields', 'expand_task_output',
    'find_mode', 'sample_prior', 'DynestyResultLogger', 'sample_nested',
    'sample_mcmc', 'find_clusters', 'mixture_product', 'sample_pmc',
    'predict_observables', 'create_mask', 'GoodnessOfFitResult', 'goodness_of_fit_tests', 'ks_test', 'create_constraint',
    'model_comparison',
    'report', 'draw_figure',
    'run', 'validate', 'list_figures', 'list_steps', 'list_step_dependencies',
]
