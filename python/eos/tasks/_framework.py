# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2022-2026 Danny van Dyk
# Copyright (c) 2024-2025 Matthew Kirk
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
import contextlib
import contextvars
import functools
import inspect
import logging
import os
import shutil
import string
import uuid

from ..ipython import __ipython__


class LogfileHandler:
    """Context manager that redirects EOS log output to a file for the duration of a ``with`` block.

    On entry, a :class:`logging.FileHandler` writing to ``<path>/<name>`` is added to the EOS logger; on
    exit, it is removed again. Used by :func:`task` to capture an invocation-specific log file.

    :param path: The directory in which the log file is created.
    :type path: str
    :param name: The name of the log file. Defaults to ``'log'``.
    :type name: str
    :param mode: The file mode passed to :class:`logging.FileHandler`, e.g. ``'w'`` to truncate or ``'a'`` to append. Defaults to ``'w'``.
    :type mode: str
    """

    def __init__(self, path, name='log', mode='w'):
        self.path = path
        self.name = name
        self.mode = mode
        self.formatter = logging.Formatter('%(asctime)-15s %(levelname)-10s %(message)s')
        self.handler = logging.FileHandler(os.path.join(path, name), mode=mode)
        self.handler.setLevel(logging.DEBUG)
        self.handler.setFormatter(self.formatter)

    def __enter__(self):
        eos.logger.addHandler(self.handler)

    def __exit__(self, type, value, traceback):
        eos.logger.removeHandler(self.handler)
        self.handler.close()


_tasks = {}
_task_outputs = {}
_task_output_fields = {}


def _check_path_component(name, value):
    if value in ('', '.', '..') or any(c in value for c in {'/', '\\', os.sep, os.altsep} - {None}) or any(c.isspace() for c in value):
        raise ValueError(f"Invalid value '{value}' for argument '{name}': expected a single path component without whitespace or backslash")


def _remove(path):
    if os.path.isdir(path) and not os.path.islink(path):
        shutil.rmtree(path)
    else:
        os.remove(path)


class _Staging:
    """Temporary directories and files for a task's outputs, each of which replaces its final counterpart once the task succeeds.

    Files written in place are backed up instead, and restored if the task fails.
    """

    def __init__(self):
        self.directories = {}
        self.files = {}
        self.backups = {}

    @staticmethod
    def _sibling(path, suffix):
        parent, name = os.path.split(path)
        stem, ext = os.path.splitext(name)
        # keep the extension, from which some writers infer the file format
        return os.path.join(parent, f'.{stem}.{uuid.uuid4().hex}.{suffix}{ext}')

    def _within_directory(self, path):
        # a path below a staged directory is written into that directory's temporary counterpart
        for directory, temporary in self.directories.items():
            if path.startswith(directory + os.sep):
                return os.path.join(temporary, os.path.relpath(path, directory))
        return None

    def _check_not_below(self, path):
        if any(other.startswith(path + os.sep) for other in (*self.directories, *self.files, *self.backups)):
            raise RuntimeError(f'Cannot stage directory \'{path}\' after staging a path below it')

    def directory(self, path):
        path = os.path.abspath(path)
        within = self._within_directory(path)
        if within is not None:
            os.makedirs(within, exist_ok=True)
            return within
        if path not in self.directories:
            self._check_not_below(path)
            temporary = self._sibling(path, 'tmp')
            os.makedirs(temporary)
            self.directories[path] = temporary
        return self.directories[path]

    def file(self, path):
        path = os.path.abspath(path)
        within = self._within_directory(path)
        if within is not None:
            os.makedirs(os.path.dirname(within), exist_ok=True)
            return within
        if path not in self.files:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            self.files[path] = self._sibling(path, 'tmp')
        return self.files[path]

    def backup(self, path):
        path = os.path.abspath(path)
        within = self._within_directory(path)
        if within is not None:
            os.makedirs(os.path.dirname(within), exist_ok=True)
            return within
        if path not in self.backups:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            previous = None
            if os.path.lexists(path):
                previous = self._sibling(path, 'old')
                os.rename(path, previous)
            self.backups[path] = previous
        return path

    def commit(self):
        moved, placed = [], []
        try:
            for path, temporary in (*self.directories.items(), *self.files.items()):
                if os.path.lexists(path):
                    previous = self._sibling(path, 'old')
                    os.rename(path, previous)
                    moved.append((path, previous))
                os.rename(temporary, path)
                placed.append(path)
        except BaseException:
            for path in reversed(placed):
                _remove(path)
            for path, previous in reversed(moved):
                os.rename(previous, path)
            self.discard()
            raise
        for _, previous in (*moved, *self.backups.items()):
            if previous is not None:
                _remove(previous)

    def discard(self):
        for temporary in (*self.directories.values(), *self.files.values()):
            if os.path.lexists(temporary):
                _remove(temporary)
        for path, previous in self.backups.items():
            if os.path.lexists(path):
                _remove(path)
            if previous is not None:
                os.rename(previous, path)


_staging = contextvars.ContextVar('_staging', default=None)


def _staged(path):
    """Return the directory into which the running task writes the output destined for ``path``."""
    staging = _staging.get()
    return path if staging is None else staging.directory(path)


def _staged_file(path):
    """Return the file into which the running task writes the output destined for ``path``."""
    staging = _staging.get()
    return path if staging is None else staging.file(path)


def _backed_up(path):
    """Return ``path``, which the running task writes in place; any previous file there is restored if the task fails."""
    staging = _staging.get()
    return path if staging is None else staging.backup(path)


def task_output_templates():
    """Return the registered task output templates without executing any task."""
    return dict(_task_outputs)


def task_output_fields():
    """Return, for each registered task, the callable that derives additional output template fields from the task's arguments."""
    return dict(_task_output_fields)


def expand_task_output(output, output_fields, arguments):
    """Expand a task output template for the given arguments, and check that each interpolated value is a single path component.

    :raises ValueError: If an interpolated value is not a single path component.
    """
    derived = output_fields(**arguments)
    fields = {**arguments, **derived}
    for _, field, spec, _ in string.Formatter().parse(output):
        if field is None:
            continue
        value = format(fields[field], spec)
        if field not in derived or value != '':
            _check_path_component(field, value)
    return output.format(**fields)

def task(name, output, mode=lambda **kwargs: 'w', modules=None, logfile=True, load_analysis_file=True, staged=True,
         output_fields=lambda **kwargs: {}):
    """Decorator that registers a function as a named EOS task.

    The decorated function is wrapped so that, on each invocation, it imports any optional ``modules``,
    resolves a string ``analysis_file`` argument into an :class:`eos.AnalysisFile`, creates the output
    directory, and (unless disabled) captures its log output into a per-invocation log file. Each value
    interpolated into ``output`` must be a single path component (i.e., must not contain white spaces, path
    component delimiters, or relative path specifications like '.' or '..'). When running
    under IPython, the task's log output is additionally shown in a collapsible output widget. The wrapped
    function is registered in the task registry under ``name`` so that it can be looked up and run by the
    ``run`` task and the command-line interface.

    :param name: The name under which the task is registered.
    :type name: str
    :param output: A format string, relative to the base directory, for the task's output directory.
        It may reference the task's arguments (e.g. ``'data/{posterior}/nested'``) and the fields returned by ``output_fields``.
    :type output: str
    :param mode: A callable that, given the task's arguments as keyword arguments, returns the file mode
        for the log file (e.g. ``'w'`` or ``'a'``). Defaults to always returning ``'w'``.
    :type mode: collections.abc.Callable
    :param modules: The names of optional Python modules to import before the task runs. Defaults to none.
    :type modules: list[str] | None
    :param logfile: Whether to capture the task's log output into a log file in the output directory. Defaults to True.
    :type logfile: bool
    :param load_analysis_file: Whether to resolve a string ``analysis_file`` argument into an
        :class:`eos.AnalysisFile` before invoking the task. Defaults to True.
    :type load_analysis_file: bool
    :param staged: Whether the output directory, and any directory the task obtains through :func:`_staged`,
        is written to a temporary directory that replaces the final one only once the task succeeds. A failed
        task thus leaves the final directories untouched. Files obtained through :func:`_staged_file` are staged,
        and those obtained through :func:`_backed_up` are restored on failure, regardless. A task invoked from within another task commits only together with the outermost one.
        Ignored if ``output`` is empty. Defaults to True.
    :type staged: bool
    :param output_fields: A callable that, given the task's arguments as keyword arguments, returns a dict of
        additional fields for ``output``; empty values are permitted. Defaults to returning no fields.
    :type output_fields: collections.abc.Callable
    :returns: A decorator that registers and returns the wrapped task function.
    """
    if modules is None:
        modules = []
    def _task(func):
        @functools.wraps(func)
        def task_wrapper(*args, **kwargs):
            # import task-specific optional modules
            try:
                import importlib as _importlib
                for module in modules:
                    func.__globals__[module] = _importlib.import_module(module)
            except ModuleNotFoundError as e:
                eos.error(f'failed to import missing module \'{module}\': {e}')
                raise e
            except ImportError as e:
                eos.error(f'failed to import module \'{module}\': {e}')
                raise e

            # extract default arguments
            _args = {
                k: v.default
                for k, v in inspect.signature(func).parameters.items()
                if v.default is not inspect.Parameter.empty
            }
            _args.update(zip(func.__code__.co_varnames, args))
            _args.update(kwargs)
            if load_analysis_file and 'analysis_file' in _args and type(_args['analysis_file']) is str:
                _args.update({ 'analysis_file': eos.AnalysisFile(_args['analysis_file'])})
            relative_output = expand_task_output(output, output_fields, _args)
            # a task invoked from within another task joins the outer task's staging
            outer = _staging.get()
            staging = outer if outer is not None else (_Staging() if output else None)
            stagedlog = failedlog = None
            token = _staging.set(staging)
            try:
                # create output directory if needed directly or for logging
                if output or logfile:
                    outputpath = f'{_args["base_directory"]}/{relative_output}'
                    if staging is None or not staged:
                        os.makedirs(outputpath, exist_ok=True)
                    else:
                        finalpath, outputpath = outputpath, staging.directory(outputpath)
                if logfile:
                    logmode = mode(**_args)
                    if staging is not None and staged:
                        stagedlog = os.path.join(outputpath, 'log')
                        parent, basename = os.path.split(os.path.normpath(finalpath))
                        failedlog = os.path.join(parent, f'.{basename}.failed.log')
                        if logmode == 'a' and os.path.isfile(os.path.join(finalpath, 'log')):
                            shutil.copy2(os.path.join(finalpath, 'log'), stagedlog)
                    # create invocation-specific log file handler
                    handler = LogfileHandler(path=outputpath, mode=logmode)
                else:
                    handler = contextlib.nullcontext()
                # use invocation-specific log file handler
                with handler:
                    iaccordion = None
                    ioutput = contextlib.nullcontext()
                    if __ipython__:
                        import ipywidgets as _ipywidgets
                        from IPython.display import display
                        ioutput = _ipywidgets.Output(layout={'height': '200px', 'overflow': 'auto'})
                        iaccordion = _ipywidgets.Accordion(children=[ioutput])
                        iaccordion.set_title(0, relative_output)
                        display(iaccordion)
                    # use invocation-specific ipython output widget (if available)
                    failure = None
                    # ipywidgets.Output suppresses whatever is raised within it, so the
                    # failure is recorded here and re-raised once it has been displayed.
                    with ioutput:
                        try:
                            result = func(**_args)
                        except BaseException as exception:
                            failure = exception
                            raise
                        if iaccordion:
                            iaccordion.selected_index = None

                    if failure is not None:
                        raise failure
            except BaseException:
                # keep the log of a failed task as a hidden sibling of its output directory
                if stagedlog is not None and os.path.isfile(stagedlog):
                    with contextlib.suppress(OSError):
                        shutil.copy2(stagedlog, failedlog)
                if staging is not None and outer is None:
                    staging.discard()
                raise
            finally:
                _staging.reset(token)

            if staging is not None and outer is None:
                staging.commit()
            if failedlog is not None and os.path.lexists(failedlog):
                os.remove(failedlog)

            return result
        _tasks[name] = task_wrapper
        _task_outputs[name] = output
        _task_output_fields[name] = output_fields
        return task_wrapper
    return _task


def _check_varied_parameters_match(analysis: eos.Analysis, data):
    # Check the parameters varied in the analysis match those of the loaded samples
    analysis_varied_params = [p.name() for p in analysis.varied_parameters]
    samples_varied_params = [p["name"] for p in data.varied_parameters]
    if analysis_varied_params != samples_varied_params:
        raise ValueError("Parameters varied in the analysis file don't match those from the loaded sample")


def _mask_output_fields(mask_name=None, **kwargs):
    return { 'mask': '' if mask_name is None else f'_mask-{mask_name}' }
