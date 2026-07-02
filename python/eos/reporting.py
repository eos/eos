# Copyright (c) 2026 Danny van Dyk
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
import os as _os
import glob as _glob
from dataclasses import dataclass
from functools import cached_property


class PosteriorData:
    r"""Provides lazy, read-only access to all data recorded for a single posterior.

    An instance corresponds to one posterior directory below ``EOS_BASE_DIRECTORY/data/``, e.g.
    ``data/CKM/``. It is the single object exposed per posterior to the Jinja2-based reporting
    framework. All recorded data (importance samples, nested-sampling results, (local) modes and their
    goodness-of-fit, and observable predictions) is loaded lazily on first access and cached, so a
    report only ever reads the files it actually references.

    :ivar name: The name of the posterior, i.e. the name of its directory below ``data/``.
    :ivar base_directory: The base directory below which the ``data/`` tree is located.

    :param name: The name of the posterior.
    :type name: str
    :param base_directory: The base directory for the storage of data files. Defaults to ``'./'``.
    :type base_directory: str, optional
    :param analysis_file: The analysis file that describes the posterior. Optional; when provided it
        allows access to the posterior's description and to a freshly created :class:`eos.Analysis`.
    :type analysis_file: eos.AnalysisFile, optional
    """

    def __init__(self, name:str, base_directory:str='./', analysis_file=None):
        self.name           = name
        self.base_directory = base_directory
        self._analysis_file = analysis_file
        self._path          = _os.path.join(base_directory, 'data', name)

    def __repr__(self):
        return f'PosteriorData(name={self.name!r}, base_directory={self.base_directory!r})'

    # -- recorded distributions -------------------------------------------------------------------

    @cached_property
    def samples(self):
        """The recorded importance samples, or ``None`` if none were recorded.

        :rtype: eos.data.ImportanceSamples | None
        """
        path = _os.path.join(self._path, 'samples')
        if not _os.path.isdir(path):
            return None
        return eos.data.ImportanceSamples(path)

    @cached_property
    def nested(self):
        """The recorded nested-sampling results, or ``None`` if none were recorded.

        :rtype: eos.data.DynestyResults | None
        """
        path = _os.path.join(self._path, 'nested')
        if not _os.path.isdir(path):
            return None
        return eos.data.DynestyResults(path)

    @cached_property
    def predictions(self):
        """The recorded observable predictions, keyed by their label.

        A prediction recorded below ``data/<name>/pred-<label>`` is exposed under the key ``<label>``.

        :rtype: dict[str, eos.data.Prediction]
        """
        result = {}
        for pred_dir in sorted(_glob.glob(_os.path.join(self._path, 'pred-*'))):
            if not _os.path.isdir(pred_dir):
                continue
            label = _os.path.basename(pred_dir)[len('pred-'):]
            result[label] = eos.data.Prediction(pred_dir)
        return result

    # -- (local) modes and goodness-of-fit --------------------------------------------------------

    @cached_property
    def modes(self):
        """The recorded (local) modes, keyed by their label.

        A mode recorded below ``data/<name>/mode-<label>`` is exposed under the key ``<label>``.

        :rtype: dict[str, eos.data.Mode]
        """
        result = {}
        for mode_dir in sorted(_glob.glob(_os.path.join(self._path, 'mode-*'))):
            if not _os.path.isdir(mode_dir):
                continue
            label = _os.path.basename(mode_dir)[len('mode-'):]
            result[label] = eos.data.Mode(mode_dir)
        return result

    @cached_property
    def mode(self):
        """The default mode, for convenient access from a report.

        Prefers the mode labelled ``'default'``; if no such mode exists but exactly one mode was
        recorded, that mode is used. Returns ``None`` if no unambiguous choice can be made.

        :rtype: eos.data.Mode | None
        """
        modes = self.modes
        if 'default' in modes:
            return modes['default']
        if len(modes) == 1:
            return next(iter(modes.values()))
        return None

    @property
    def has_mode(self):
        """Whether an unambiguous mode is available.

        :rtype: bool
        """
        return self.mode is not None

    # -- link back to the analysis file -----------------------------------------------------------

    @cached_property
    def description(self):
        """The description of the posterior taken from the analysis file, if available.

        :rtype: eos.analysis_file_description.PosteriorDescription | None
        """
        if self._analysis_file is None:
            return None
        return self._analysis_file.posteriors.get(self.name, None)

    @cached_property
    def analysis(self):
        """A freshly created :class:`eos.Analysis` for the posterior, if an analysis file is available.

        The analysis is only constructed on first access, so a report that never references it pays no
        cost for building it.

        :rtype: eos.Analysis | None
        """
        if self._analysis_file is None or self.name not in self._analysis_file.posteriors:
            return None
        return self._analysis_file.analysis(self.name)

    # -- figures ----------------------------------------------------------------------------------

    def corner_figure(self, output_directory:str, variables:list=None, format:str='pdf', label:str=None):
        r"""Draw a corner figure of the recorded importance samples and return its relative path.

        Reuses :class:`eos.figure.CornerFigure` to produce a triangular arrangement of the 1D and 2D
        marginal posterior densities. The figure is written to ``<output_directory>/corner-<name>.<format>``
        and the path relative to ``output_directory`` is returned, ready to be embedded in a report.

        :param output_directory: The directory into which the figure is written. Created if necessary.
        :type output_directory: str
        :param variables: The names of the parameters to show. Defaults to ``None``, in which case all
            recorded parameters are shown.
        :type variables: list[str], optional
        :param format: The file extension / format of the figure. Defaults to ``'pdf'``.
        :type format: str, optional
        :param label: The label used for the samples in the figure. Defaults to the posterior's name.
        :type label: str, optional
        :returns: The path to the figure relative to ``output_directory``.
        :rtype: str
        :raises RuntimeError: If no importance samples were recorded for the posterior.
        """
        if self.samples is None:
            raise RuntimeError(f'No importance samples recorded for posterior \'{self.name}\'; cannot draw a corner figure')

        from eos.analysis_file_description import AnalysisFileContext

        _os.makedirs(output_directory, exist_ok=True)

        filename = f'corner-{self.name}.{format}'
        output   = _os.path.join(output_directory, filename)

        figure = eos.figure.FigureFactory.from_dict(
            type='corner',
            contents=[{
                'path':  _os.path.join('data', self.name, 'samples'),
                'label': label if label is not None else self.name,
            }],
            variables=variables,
        )
        figure.draw(context=AnalysisFileContext(base_directory=self.base_directory), output=output)

        return filename


class AnalysisData:
    r"""Provides lazy, disk-driven access to all data recorded below ``EOS_BASE_DIRECTORY/data/``.

    An instance is the single entry point exposed to the Jinja2-based reporting framework. Its
    :attr:`posteriors` are discovered by scanning the ``data/`` directory for subdirectories; each
    such directory is exposed as a :class:`PosteriorData`. The scan itself reads no data files -- each
    :class:`PosteriorData` loads its own contents lazily on first access -- so a report only ever reads
    what it references.

    The container is dict-like: iterating over it yields the :class:`PosteriorData` instances, ``len``
    reports their number, and indexing/containment tests operate on posterior names.

    :ivar base_directory: The base directory below which the ``data/`` tree is located.

    :param base_directory: The base directory for the storage of data files. Defaults to ``'./'``.
    :type base_directory: str, optional
    :param analysis_file: The analysis file that describes the posteriors. Optional; when provided it
        is forwarded to each :class:`PosteriorData` so that descriptions and analyses become available.
    :type analysis_file: eos.AnalysisFile, optional
    """

    def __init__(self, base_directory:str='./', analysis_file=None):
        self.base_directory = base_directory
        self._analysis_file = analysis_file
        self._data_directory = _os.path.join(base_directory, 'data')

    def __repr__(self):
        return f'AnalysisData(base_directory={self.base_directory!r}, posteriors={self.names!r})'

    @cached_property
    def posteriors(self):
        """The recorded posteriors, keyed by name and discovered from the ``data/`` directory.

        Every immediate subdirectory of ``data/`` is treated as a posterior and exposed as a
        :class:`PosteriorData`. Returns an empty mapping if ``data/`` does not exist.

        :rtype: dict[str, PosteriorData]
        """
        result = {}
        if not _os.path.isdir(self._data_directory):
            return result
        for entry in sorted(_os.listdir(self._data_directory)):
            if not _os.path.isdir(_os.path.join(self._data_directory, entry)):
                continue
            result[entry] = PosteriorData(entry, base_directory=self.base_directory, analysis_file=self._analysis_file)
        return result

    @property
    def names(self):
        """The names of the recorded posteriors, in sorted order.

        :rtype: list[str]
        """
        return list(self.posteriors.keys())

    def __iter__(self):
        return iter(self.posteriors.values())

    def __len__(self):
        return len(self.posteriors)

    def __getitem__(self, name):
        return self.posteriors[name]

    def __contains__(self, name):
        return name in self.posteriors
