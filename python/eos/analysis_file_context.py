#!/usr/bin/python
# vim: set sw=4 sts=4 et tw=120 :

# Copyright (c) 2024-2026 Danny van Dyk
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

from dataclasses import dataclass
import os

@dataclass(kw_only=True)
class AnalysisFileContext:
    """
    Auxiliary class used to resolve paths relative to an analysis file's base directory.

    Instances are passed to consumers of an analysis file (e.g. figures and tasks) so that relative data
    paths can be resolved against a common base directory.

    :param base_directory: The base directory against which relative paths are resolved. Defaults to the current working directory.
    :type base_directory: str
    """
    base_directory:str='./'

    def __post_init__(self):
        if not os.path.exists(self.base_directory):
            raise ValueError(f'Base directory \'{self.base_directory}\' does not exist')

        if not os.path.isdir(self.base_directory):
            raise ValueError(f'Base directory \'{self.base_directory}\' is not a directory')

    def data_path(self, relative_path:str):
        """Resolve a path relative to :attr:`base_directory` into an absolute path.

        :param relative_path: A path relative to the base directory.
        :type relative_path: str
        :returns: The corresponding absolute path.
        :rtype: str
        """
        return os.path.abspath(os.path.join(self.base_directory, relative_path))
