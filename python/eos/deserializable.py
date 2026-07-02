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

import os as _os
import yaml as _yaml

class Deserializable:
    """
    Utility class to create instances of classes from YAML data.
    """
    def __init__(self):
        pass

    @classmethod
    def from_yaml(cls, yaml_data:str):
        kwargs = _yaml.safe_load(yaml_data)
        if not isinstance(kwargs, dict):
            raise ValueError('The provided YAML data does not describe a mapping')
        return cls.from_dict(**kwargs)

    @classmethod
    def from_yaml_file(cls, path:str):
        """Load and parse a single YAML file, then deserialize it via :meth:`from_dict`.

        Centralizes the checks and error handling shared by all on-disk descriptions: the file must
        exist and contain a top-level YAML mapping. The parsed mapping is forwarded to :meth:`from_dict`,
        so any per-class normalization and validation applies unchanged.

        :param path: Path to the YAML file to load.
        :type path: str
        :returns: The deserialized instance.
        :raises RuntimeError: If the file does not exist, is not a file, or does not contain a mapping.
        """
        if not _os.path.exists(path) or not _os.path.isfile(path):
            raise RuntimeError(f'Description file {path} does not exist or is not a file')
        with open(path) as f:
            data = _yaml.safe_load(f)
        if not isinstance(data, dict):
            raise RuntimeError(f'Description file {path} does not contain a YAML mapping')
        return cls.from_dict(**data)

    @staticmethod
    def make(cls, **kwargs):
        if not isinstance(cls, type):
            raise ValueError(f'Argument cls=\'{cls}\' is not a type')

        if not issubclass(cls, Deserializable):
            raise ValueError(f'Class {cls.__name__} is not a subclass of Deserializable')

        try:
            result = cls(**kwargs)
        except Exception as e:
            raise ValueError(f'When creating {cls.__name__} from {kwargs}: {e}')

        return result

    @classmethod
    def from_dict(cls, **kwargs):
        return Deserializable.make(cls, **kwargs)
