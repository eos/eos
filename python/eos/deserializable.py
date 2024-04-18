# Copyright (c) 2024 Danny van Dyk
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
        return cls.from_dict(**kwargs)

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
