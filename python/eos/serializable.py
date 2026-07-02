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

import dataclasses as _dc
import yaml as _yaml


# Top-level packages whose objects must never be serialized. On-disk descriptions must hold only
# native Python types; a value from one of these packages (e.g. a numpy array or scalar) is rejected
# rather than being silently coerced or emitted with a library-specific YAML tag. Extend this list to
# exclude further libraries.
_EXCLUDED_MODULES = ['numpy']


class _SimpleDumper(_yaml.SafeDumper):
    """A :class:`yaml.SafeDumper` restricted to simple, native Python types.

    Any object whose defining top-level module is listed in :data:`_EXCLUDED_MODULES` is rejected with a
    clear error instead of being serialized. The check is by module name, so it needs no import of the
    excluded libraries and catches values at any depth of the structure.
    """
    def represent_data(self, data):
        top_module = type(data).__module__.split('.', 1)[0]
        if top_module in _EXCLUDED_MODULES:
            raise TypeError(
                f'refusing to serialize the object of type \'{type(data).__module__}.{type(data).__name__}\' '
                f'from the excluded module \'{top_module}\'; convert it to a native Python type '
                f'(e.g. via float(), int(), or .tolist()) first'
            )
        return super().represent_data(data)


class Serializable:
    """
    Utility mixin to serialize instances of classes to YAML data.

    This is the inverse of :class:`eos.deserializable.Deserializable`. Classes that participate in
    write/read symmetry inherit from both; classes that are only ever read from disk inherit
    :class:`Deserializable` alone.
    """
    def to_dict(self):
        """Serialize this instance into a dictionary suitable for YAML output.

        The default implementation requires the instance to be a :func:`dataclasses.dataclass` and
        returns :func:`dataclasses.asdict`. Classes whose on-disk representation differs from their
        field layout (e.g. keys that are not valid Python identifiers) override this method to emit
        the canonical on-disk keys. It is the inverse of :meth:`Deserializable.from_dict`.

        :returns: The dictionary representation of this instance.
        :rtype: dict
        """
        if not _dc.is_dataclass(self):
            raise TypeError(f'{type(self).__name__} is not a dataclass and must override to_dict()')
        return _dc.asdict(self)

    def to_yaml_file(self, path:str):
        """Serialize this instance via :meth:`to_dict` and write it to a YAML file.

        :param path: Path to the YAML file to write.
        :type path: str
        """
        with open(path, 'w') as f:
            _yaml.dump(self.to_dict(), f, Dumper=_SimpleDumper, default_flow_style=False)
