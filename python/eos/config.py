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

# The '_eos_data' package only exists in a wheel install (it carries the
# reference/constraint/parameter YAML files bundled with the wheel); its
# presence at import time is therefore what distinguishes a wheel install
# from a source install. Use find_spec() rather than an actual import so
# that checking for it never triggers loading it.
import importlib.util as _importlib_util
is_wheel = _importlib_util.find_spec('_eos_data') is not None
