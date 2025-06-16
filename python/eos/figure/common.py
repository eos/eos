# Copyright (c) 2025 Danny van Dyk
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

from dataclasses import dataclass, field
from eos.deserializable import Deserializable

import numpy as np

@dataclass(kw_only=True)
class Range(Deserializable):
    r""" Collects the relevant information to generate a range of floating point values.

    :param min: Minimal value.
    :type min: float
    :param max: Maximal value.
    :type max: float
    :param num: Number of values to generate in the range.
    :type num: int
    """
    min:float
    max:float
    num:int

    def __post_init__(self):
        if self.min >= self.max:
            eos.error(f"Range min ({self.min}) must be smaller than max ({self.max})")
            raise ValueError("Invalid range")

    @property
    def values(self):
        r""" Generate the range of values.

        :return: List of floating point values.
        :rtype: list[float]
        """
        return np.linspace(self.min, self.max, self.num).tolist()
