# Copyright (c) 2023-2024 Danny van Dyk
# Copyright (c) 2023 Philip Lueghausen
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

import copy as _copy
import eos
import numpy as _np
import os
import yaml as _yaml


@dataclass(kw_only=True)
class Item(Deserializable):
    r"""Base class for items to be drawn in a figure."""

    color:str=field(default='black')
    label:str=field(default=None)
    linestyle:str=field(default='solid')

    def prepare(self):
        "Prepare the item for drawing"
        raise NotImplementedError

    def draw(self, ax, **kwargs):
        "Draw the item on the axes"
        raise NotImplementedError


@dataclass(kw_only=True)
class ObservableItem(Item):
    r"""Represents an observable as a function of a variable.
    The variable can be either a kinematic variable or a parameter.

    :param observable: EOS qualified name of an observable
    :type observable: str
    :param variable: Name of a kinematic variable or of a parameter
    :type variable: str
    :param xvalues: Values to be used as sampling points
    :type xvalues: list of float numbers
    :param fixed_kinematics: Names and values of fixed kinematic variables
    :type fixed_kinematics: dict
    :param fixed_parameters: Names and values of fixed parameters
    :type fixed_parameters: dict
    :param fixed_parameters_from_file: Path to file that contains names and values of fixed parameters in the YAML format
    :type fixed_parameters_from_file: str
    :param options: Names and values of options passed to the EOS observable
    :type options: dict
    :param label: Label used in the legend
    :type label: str
    """

    observable:eos.QualifiedName
    variable:str
    xvalues:list[float]
    fixed_kinematics:dict=None
    fixed_parameters:dict=None
    fixed_parameters_from_file:str=None
    options:dict=None
    label:str=None

    def __post_init__(self):
        eos.info(f'Handling item to plot {self.observable}')
        self._observable_entry = eos.Observables._get_obs_entry(self.observable)
        valid_kinematic_variables = {kv for kv in self._observable_entry.kinematic_variables()}

        # Create kinematics
        self._kinematics = eos.Kinematics()
        if self.fixed_kinematics:
            for k, v in self.fixed_kinematics.items():
                if k not in valid_kinematic_variables:
                    raise ValueError(f"Kinematic variable '{k}' does not match any of the declared kinematic variables '{self.observable}': {valid_kinematic_variables.__repr__()}")
                self._kinematics.declare(k, v)

        # Create parameters
        self._parameters = eos.Parameters.Defaults()
        if self.fixed_parameters_from_file:
            eos.warn('Overriding parameters from file')
            self._parameters.override_from_file(self.fixed_parameters_from_file)

        if self.fixed_parameters:
            if self.fixed_parameters_from_file:
                eos.warn('Overriding values read from \'parameters-from-file\' with explicit values in \'parameters\'')

            for key, value in self.fixed_parameters.items():
                self._parameters.set(key, value)

        # Declare variable that is either parameter or kinematic
        self._variable = None

        # Does the variable correspond to one of the kinematic variables?
        if self.variable in valid_kinematic_variables:
            self._variable = self._kinematics.declare(self.variable, _np.nan)
        else:
            # Is the variable name a QualifiedName?
            try:
                qn = eos.QualifiedName(self.variable)
                # Continues only if no failure occures
                self._variable = self.parameters[qn]
                self._variable.set(_np.nan)
            except RuntimeError:
                raise ValueError(f"Value of 'variable' for observable '{self.observable}' is neither a valid kinematic variable nor a valid parameter: '{self.variable}")

        # Create options
        self._options = eos.Options()
        if self.options:
            for key, value in self.options.items():
                self._options.declare(key, value)

        self._observable = eos.Observable.make(self.observable, self._parameters, self._kinematics, self._options)


    def prepare(self):
        "Evaluate the observable at the sample points"
        self.yvalues = _np.empty((len(self.xvalues),))
        for i, x in enumerate(self.xvalues):
            self._variable.set(x)
            self.yvalues[i] = self._observable.evaluate()


    def draw(self, ax, **kwargs):
        "Draw a line plot of the observable"
        ax.plot(self.xvalues, self.yvalues, label=self.label, **kwargs)


@dataclass(kw_only=True)
class OneDimensionalHistogramItem(Item):
    r"""Represents a one-dimensional histogram.

    :param datafile: Path to the file that contains the histogram data.
    :type datafile: str
    :param variable: Name of the datafile's variable that shall be histogrammed.
    :type variable: str
    :param bins: Number of histogram bins.
    :type bins: int
    """

    datafile:str
    variable:str
    bins:int=field(default=50)

    def __post_init__(self):
        os.path.exists(self.datafile) or eos.error(f"Data file '{self.datafile}' does not exist")
        abspath = os.path.abspath(self.datafile)
        name = os.path.split(abspath)[-1]
        if name == 'samples':
            self._datafile = eos.data.ImportanceSamples(abspath)
        elif name.startswith('pred-'):
            self._datafile = eos.data.Prediction(abspath)
        else:
            eos.error(f"Data file '{self.datafile}' has an unsupported format")
            raise NotImplementedError

        if self.variable not in self._datafile.lookup_table:
            raise ValueError(f"Data file '{self.datafile}' does not contain samples of variable '{self.variable}'")

        if self.bins < 2:
            raise ValueError(f"Number of bins '{self.bins}' is smaller than 2")

    def prepare(self):
        "Prepare the histogram for drawing"
        pass

    def draw(self, ax):
        "Draw the histogram"
        idx = self._datafile.lookup_table[self.variable]
        ax.hist(self._datafile.samples[:, idx], weights=self._datafile.weights, bins=self.bins, label=self.label)

@dataclass
class TwoDimensionalHistogramItem(Item):
    r"""Represents a two-dimensional histogram.

    :param datafile: Path to the file that contains the histogram data.
    :type datafile: str
    :param variable: Name of the datafile's variable that shall be histogrammed.
    :type variables: tuple[str, str]
    :param bins: Number of histogram bins.
    :type bins: int
    """

    datafile:str
    variables:tuple[str, str]
    bins:int=field(default=50)

    def __post_init__(self):
        os.path.exists(self.datafile) or eos.error(f"Data file '{self.datafile}' does not exist")
        abspath = os.path.abspath(self.datafile)
        name = os.path.split(abspath)[-1]
        if name == 'samples':
            self._datafile = eos.data.ImportanceSamples(abspath)
        elif name.startswith('pred-'):
            self._datafile = eos.data.Prediction(abspath)
        else:
            eos.error(f"Data file '{self.datafile}' has an unsupported format")
            raise NotImplementedError

        for v in self.variables:
            if v not in self._datafile.lookup_table:
                raise ValueError(f"Data file '{self.datafile}' does not contain samples of variable '{v}'")

        if self.bins < 2:
            raise ValueError(f"Number of bins '{self.bins}' is smaller than 2")

    def prepare(self):
        pass

    def draw(self, ax):
        xidx = self._datafile.lookup_table[self.variables[0]]
        yidx = self._datafile.lookup_table[self.variables[1]]
        ax.hist2d(self._datafile.samples[:, xidx], self._datafile.samples[:, yidx],
                  bins=self.bins, cmin=1, label=self.label)


class ItemFactory:
    registry = {
        'observable': ObservableItem,
        'histogram1D': OneDimensionalHistogramItem,
        'histogram2D': TwoDimensionalHistogramItem,
    }

    @staticmethod
    def from_yaml(yaml_data:str):
        kwargs = _yaml.safe_load(yaml_data)
        return ItemFactory.from_dict(**kwargs)

    @staticmethod
    def from_dict(**kwargs):
        "Factory method to create an item from a dictionary"
        if 'type' not in kwargs:
            raise ValueError(f"Content item { kwargs } does not contain element 'type'")

        if kwargs['type'] not in ItemFactory.registry:
            raise ValueError(f'Unknown content item type: {kwargs["type"]}')

        item_type = kwargs.pop('type')

        return ItemFactory.registry[item_type].from_dict(**kwargs)
