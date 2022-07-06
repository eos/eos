import eos
import logging
import unittest

logger = logging.getLogger(__name__)

class ImplementationError(Exception):
    pass

class BaseSpec:
    """
    Represents the specification of something.
    It contains a validation accompanied by a docstring.
    The Specification does not need to be initialised with a value,
    but a value can be set and the object tracks if the value was
    used.
    """

    def __init__(self, name, default = None, optional = False):

        self.name = name
        self.default = default
        self.was_used = False
        self.optional = optional
        self.owner = None # a SpecCollection

        self._is_set = False
        self._value = None
        # might want to force default value other than None of spec is optional

    @property
    def value(self):
        # is treated like an attribute via property dec
        self.was_used = True
        if (self._is_set == False) and (self.optional == True):
            return self.default
        else:
            return self._value

    @value.setter
    def value(self, x):
        self._value = x
        self._is_set = True

    def validate(self):
        """
        Raises error if value is not valid.
        The error should suggest valid values.

        The docstring describes the value and appears in the online documentation.
        """
        # perform checks
        if (self._is_set == False) and (self.optional == False):
            raise ImplementationError("Cannot validate value that is not set")
        pass


class BaseSpecTest(unittest.TestCase):

    def test_was_used(self):
        s = BaseSpec(name ='name')
        s.value = 3.1415

        self.assertFalse(s.was_used)
        value = s.value
        self.assertIs(value, 3.1415)
        self.assertTrue(s.was_used)

    def test_default(self):
        ref = 'hi mom'
        s = BaseSpec(name='name', optional=True, default=ref)
        # don't set value
        self.assertEqual(s.value, ref)

    def test_validate_before_set(self):
        s = BaseSpec(name='name')
        with self.assertRaises(ImplementationError):
            s.validate()


class SpecCollection:
    """
    Represents a collection of specifications.

    This should be used with a context manager:
    ```
    with mySpecs as specs:
        ...use specs...
    # checked that specs were all accessed
    ```
    """

    def __init__(self, template):
        """
        template must be a list of Specifications derived from BaseSpec,
        like [BaseSpec( ), BaseSpec( ),   ]
        """
        self.specs_dict = {spec.name: spec for spec in template}

        # set self as owner to provide Specs access to each other
        for _, spec in self.specs_dict.items():
            spec.owner = self

    def set_and_validate(self, values):
        "Validate and set user-provided values"

        # check if all mand. keys are provided
        mand_keys_missing = []
        for name in self.names:
            if not name in values:
                if self[name].optional == False:
                    mand_keys_missing.append(name)
        if not mand_keys_missing == []:
            raise ValueError(f"Mandatory keys are missing: {mand_keys_missing}")

        # set all provided values
        for name in values.keys():
            if name in self.names:
                self[name].value = values[name]

        # validate all values afterwards;
        # dependencies can be injected since the Specs a SpecCollection
        # owns are aware of the owner
        for _, spec in self.specs_dict.items():
            spec.validate()

    def __enter__(self):
        # want to use context manager like: `with SpecCollection as specs``
        return self

    def __exit__(self, type, value, traceback):
        if type is None: # no exception occured within the context manager
            # check if implementation actually uses all parameters
            unused_specs = []
            for name, spec in self.specs_dict.items():
                if spec.was_used == False:
                    unused_specs.append(name)
            if not unused_specs == []:
                raise ImplementationError(f"Specifications were not used: {unused_specs}")
        else: # re-raise occured exceptions
            return False

    @property
    def names(self):
        return self.specs_dict.keys()

    def __getitem__(self, key):
        try:
            return self.specs_dict[key]
        except KeyError as e:
            raise ImplementationError(f"Unknown spec: {key}")


class SpecCollectionTest(unittest.TestCase):

    def test_access_unknown_key(self):
        "Get unknown spec"
        template = [
            BaseSpec(name=f'mand0', optional=False)
        ]
        s = SpecCollection(template)

        with self.assertRaises(ImplementationError):
            s['something'].value

    def test_valid_set_values(self):
        values = {
            'key0': 1.0,
            'key1': 2.0,
        }
        template = [
            BaseSpec(name='key0'),
            BaseSpec(name='key1'),
            BaseSpec(name='opt', optional=True),
        ]
        s = SpecCollection(template)
        s.set_and_validate(values) # should not raise

    def test_invalid_mand_keys_provided(self):
        values = {
            'mand0': None,
            'opt': None,
        }
        template = [
            BaseSpec(name=f'mand0', optional=False),
            BaseSpec(name=f'mand1', optional=False),
            BaseSpec(name=f'mand2', optional=False),
            BaseSpec(name=f'opt', optional=True),
        ]
        s = SpecCollection(template)
        try:
            s.set_and_validate(values)
        except Exception as e:
            self.assertEqual(str(e), "Mandatory keys are missing: ['mand1', 'mand2']")

    def test_valid_mand_keys_provided(self):
        values = {
            'mand0': None,
            'mand1': None,
            'mand2': None,
        }
        template = [
            BaseSpec(name=f'mand0', optional=False),
            BaseSpec(name=f'mand1', optional=False),
            BaseSpec(name=f'mand2', optional=False),
            BaseSpec(name=f'opt', optional=True, default="hi"),
        ]
        s = SpecCollection(template)
        s.set_and_validate(values) # should not raise

    def test_all_specs_used(self):
        "Access all mandatory values: must not raise"
        template = (
            BaseSpec(name='observable'),
            BaseSpec(name='variable'),
        )
        s = SpecCollection(template)
        with s as specs:
            for name in specs.names:
                specs[name].value
        # should not raise since all values were accessed

    def test_not_all_specs_used(self):
        "Leave mandatory spec unaccessed: must raise"
        template = (
            BaseSpec(name='observable'),
            BaseSpec(name='variable'),
        )
        s = SpecCollection(template)
        with self.assertRaises(ImplementationError):
            with s as specs:
                for name in specs.names:
                    pass # do not use spec, must raise


class ObservableSpec(BaseSpec):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def validate(self):
        super().validate()
        "``string`` -- The name of the observable; see `the complete list of observables <../observables.html>`_."
        try:
            # raises if the name is invalid
            eos.Observables._get_obs_entry(self._value)
            pass
        except Exception as e:
            raise e


class ObservableSpecTest(unittest.TestCase):

    def test_valid_name(self):
        s = ObservableSpec(name = 'name')
        s.value = 'B->Dlnu::dBR/dq2;l=mu'
        s.validate()
        # should not raise error

    def test_invalid_name(self):
        s = ObservableSpec(name = 'name')
        s.value = 'Not->Valid::Observable'
        with self.assertRaises(ValueError):
            s.validate()

    def test_value_not_set(self):
        s = ObservableSpec(name = 'name')
        with self.assertRaises(ImplementationError):
            s.validate()


class VariableSpec(BaseSpec):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._type = None # updated on validate with either 'kinematic' or 'parameter'

    def validate(self):
        "``string`` -- Must be a valid kinematic variable or a parameter."

        super().validate()

        try:
            obs_name = self.owner['observable'].value
        except:
            raise ImplementationError("'VariableSpec' requires 'observable' to part of the SpecCollection")

        obs_entry = eos.Observables._get_obs_entry(obs_name)
        valid_kin_vars = [kv for kv in obs_entry.kinematic_variables()]

        if self._value in valid_kin_vars:
            self._type = 'kinematic'

        else: # if not, is the variable name a QualifiedName?
            try:
                qn = eos.QualifiedName(self._value)
                # Does not raise? Then
                self._type = 'parameter'

            except RuntimeError as e:
                raise ValueError(f"Value of 'variable' for observable '{obs_name}'"
                                 f" is neither a valid kinematic variable nor parameter: '{self._value}'") from e

    @property
    def type(self):
        if self._type is None:
            raise RuntimeError("No value set")
        return self._type


class VariableSpecTest(unittest.TestCase):

    def setUp(self):
        self.s = SpecCollection(template=[
            VariableSpec(name = 'variable'),
            ObservableSpec(name = 'observable')
        ])

    def test_valid_kinematic(self):
        item = {
            'observable': 'B->Dlnu::dBR/dq2;l=mu',
            'variable': 'q2'
        }
        self.s.set_and_validate(item)
        self.assertIs(self.s['variable'].type, 'kinematic')

    def test_invalid_kinematic(self):
        item = {
            'observable': 'B->Dlnu::dBR/dq2;l=mu',
            'variable': 'q3'
        }

        with self.assertRaises(ValueError):
            self.s.set_and_validate(item)

    def test_valid_parameter(self):
        item = {
            'observable': 'B->Dlnu::dBR/dq2;l=mu',
            'variable': 'mass::tau'
        }
        self.s.set_and_validate(item)
        self.assertIs(self.s['variable'].type, 'parameter')

    @unittest.skip("TO DO: Should this raise an error?") # The pattern is correct, but the parameter is not defined.
    def test_invalid_parameter(self):
        item = {
            'observable': 'B->Dlnu::dBR/dq2;l=mu',
            'variable': 'mass::feather'
        }
        self.s.set_and_validate(item)
        self.assertIs(self.s['variable'].type, 'parameter')


# class KinematicsSpec(BaseSpec):

#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)
#         self.observable_name = None # Needs to be set after construction, before validate_set

#     def validate_set(self, value):
#         "``dict`` like ``{kinematic: value}`` -- ``kinematic`` must be a valid name, ``value`` must be a float"
#         if self.observable_name == None:
#             raise ImplementationError("Attribute 'observable_name' not set")

#         obs_name = self.observable_name
#         obs_entry = eos.Observables._get_obs_entry(obs_name)
#         valid_kin_vars = [kv for kv in obs_entry.kinematic_variables()]

#         if not isinstance(value, dict):
#             raise ValueError(f"'{self.name}' must be a dictionary of form {{'kinematic_name': value}}")

#         invalid_names = []
#         for name in value.keys():
#             if name not in valid_kin_vars:
#                 invalid_names.append(name)
#                 # might also want to test that arguments are floats

#         if not invalid_names == []:
#             raise ValueError(
#                     f"Kinematic quantity '{name}' does not match known ones "
#                     f"for observable '{obs_name}': {valid_kin_vars}")
#         else:
#             self._value = value


# class KinematicsSpecTest(unittest.TestCase):

#     def test_valid(self):
#         s = KinematicsSpec(name = 'name')
#         ref = {'q2': 1.0}
#         s.observable_name = 'B->Dlnu::dBR/dq2;l=mu'
#         s.validate_set(ref)
#         self.assertEqual(s.value, ref)

#     def test_invalid(self):
#         s = KinematicsSpec(name = 'name')
#         s.observable_name = 'B->Dlnu::dBR/dq2;l=mu'
#         with self.assertRaises(ValueError):
#             s.validate_set({'q2': 1.0, 'q3': 5.0})

#     def test_wrong_argument_no_dict(self):
#         s = KinematicsSpec(name = 'kinematics')
#         try:
#             not_a_dict = 1.0
#             s.observable_name = 'B->Dlnu::dBR/dq2;l=mu'
#             s.validate_set(not_a_dict)
#         except Exception as e:
#             err_str = r"'kinematics' must be a dictionary of form {'kinematic_name': value}"
#             self.assertEqual(str(e), err_str)


# class ParametersSpec(BaseSpec):

#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)
#         self.observable_name = None

#     def validate_set(self, value):
#         "``dict`` like ``{parameter: value}`` -- ``parameter`` must be a valid name, ``value`` must be a float"

#         if self.observable_name == None:
#             raise ImplementationError("Attribute 'observable_name' not set")

#         if not isinstance(value, dict):
#             raise ValueError(f"'{self.name}' must be a dictionary of form {{'parameter_name': value}}")

#         invalid_pars = []
#         for par in value.keys():
#             try:
#                 qn = eos.QualifiedName(par)
#             except:
#                 invalid_pars.append(par)

#         if not invalid_pars == []:
#             raise ValueError(f"'{self.name}' for observable '{self.observable_name}'"
#                              f" contain invalid names: {invalid_pars}")
#         else:
#             self._value = value


# class ParametersSpecTest(unittest.TestCase):

#     def test_valid(self):
#         s = ParametersSpec(name='parameters')
#         s.observable_name = 'some_observable'
#         ref = {'mass::tau': 1.0}
#         s.validate_set(ref)
#         self.assertEqual(s.value, ref)

#     def test_invalid(self):
#         s = ParametersSpec(name='parameters')
#         s.observable_name = 'some_observable'
#         ref = {'mass::tau': 1.0, 'nada': 0.5}
#         ref_err = r"'parameters' for observable 'some_observable' contain invalid names: ['nada']"
#         try:
#             s.validate_set(ref)
#         except Exception as e:
#             self.assertEqual(str(e), ref_err)

#     def test_wrong_argument_no_dict(self):
#         s = ParametersSpec(name='parameters')
#         s.observable_name = 'B->Dlnu::dBR/dq2;l=mu'
#         try:
#             not_a_dict = 1.0
#             s.validate_set(not_a_dict)
#         except Exception as e:
#             err_str = r"'parameters' must be a dictionary of form {'parameter_name': value}"
#             self.assertEqual(str(e), err_str)


# class RangeSpec(BaseSpec):

#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)

#     def validate_set(self, value):
#         "(*list* or *tuple* of two *float* like [min, max]) -- Range of the x-axis"
#         try:
#             if len(value) == 2:
#                 self._value = value
#             else:
#                 raise ValueError()
#         except:
#             raise ValueError(f"{self.name} must be a list or tuple of float")


# class RangeSpecTests(unittest.TestCase):

#     def test_valid1(self):
#         "Pass a list with two items"
#         s = RangeSpec(name = 'range')
#         ref = [1.0, 2.0]
#         s.validate_set(ref) # must not raise error
#         self.assertEqual(s.value, ref)

#     def test_valid2(self):
#         "Pass a tuple with two items"
#         s = RangeSpec(name = 'range')
#         ref = (1.0, 2.0)
#         s.validate_set(ref) # must not raise error
#         self.assertEqual(s.value, ref)

#     def test_invalid(self):
#         "Pass a value instead of list"
#         s = RangeSpec(name = 'range')
#         with self.assertRaises(ValueError):
#             s.validate_set(1.0)

#     def test_invalid2(self):
#         "Pass a list with three items"
#         s = RangeSpec(name = 'range')
#         with self.assertRaises(ValueError):
#             s.validate_set([1.0, 2.0, 3.0])


# class ParameterFileSpec(BaseSpec):

#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)

#     def validate_set(self, value):
#         "``str`` -- the path of a hdf5 file containing parameters"
#         if not isinstance(value, str):
#             raise ValueError(f"{self.name} must be a string")

# class ParameterFileSpecTests(unittest.TestCase):

#     def test_valid(self):
#         s = ParameterFileSpec(name = 'parameters-from-file')
#         s.validate_set('some/path.file')

#     def test_valid(self):
#         s = ParameterFileSpec(name = 'parameters-from-file')
#         with self.assertRaises(ValueError):
#             s.validate_set(3.1415)


if __name__ == '__main__':
    unittest.main(verbosity=7)

    # template = (
    #     ObservableSpec(name='observable'),
    #     VariableSpec(name='variable'),
    #     KinematicsSpec(name='kinematics', optional=True),
    #     RangeSpec(name='range')
    # )
    # s = Specs(template)
    # for name in s.names:
    #     opt = ''
    #     if s[name].optional:
    #         opt = ' (optional) '
    #     print(f"{name}{opt}, {s[name].validate_set.__doc__}")
