import unittest
import eos

class StaticMethodTests(unittest.TestCase):

    def test_assert_valid_name_0(self):
        "name is valid"

        valid_name = 'B->Dlnu::dBR/dq2;l=e,q=d'
        eos.Observables.assert_valid_name(valid_name)


    def test_assert_valid_name_1(self):
        "name is invalid"

        invalid_name = 'prefix::Philipp'
        with self.assertRaisesRegex(ValueError, "Observable with name 'prefix::Philipp' is not known"):
            eos.Observables.assert_valid_name(invalid_name)


    def test_assert_valid_name_and_kinematic_0(self):
        "name and kinamtic name are valid"

        valid_name = 'B->Dlnu::dBR/dq2;l=e,q=d'
        valid_kinematic = 'q2'
        eos.Observables.assert_valid_name_and_kinematic(valid_name, valid_kinematic)


    def test_assert_valid_name_and_kinematic_1(self):
        "name and kinamtic name are valid"

        valid_name = 'B->Dlnu::dBR/dq2;l=e,q=d'
        invalid_kinematic = 'qToTheTwo'
        with self.assertRaisesRegex(ValueError,
            "Kinematic variable 'qToTheTwo' does not match known ones: \['q2'\]"):
            eos.Observables.assert_valid_name_and_kinematic(valid_name, invalid_kinematic)


if __name__ == '__main__':
    unittest.main(verbosity=5)

