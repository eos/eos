import unittest
import eos

class StaticMethodTests(unittest.TestCase):

    def test_get_obs_entry_0(self):
        "name is valid"

        valid_name = 'B->Dlnu::dBR/dq2;l=e,q=d'
        eos.Observables._get_obs_entry(valid_name)


    def test_get_obs_entry_1(self):
        "name is invalid"

        invalid_name = 'prefix::Philipp'
        with self.assertRaisesRegex(ValueError, "Observable with name 'prefix::Philipp' is not known"):
            eos.Observables._get_obs_entry(invalid_name)


if __name__ == '__main__':
    unittest.main(verbosity=5)

