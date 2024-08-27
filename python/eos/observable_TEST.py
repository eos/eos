import unittest
import eos

class ClassOperatorTests(unittest.TestCase):

    def test_get_obs_entry_0(self):
        "name is valid"

        valid_name = 'B->Dlnu::dBR/dq2;l=e,q=d'
        eos.Observables()[valid_name]


    def test_get_obs_entry_1(self):
        "name is invalid"

        invalid_name = 'prefix::Philipp'
        with self.assertRaisesRegex(RuntimeError, "Unknown Observable Error: 'prefix::Philipp' not known"):
            eos.Observables()[invalid_name]


if __name__ == '__main__':
    unittest.main(verbosity=5)
