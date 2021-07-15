import unittest
import eos

class StaticMethodTests(unittest.TestCase):

    def test_assert_valid_name_0(self):
        "Value of 'name' is valid"

        eos.Parameters._assert_valid_name('mass::mu')


    def test_assert_valid_name_1(self):
        "Value of 'name' is invalid (for small-scale physics)"

        with self.assertRaisesRegex(ValueError, "Parameter name 'mass::your_mum' is not known"):
            eos.Parameters._assert_valid_name('mass::your_mum')


if __name__ == '__main__':
    unittest.main(verbosity=5)

