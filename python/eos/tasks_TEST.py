import unittest
import eos
from pathlib import Path

class ClassMethodTests(unittest.TestCase):

    _analysis_file = str(Path(__file__).parent / "analysis_file_TEST.d/valid-analysis-file.yaml")

    def test_validate_task(self):
        eos.tasks.validate(self._analysis_file)

    def test_list_steps_task(self):
        steps = eos.tasks.list_steps(self._analysis_file)
        self.assertEqual(steps, {'CKM-all,WET-all.sample', 'CKM-all.corner-plot', 'WET-all.mode,corner-plot'})


if __name__ == '__main__':
    unittest.main(verbosity=5)
