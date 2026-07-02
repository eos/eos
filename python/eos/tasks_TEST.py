import os
# The report task draws a corner figure; select a head-less backend before eos (and thus
# matplotlib) is imported.
os.environ.setdefault('MPLBACKEND', 'Agg')

import shutil
import tempfile
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


class ReportTaskTests(unittest.TestCase):

    _fixture = Path(__file__).parent / "reporting_TEST.d"

    def test_report_task(self):
        "Render a report from the reporting_TEST.d fixtures (without generating a PDF)."
        # Copy the fixture into a temporary directory and render there, so that neither the
        # intermediate file, the drawn figures, nor the task's log pollute the source tree.
        workdir = tempfile.mkdtemp(prefix='eos-report-task-')
        self.addCleanup(shutil.rmtree, workdir, ignore_errors=True)
        base = os.path.join(workdir, 'base')
        shutil.copytree(self._fixture, base)

        cwd = os.getcwd()
        self.addCleanup(os.chdir, cwd)
        os.chdir(base)

        eos.tasks.report('analysis.yaml', 'report.md.jinja', base_directory='.', generate_pdf=False)

        with open(os.path.join(base, 'report.md')) as f:
            rendered = f.read()

        # Both posteriors are discovered from the fixture's data/ tree.
        self.assertIn('## Posterior `CKM`', rendered)
        self.assertIn('## Posterior `FF`', rendered)
        # 'CKM' has a recorded mode, hence a goodness-of-fit table with a local p-value and a total row.
        self.assertIn('`B->pilnu::BR`', rendered)
        self.assertIn('**total**', rendered)
        # 'FF' has no mode, hence the goodness-of-fit fallback.
        self.assertIn('No goodness-of-fit information has been recorded', rendered)
        # A corner figure is produced for 'CKM' and referenced from the report.
        self.assertIn('figures/corner-CKM.pdf', rendered)
        self.assertTrue(os.path.isfile(os.path.join(base, 'figures', 'corner-CKM.pdf')))


if __name__ == '__main__':
    unittest.main(verbosity=5)
