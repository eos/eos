import os
# The report task draws a corner figure; select a head-less backend before eos (and thus
# matplotlib) is imported.
os.environ.setdefault('MPLBACKEND', 'Agg')

import shutil
import tempfile
import unittest
import eos
import yaml
from pathlib import Path


class ReportTaskTests(unittest.TestCase):

    _fixture = Path(__file__).parent.parent / "reporting_TEST.d"

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
        # Parameter names are typeset as math; emitting the raw LaTeX as text breaks the PDF.
        self.assertIn('- `CKM::abs(V_ub)` ($|V_{ub}|$)', rendered)
        self.assertIn('- `B->pi::f_+(0)@BCL2008` ($f_+^{B\\to \\pi}(0)$)', rendered)
        # A corner figure is produced for 'CKM' and referenced from the report.
        self.assertIn('figures/corner-CKM.pdf', rendered)
        self.assertTrue(os.path.isfile(os.path.join(base, 'figures', 'corner-CKM.pdf')))

    def test_report_task_load_yaml(self):
        "Render a template that reads a YAML file below the base directory."
        workdir = tempfile.mkdtemp(prefix='eos-report-task-')
        self.addCleanup(shutil.rmtree, workdir, ignore_errors=True)
        base = os.path.join(workdir, 'base')
        shutil.copytree(self._fixture, base)
        os.makedirs(os.path.join(base, 'results'))
        with open(os.path.join(base, 'results', 'accuracy.yaml'), 'w') as f:
            yaml.safe_dump({ 'max_relative_error': 1.5e-15 }, f)
        with open(os.path.join(base, 'accuracy.md.jinja'), 'w') as f:
            f.write("error: {{ '%.1e' | format(load_yaml('results/accuracy.yaml').max_relative_error) }}\n")

        cwd = os.getcwd()
        self.addCleanup(os.chdir, cwd)
        os.chdir(base)

        eos.tasks.report('analysis.yaml', 'accuracy.md.jinja', base_directory='.', generate_pdf=False)

        with open(os.path.join(base, 'accuracy.md')) as f:
            self.assertEqual(f.read().strip(), 'error: 1.5e-15')


if __name__ == '__main__':
    unittest.main(verbosity=5)
