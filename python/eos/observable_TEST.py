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



class FilterTests(unittest.TestCase):

    _QN = eos.QualifiedName('B->pi::f_+(q2)@BCL2008')

    def test_no_filter(self):
        "Without a filter every entry is accepted."

        self.assertTrue(eos.Observables().filter_entry(self._QN))

    def test_prefix(self):
        "The prefix filter matches a substring of the prefix part."

        self.assertTrue (eos.Observables(prefix='B->pi').filter_entry(self._QN))
        self.assertFalse(eos.Observables(prefix='B->D').filter_entry(self._QN))

    def test_name(self):
        "The name filter matches a substring of the name part."

        self.assertTrue (eos.Observables(name='f_+').filter_entry(self._QN))
        self.assertFalse(eos.Observables(name='f_0').filter_entry(self._QN))

    def test_suffix(self):
        "The suffix filter matches a substring of the suffix part."

        self.assertTrue (eos.Observables(suffix='BCL2008').filter_entry(self._QN))
        self.assertFalse(eos.Observables(suffix='BSZ2015').filter_entry(self._QN))

    def test_combined(self):
        "All provided filters must match."

        self.assertTrue (eos.Observables(prefix='B->pi', name='f_+', suffix='BCL2008').filter_entry(self._QN))
        self.assertFalse(eos.Observables(prefix='B->pi', name='f_+', suffix='BSZ2015').filter_entry(self._QN))


class ReprHTMLTests(unittest.TestCase):

    def test_complete_list(self):
        "The complete list of observables renders with one row per observable."

        result = eos.Observables()._repr_html_()

        self.assertIn('<th rowspan="2">qualified name</th>', result)
        self.assertIn('<tt>B->Dlnu::dBR/dq2</tt>', result)

    def test_filtered_list(self):
        "A filter removes both the entries it rejects and the groups left empty."

        result = eos.Observables(prefix='B->Dlnu', name='dBR/dq2')._repr_html_()

        self.assertIn   ('<tt>B->Dlnu::dBR/dq2</tt>',  result)
        self.assertNotIn('<tt>B->pilnu::dBR/dq2</tt>', result)

    def test_empty_list(self):
        "A filter that rejects every observable yields a table without rows."

        result = eos.Observables(prefix='nosuchprefix')._repr_html_()

        self.assertNotIn('<tt>', result)

    def test_kinematic_variables(self):
        "The kinematic variables of an entry are listed, or replaced by a dash."

        entry = eos.Observables()['B->Dlnu::dBR/dq2;l=e,q=d']
        kvs   = [str(kv) for kv in entry.kinematic_variables()]

        self.assertGreater(len(kvs), 0)

        result = eos.Observables(prefix='B->Dlnu', name='dBR/dq2')._repr_html_()
        self.assertIn('<br>'.join([f'<tt>{kv}</tt>' for kv in kvs]), result)
        self.assertIn('&mdash;', eos.Observables(prefix='0->Kpi', name='Delta_CT')._repr_html_())

    def test_options(self):
        "The options of an entry are listed with their allowed and default values."

        entry  = eos.Observables()['B->Dlnu::dBR/dq2;l=e,q=d']
        result = eos.Observables(prefix='B->Dlnu', name='dBR/dq2')._repr_html_()

        for option in entry.options():
            self.assertIn(f'<tt>{option.key}</tt>', result)
            for allowed_value in option.allowed_values:
                self.assertIn(f'<tt>{allowed_value}</tt>', result)

    def test_showall(self):
        "Observables without a symbol are shown only when 'showall' is requested."

        entry = eos.Observables()['B->D::f_+[s^1/s^0](q2)']
        self.assertEqual(entry.latex(), '')

        self.assertNotIn('<tt>B->D::f_+[s^1/s^0](q2)</tt>', eos.Observables(prefix='B->D')._repr_html_())
        self.assertIn   ('<tt>B->D::f_+[s^1/s^0](q2)</tt>', eos.Observables(prefix='B->D', showall=True)._repr_html_())


if __name__ == '__main__':
    unittest.main(verbosity=5)
