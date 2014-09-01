#! /usr/bin/env python

## Copyright (c) 2014 Stephan Jahn
##
## This file is part of the EOS project. EOS is free software;
## you can redistribute it and/or modify it under the terms of the GNU General
## Public License version 2, as published by the Free Software Foundation.
##
## EOS is distributed in the hope that it will be useful, but WITHOUT ANY
## WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, write to the Free Software Foundation, Inc., 59 Temple
## Place, Suite 330, Boston, MA  02111-1307  USA

import eos
from eos import EOSError
import numpy as np
import unittest

class PosteriorTest(unittest.TestCase):
    def test_posterior_construction(self):
        # constraints to be considered
        constraints = ["B->X_sll::BR[1.0,6.0]@Belle-2005A",
                       "B^0->K^*0gamma::S_K@Belle-2006",
                       "B->X_sgamma::E_1[1.8]+E_2[1.8]@Belle-2008"]

        # Parameters to be scanned
        priors = []
        priors.append(  eos.LogPrior.Flat("QCD::alpha_s(MZ)", range_min=-15, range_max=15)  )
        priors.append(  eos.LogPrior.Gauss("CKM::A", range_min=1.0, range_max=2.834,
                                                lower=1.04-0.1, central=1.04, upper=1.04+0.01 ) )
        priors.append(  eos.LogPrior.Gauss("QCD::mu_t", range_min=0.774, range_max=0.834,
                                                lower=0.804-0.01, central=0.804, upper=0.804+0.01 ) )
        priors.append(  eos.LogPrior.LogGamma("QCD::mu_b", range_min=1, range_max=2,
                                                lower=1.5-0.2, central=1.5, upper=1.5+0.1 ) )

        options = {"scan-mode": "cartesian", "model": "WilsonScan", "form-factors": "BZ2004"}

        # should be a valid call
        ana = eos.Analysis(constraints, priors, options)

    def test_constraint_options(self):
        # add some local option
        constraints = ["B->X_sll::BR[1.0,6.0]@Belle-2005A",
                       ("B^0->K^*0gamma::S_K@Belle-2006", {"form-factors": "KMPW2010"}),
                       "B->X_sgamma::E_1[1.8]+E_2[1.8]@Belle-2008",
                       ("B->K^*ll::BRavg@LowRecoil", (2e-7, 3e-7, 4e-7), 1, {"s_min": 15., "s_max": 16.}, {"form-factors": "KMPW2010", "q": "d", "l": "mu"}),
                       eos.Constraint("B^0_s->mu^+mu^-::BR@CMS-2013B")]

        # Parameters to be scanned
        priors = []
        priors.append(eos.LogPrior.Flat("QCD::alpha_s(MZ)", range_min=-15, range_max=15))
        priors.append(eos.LogPrior.Gauss("CKM::A", range_min=1.0, range_max=2.834,
                                         lower=1.04-0.1, central=1.04, upper=1.04+0.01))

        global_options = {"scan-mode": "cartesian", "model": "WilsonScan", "form-factors": "BZ2004"}

        # should be a valid call
        ana = eos.Analysis(constraints, priors, global_options)

        self.assertEqual(ana.constraints[0].options["form-factors"], "BZ2004")
        self.assertEqual(ana.constraints[1].options["form-factors"], "KMPW2010")
        self.assertEqual(ana.constraints[2].options["form-factors"], "BZ2004")

    def test_error_messages(self):
        params_OK       = [eos.LogPrior.Flat("QCD::alpha_s(MZ)", range_min=-15, range_max=15),
                           eos.LogPrior.Gauss("CKM::A", range_min=1.0, range_max=2.834,
                                               lower=1.04-0.01, central=1.04, upper=1.04+0.01 )]
        one_param_twice = [eos.LogPrior.Flat("QCD::alpha_s(MZ)", range_min=-15, range_max=15),
                           eos.LogPrior.Gauss("QCD::mu_t", range_min=0.774, range_max=0.834,
                                               lower=0.804-0.01, central=0.804, upper=0.804+0.01 ),
                           eos.LogPrior.LogGamma("QCD::mu_t", range_min=0, range_max=2,
                                                  lower=0.5-0.2, central=0.5, upper=0.5+0.1 )]
        unknown_param   = [eos.LogPrior.Flat("QCD::alpha_s(MZ)", range_min=-15, range_max=15),
                           eos.LogPrior.Gauss("QCD::mu_t", range_min=0.774, range_max=0.834,
                                               lower=0.804-0.01, central=0.804, upper=0.804+0.01 ),
                           eos.LogPrior.LogGamma("i_do_not_exist", range_min=0, range_max=2,
                                                  lower=0.5-0.2, central=0.5, upper=0.5+0.1 )]
        invalid_range   = [eos.LogPrior.Flat("QCD::alpha_s(MZ)", range_min=-15, range_max=15),
                           eos.LogPrior.Gauss("QCD::mu_t", range_min=0.774, range_max=0.834,
                                               lower=0.804-0.01, central=0.804, upper=0.804+0.01 ),
                           eos.LogPrior.LogGamma("CKM::A", range_min=9, range_max=2,
                                                  lower=0.5-0.2, central=0.5, upper=0.5+0.1 )]
        unknown_type    = ["CKM:A"]
        constraints_OK  = ["B->X_sll::BR[1.0,6.0]@Belle-2005A"]
        unknown_constr  = ["B->K@never"]

        self.assertRaisesRegexp(EOSError, "parameter.*QCD::mu_t.*twice", eos.Analysis, constraints_OK, one_param_twice)
        self.assertRaisesRegexp(EOSError, "Unknown.*i_do_not_exist", eos.Analysis, constraints_OK, unknown_param)
        self.assertRaisesRegexp(EOSError, "Range Error.*CKM::A", eos.Analysis, constraints_OK, invalid_range)
        self.assertRaisesRegexp(EOSError, "Unknown prior type*", eos.Analysis, constraints_OK, unknown_type)
        self.assertRaisesRegexp(EOSError, "'B->K@never'.*unknown", eos.Analysis, unknown_constr, params_OK)

    def test_call(self):
        constraints = ["B^0_s->mu^+mu^-::BR@CMS-2013B"]
        priors = [eos.LogPrior.Flat("Re{c10}", range_min=-15, range_max=15)]
        options = {"scan-mode": "cartesian", "model": "WilsonScan", "form-factors": "BZ2004"}

        ana = eos.Analysis(constraints, priors, options)

        # store result in array, make sure in-place modification
        # is propagated to eos
        x = np.array([2.3])
        result = ana(x)
        self.assertAlmostEqual(result, 13.614259950362, places=12)

        # better fit near SM value of -4.17
        x[0] = -4.1
        self.assertGreater(ana(x), result)

        self.assertRaisesRegexp(AssertionError, 'must.*\(1,\).*got.*\(0,\)',  ana, [])
        self.assertRaisesRegexp(AssertionError, 'must.*\(1,\).*got.*\(2,\)',  ana, [1,5.4])


        priors = [eos.LogPrior.Gauss("Re{c10}", range_min=-15, range_max=15,
                                      lower=0, upper=5, central=2)]

        ana = eos.Analysis(constraints, priors, options)

        self.assertAlmostEqual(ana([2.3]), 14.992913853602, places=12)


        priors = [eos.LogPrior.LogGamma("Re{c10}", range_min=-15, range_max=15,
                                         lower=0, upper=5, central=2)]

        ana = eos.Analysis(constraints, priors, options)

        self.assertAlmostEqual(ana([2.3]), 15.210677762849, places=12)

class PriorTest(unittest.TestCase):
    def test_repr(self):
        flat_prior_description     = eos.LogPrior.Flat    ("CKM::A", range_min=-15, range_max=15)
        gauss_prior_description    = eos.LogPrior.Gauss   ("QCD::mu_t", range_min=1.0, range_max=2.834,
                                                           lower=1.04-0.01, central=1.04, upper=1.04+0.01 )
        loggamma_prior_description = eos.LogPrior.LogGamma("Re{c10}", range_min=0, range_max=2,
                                                           lower=0.5-0.2, central=0.5, upper=0.5+0.1 )

        self.assertEqual(repr(flat_prior_description),     str(flat_prior_description))
        self.assertEqual(repr(gauss_prior_description),    str(gauss_prior_description))
        self.assertEqual(repr(loggamma_prior_description), str(loggamma_prior_description))

        self.assertEqual(repr(flat_prior_description),     'Prior "Flat" for parameter "CKM::A"')
        self.assertEqual(repr(gauss_prior_description),    'Prior "Gauss" for parameter "QCD::mu_t"')
        self.assertEqual(repr(loggamma_prior_description), 'Prior "LogGamma" for parameter "Re{c10}"')


        priors = [flat_prior_description, gauss_prior_description, loggamma_prior_description]
        constraints = ["B^0_s->mu^+mu^-::BR@CMS-2013B"]
        options = {"scan-mode": "cartesian"}
        ana = eos.Analysis(constraints, priors, options)

        for x,y in zip(ana.constraints, constraints):
            self.assertEqual(x.name, y)
        self.assertEqual(ana.priors, priors)
        self.assertEqual(ana.global_options, options)

        # deep copy, not just reference?
        self.assertNotEqual(id(ana.priors), id(priors))
        self.assertNotEqual(id(ana.constraints), id(constraints))
        self.assertNotEqual(id(ana.global_options), id(options))

    def test_nuisance(self):
        gauss_prior_description = eos.LogPrior.Gauss("QCD::mu_t", range_min=1.0, range_max=2.834,
                                                           lower=1.04-0.01, central=1.04, upper=1.04+0.01,
                                                           nuisance=True)
        flat_prior_description = eos.LogPrior.Flat("Re{C9}", range_min=-15, range_max=15)
        self.assert_(gauss_prior_description.nuisance)
        self.assertFalse(flat_prior_description.nuisance)

class SetNumberOfIntegrationPointsTest(unittest.TestCase):
    def test_set_and_get_n(self):
        # constraints to be considered
        constraints = ["B^+->K^+mu^+mu^-::BR[1.00,6.00]@BaBar-2012"]

        # Parametereters to be scanned
        priors = [eos.LogPrior.Flat("mass::c"       , range_min=1, range_max=15),
                  eos.LogPrior.Flat("mass::b(MSbar)", range_min=1, range_max=15)]

        options = {"scan-mode": "cartesian", "model": "WilsonScan", "form-factors": "KMPW2010"}

        ana = eos.Analysis(constraints, priors, options)

        intial_integrate_n = eos.get_integrate_n()
        self.assertEqual(intial_integrate_n, 64)

        result1 = ana([1.,2.])

        final_integrate_n = 16
        eos.set_integrate_n(final_integrate_n)
        self.assertEqual( final_integrate_n, eos.get_integrate_n() )

        result2 = ana([1.,2.])

        self.assertNotEqual(result1, result2)

    def test_python_errors(self):
        with self.assertRaises(EOSError):
            eos.set_integrate_n(5)

        with self.assertRaises(EOSError):
            eos.set_integrate_n(4)

        with self.assertRaises(EOSError):
            eos.set_integrate_n(84)

if __name__ == '__main__':
    unittest.main()
