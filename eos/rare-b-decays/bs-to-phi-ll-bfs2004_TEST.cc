/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021 Méril Reboud
 * Copyright (c) 2025 Danny van Dyk
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <test/test.hh>
#include <eos/maths/complex.hh>
#include <eos/observable.hh>
#include <eos/rare-b-decays/bs-to-phi-ll.hh>

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

using namespace test;
using namespace eos;

class BsToPhiDileptonBFS2004NaiveTest :
    public TestCase
{
    public:
    BsToPhiDileptonBFS2004NaiveTest() :
        TestCase("bs_to_phi_dilepton_BFS2004_naive_test")
    {
    }

    virtual void run() const
    {
        {
            Parameters p = Parameters::Defaults();

            Options oo
            {
                { "model"_ok, "WET" },
                { "scan-mode"_ok, "cartesian" },
                { "tag"_ok, "BFS2004" },
                { "qcdf-integrals"_ok, "mixed" },
                { "form-factors"_ok, "BSZ2015" },
                { "l"_ok, "mu" },
                { "q"_ok, "s" }
            };

            static const double eps = 1e-3;
            static const double q2 = 6.0;

            BsToPhiDilepton d(p, oo);
            auto amps = d.amplitudes(q2);

            TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_left,  complex<double>(-1.23979e-10,3.78483e-15)  , eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_right, complex<double>(7.05843e-12,3.78483e-15)   , eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_left,  complex<double>(-5.56392e-11,1.80598e-12)  , eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_right, complex<double>(2.50889e-11,1.80598e-12)   , eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_left,  complex<double>(4.85643e-11,-1.66461e-12)  , eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_right, complex<double>(-2.31595e-11,-1.66461e-12) , eps);
       }
    }
} bs_to_phi_dilepton_BFS2004_naive_test;
