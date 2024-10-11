/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Matthew Kirk
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
#include <eos/form-factors/parametric-kkrvd2024.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class ParametricKKRvD2024Test :
    public TestCase
{
    public:
        ParametricKKRvD2024Test() :
            TestCase("parametric_KKRvD2024_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-7;

            // t0 = -1
            {

                Parameters p = Parameters::Defaults();
                p["mass::pi^+"]                    =  0.13957;
                p["0->pipi::t_0@KKRvD2024"]         = -1.0;
                p["0->pipi::b_(+,1)^2@KKRvD2024"]   = -0.0182238;
                p["0->pipi::b_(+,1)^3@KKRvD2024"]   = -0.0225337;
                p["0->pipi::b_(+,1)^4@KKRvD2024"]   = -0.00731115;
                p["0->pipi::b_(+,1)^5@KKRvD2024"]   = -0.00483614;
                p["0->pipi::b_(+,1)^6@KKRvD2024"]   =  0;
                p["0->pipi::b_(+,1)^7@KKRvD2024"]   =  0;
                p["0->pipi::b_(+,1)^8@KKRvD2024"]   =  0;
                p["0->pipi::b_(+,1)^9@KKRvD2024"]   =  0;
                p["0->pipi::M_(+,1)@KKRvD2024"]     =  0.760895;
                p["0->pipi::Gamma_(+,1)@KKRvD2024"] =  0.146155;

                /* 0->PP factory */
                {
                    std::shared_ptr<FormFactors<VacuumToPP>> ff = FormFactorFactory<VacuumToPP>::create("0->pipi::KKRvD2024", p, Options{ });

                    TEST_CHECK(nullptr != ff);
                }

                /* f_+ at timelike q2 > 0.0 */
                {
                    KKRvD2024FormFactors<VacuumToPiPi> ff(p, Options{ });

                    const auto chi = 0.00683918; // GeV^-2, at Q^2 = 1 GeV^2 using [BL:1998A] Sec VI.A

                    TEST_CHECK_NEARLY_EQUAL(real(ff.z( 0.0)),                     -0.576215878,   eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.z( 0.0)),                      0.0,           eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.z(+0.1)),                     -0.959852981,   eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.z(+0.1)),                     -0.280503573,   eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.z(+0.5)),                     -0.437225519,   eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.z(+0.5)),                     -0.899351903,   eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitilde_p(ff.z( 0.0), chi)),  0.1036646321,  eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitilde_p(ff.z( 0.0), chi)),  0.0,           eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitilde_p(ff.z(+0.1), chi)),  0.0774207417,  eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitilde_p(ff.z(+0.1), chi)), -0.0086141281,  eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.phitilde_p(ff.z(+0.5), chi)),  0.0526299641,  eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.phitilde_p(ff.z(+0.5), chi)), -0.0558825275,  eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.b_0(),                              0.0901897577,  eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.b_1(),                             -0.0451794361,  eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p( 0.0)),                    1.0,           eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p( 0.0)),                    0.0,           eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(+0.1)),                    1.263408016,   eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(+0.1)),                    0.014150225,   eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(+0.5)),                    3.60167333,    eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(+0.5)),                    3.76175784,    eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.dFdq2_q2eq0(),                      1.93155736,    eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.r_pi_squared(),                     0.451265116,   eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.re_residue_rho(),                  -0.686247061,   eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.im_residue_rho(),                   0.272035152,   eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.re_residue_rho_q2(),               -0.708332619,   eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.im_residue_rho_q2(),                0.135266983,   eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.saturation(),                       0.444609294,   eps);
                }
            }
        }
} parametric_KKRvD2024_test;
