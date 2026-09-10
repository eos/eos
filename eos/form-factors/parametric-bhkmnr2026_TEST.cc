/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Fatemeh Nouri
 * Copyright (c) 2026 Méril Reboud
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

#include <eos/form-factors/parametric-bhkmnr2026.hh>

#include <test/test.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class ParametricBHKMNR2026Test : public TestCase
{
    public:
        ParametricBHKMNR2026Test() :
            TestCase("parametric_BHKMNR2026_test")
        {
        }

        virtual void
        run() const
        {
            static const double eps = 1e-7;

            {
                Parameters p                           = Parameters::Defaults();
                p["mass::pi^+"]                        = 0.13957;
                p["0->pipi::s_0@BHKMNR2026"]           = 0.0;
                p["0->pipi::s_in@BHKMNR2026"]          = 1.0;
                p["0->pipi::a_(+,1)^4@BHKMNR2026"]     = 0.20;
                p["0->pipi::a_(+,1)^5@BHKMNR2026"]     = 0.12;
                p["0->pipi::a_(+,1)^6@BHKMNR2026"]     = 0.07;
                p["0->pipi::a_(+,1)^7@BHKMNR2026"]     = 0.02;
                p["0->pipi::a_(+,1)^8@BHKMNR2026"]     = 0;
                p["0->pipi::a_(+,1)^9@BHKMNR2026"]     = 0;
                p["0->pipi::a_(+,1)^10@BHKMNR2026"]    = 0;
                p["0->pipi::a_(+,1)^11@BHKMNR2026"]    = 0;
                p["0->pipi::a_(+,1)^12@BHKMNR2026"]    = 0;
                p["0->pipi::M_(+,1,0)@BHKMNR2026"]     = 0.760895;
                p["0->pipi::Gamma_(+,1,0)@BHKMNR2026"] = 0.146155;

                /* 0->PP factory */
                {
                    std::shared_ptr<FormFactors<VacuumToPP>> ff = FormFactorFactory<VacuumToPP>::create("0->pipi::BHKMNR2026", p, Options{});

                    TEST_CHECK(nullptr != ff);
                }


                {
                    Options o{
                        { "n-resonances-I1"_ok, "1"_ov },
                        {               "I"_ok, "1"_ov }
                    };
                    BHKMNR2026FormFactors<VacuumToPiPi> ff(p, o);


                    TEST_CHECK_NEARLY_EQUAL(real(ff.psi(0.0)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.psi(0.0)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.psi(0.077919139600000)), -0.24972020, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.psi(0.077919139600000)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.psi(0.5)), -0.02968043, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.psi(0.5)), 0.55209340, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.psi(1.0)), 0.12710102, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.psi(1.0)), 0.99188978, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.psi(complex<double>(0.5, 0.5))), 0.22651223, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.psi(complex<double>(0.5, 0.5))), -0.44265917, eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.psi21(0.0)), -1.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.psi21(0.0)), 0.00000000, eps);


                    TEST_CHECK_NEARLY_EQUAL(real(ff.P(complex<double>(0, 0))), 2.53121472, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.P(complex<double>(0, 0))), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.P(complex<double>(0.5, 0))), 0.52719604, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.P(complex<double>(0.5, 0))), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.P(complex<double>(1, 0))), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.P(complex<double>(1, 0))), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.P(complex<double>(0.5, 0.5))), -0.57992579, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.P(complex<double>(0.5, 0.5))), -0.91396391, eps);


                    TEST_CHECK_NEARLY_EQUAL(real(ff.dPdpsi(0.0)), -3.38006038, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dPdpsi(0.0)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dPdpsi(0.5)), -2.59667560, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dPdpsi(0.5)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dPdpsi(complex<double>(0.5, 0.5))), 0.29304014, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dPdpsi(complex<double>(0.5, 0.5))), 4.22501980, eps);


                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms_I1(0, 0.0)), -3.38006038, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms_I1(0, 0.0)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms_I1(1, 0.0)), 2.53121472, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms_I1(1, 0.0)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms_I1(0, 0.5)), -2.59667560, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms_I1(0, 0.5)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms_I1(1, 0.5)), -0.77114178, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms_I1(1, 0.5)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms_I1(1, complex<double>(0.5, 0.5))), -2.54591562, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms_I1(1, complex<double>(0.5, 0.5))), 1.34506606, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms_I1(2, complex<double>(0.5, 0.5))), -1.77847178, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms_I1(2, complex<double>(0.5, 0.5))), -1.34736963, eps);


                    const auto constrained_a = ff.constrained_a_fp_I1();
                    TEST_CHECK_NEARLY_EQUAL(constrained_a[0], 0.39506723, eps);
                    TEST_CHECK_NEARLY_EQUAL(constrained_a[1], -0.08512273, eps);
                    TEST_CHECK_NEARLY_EQUAL(constrained_a[2], 0.27898849, eps);
                    TEST_CHECK_NEARLY_EQUAL(constrained_a[3], -0.13769399, eps);


                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(0.0)), 1.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(0.0)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(0.5)), 2.89045313, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(0.5)), 4.26305734, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(complex<double>(0.5, 0.5))), 0.14423426, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(complex<double>(0.5, 0.5))), 0.94360544, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(1.0)), -1.45888696, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(1.0)), 0.04846081, eps);


                    TEST_CHECK_NEARLY_EQUAL(real(ff.partial_wave(0.0)), -0.21320718, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.partial_wave(0.0)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.partial_wave(0.077919139600000)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.partial_wave(0.077919139600000)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.partial_wave(complex<double>(0.5, 0.5))), 0.02966105, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.partial_wave(complex<double>(0.5, 0.5))), 0.29786040, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.partial_wave(1.0)), -0.03455458, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.partial_wave(1.0)), 0.00114782, eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.scattering_length_parameters()[0]), -5.35424361, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.scattering_length_parameters()[0]), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.scattering_length_parameters()[1]), -8.91801603, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.scattering_length_parameters()[1]), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.scattering_length_parameters()[2]), 81.96437920, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.scattering_length_parameters()[2]), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.scattering_length_parameters()[3]), 997.64746761, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.scattering_length_parameters()[3]), 0.00000000, eps);

                    // needed for charged pion radius
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_11(0.0)), -1.55081502, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_11(0.0)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_11(0.077919139600000)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_11(0.077919139600000)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_11(complex<double>(0.5, 0.5))), -2.20224209, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_11(complex<double>(0.5, 0.5))), -3.43995730, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_11(1.0)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_11(1.0)), 0.00000000, eps);


                    TEST_CHECK_RELATIVE_ERROR(ff.saturation(), 0.398058745, eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p_of_psi(complex<double>(0.5, 0.5))), -0.15313783, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p_of_psi(complex<double>(0.5, 0.5))), -0.31712693, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p_of_psi(complex<double>(0.7, 0.3))), 0.00458081, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p_of_psi(complex<double>(0.7, 0.3))), -0.13037101, eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.abs2_f_p_of_psi(0.5, 0.5), 0.12402069, eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.abs2_f_p_of_psi(0.7, 0.3), 0.01701758, eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.abs2_f_p_of_psi(0.0, 0.0), 1.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.abs2_f_p_of_psi(-1.0, 0.0), 0.00000000, eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.arg_f_p_of_psi(0.5, 0.5), -2.02066352, eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.arg_f_p_of_psi(0.7, 0.3), -1.53567402, eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.arg_f_p_of_psi(0.0, 0.0), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.arg_f_p_of_psi(-1.0, 0.0), 0.00000000, eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.re_residue_rho(), -0.21726707, eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.im_residue_rho(), 0.35450739, eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.root_penalty(), 1.00000000, eps);
                }


                p["0->pipi::a_(+,0)^4@BHKMNR2026"]     = 0.30;
                p["0->pipi::a_(+,0)^5@BHKMNR2026"]     = -0.12;
                p["0->pipi::M_(+,0,0)@BHKMNR2026"]     = 0.782;
                p["0->pipi::Gamma_(+,0,0)@BHKMNR2026"] = 0.010;

                {
                    Options o{
                        { "n-resonances-I1"_ok,   "1"_ov },
                        { "n-resonances-I0"_ok,   "1"_ov },
                        {               "I"_ok, "1|0"_ov }
                    };
                    BHKMNR2026FormFactors<VacuumToPiPi> ff(p, o);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.Q(0.0)), 2.63724196, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.Q(0.0)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.Q(0.5)), 1.58692260, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.Q(0.5)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.Q(complex<double>(0.5, 0.5))), 0.96123390, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.Q(complex<double>(0.5, 0.5))), -1.26672457, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.Q(1.0)), 0.72405168, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.Q(1.0)), 0.00000000, eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.dQdpsi(0.0)), -0.01344323, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dQdpsi(0.0)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dQdpsi(0.5)), -2.52319095, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dQdpsi(0.5)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dQdpsi(complex<double>(0.5, 0.5))), -1.75330112, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dQdpsi(complex<double>(0.5, 0.5))), 3.12056474, eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms_I0(0, 0.0)), -0.01344323, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms_I0(0, 0.0)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms_I0(1, 0.0)), 2.63724196, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms_I0(1, 0.0)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms_I0(0, 0.5)), -2.52319095, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms_I0(0, 0.5)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms_I0(1, 0.5)), 0.32532713, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms_I0(1, 0.5)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms_I0(1, complex<double>(0.5, 0.5))), -1.47569903, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms_I0(1, complex<double>(0.5, 0.5))), -0.58309276, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_terms_I0(2, complex<double>(0.5, 0.5))), 0.66767610, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_terms_I0(2, complex<double>(0.5, 0.5))), -1.18214123, eps);

                    const auto constrained_a = ff.constrained_a_fp_I0();
                    TEST_CHECK_NEARLY_EQUAL(constrained_a[0], 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(constrained_a[1], -0.00237374, eps);
                    TEST_CHECK_NEARLY_EQUAL(constrained_a[2], 0.44717112, eps);
                    TEST_CHECK_NEARLY_EQUAL(constrained_a[3], 1.25746412, eps);


                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(0.0)), 1.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(0.0)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(0.5)), -8.71739159, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(0.5)), 13.7120901, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(complex<double>(0.5, 0.5))), 0.64534949, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(complex<double>(0.5, 0.5))), 0.43305815, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p(1.0)), -1.66223092, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p(1.0)), 3.41672815, eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.partial_wave(0.0)), -0.29680620, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.partial_wave(0.0)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.partial_wave(0.077919139600000)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.partial_wave(0.077919139600000)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.partial_wave(complex<double>(0.5, 0.5))), -0.05671173, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.partial_wave(complex<double>(0.5, 0.5))), 0.49602838, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.partial_wave(1.0)), -0.40967459, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.partial_wave(1.0)), 0.84208921, eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.scattering_length_parameters()[0]), -7.06958304, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.scattering_length_parameters()[0]), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.scattering_length_parameters()[1]), -2.98208882, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.scattering_length_parameters()[1]), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.scattering_length_parameters()[2]), 234.44727345, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.scattering_length_parameters()[2]), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.scattering_length_parameters()[3]), 1080.91461360, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.scattering_length_parameters()[3]), 0.00000000, eps);

                    // needed for charged pion radius
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_11(0.0)), -1.55707516, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_11(0.0)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_11(0.077919139600000)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_11(0.077919139600000)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_11(complex<double>(0.5, 0.5))), -0.31526741, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.dfdpsi_11(complex<double>(0.5, 0.5))), 2.19916563, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_11(1.0)), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.dfdpsi_11(1.0)), 0.00000000, eps);


                    TEST_CHECK_RELATIVE_ERROR(ff.saturation(), 85.25893239, eps);

                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p_of_psi(complex<double>(0.5, 0.5))), 0.11417128, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p_of_psi(complex<double>(0.5, 0.5))), -0.57836824, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(ff.f_p_of_psi(complex<double>(0.7, 0.3))), 0.08283404, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(ff.f_p_of_psi(complex<double>(0.7, 0.3))), -0.23132563, eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.abs2_f_p_of_psi(0.5, 0.5), 0.34754491, eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.abs2_f_p_of_psi(0.7, 0.3), 0.06037303, eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.abs2_f_p_of_psi(0.0, 0.0), 1.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.abs2_f_p_of_psi(-1.0, 0.0), 0.00000000, eps);

                    TEST_CHECK_NEARLY_EQUAL(ff.arg_f_p_of_psi(0.5, 0.5), -1.37589970, eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.arg_f_p_of_psi(0.7, 0.3), -1.22693782, eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.arg_f_p_of_psi(0.0, 0.0), 0.00000000, eps);
                    TEST_CHECK_NEARLY_EQUAL(ff.arg_f_p_of_psi(-1.0, 0.0), 0.00000000, eps);
                }
            }
        }
} parametric_BHKMNR2026_test;
