/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Danny van Dyk
 * Copyright (c) 2022 Stephan Kürten
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
#include <eos/b-decays/b-to-3l-nu.hh>

using namespace test;
using namespace eos;

class BToThreeLeptonsNeutrinoTest :
    public TestCase
{
    public:
        BToThreeLeptonsNeutrinoTest() :
            TestCase("b_to_3l_nu_test")
        {
        }

        virtual void run() const
        {

            Parameters p = Parameters::Defaults();

            {

                BToThreeLeptonsNeutrino three_lnu(p, Options{ });

                // Tested against an IPython/Jupyter notebook implementation of Stephan Kürten
                //TODO: do these tests for the branching ratio instead of decay width since I dont export decay widths anymore
                TEST_CHECK_NEARLY_EQUAL(4.4169722387470174e-11, three_lnu.double_differential_branching_ratio(1.2, 0.6), 1e-15);
                TEST_CHECK_NEARLY_EQUAL(4.4899433542534617e-11, three_lnu.double_differential_branching_ratio(1.2, 0.9), 1e-15);
                TEST_CHECK_NEARLY_EQUAL(1.9615874219144387e-11, three_lnu.double_differential_branching_ratio(1.5, 0.6), 1e-15);
                TEST_CHECK_NEARLY_EQUAL(1.991304494545849e-11 , three_lnu.double_differential_branching_ratio(1.5, 0.9), 1e-15);
                TEST_CHECK_NEARLY_EQUAL(7.978355662846529e-09 , three_lnu.double_differential_branching_ratio(0.6, 1.2), 1e-13);
                TEST_CHECK_NEARLY_EQUAL(1.6747424093160892e-10, three_lnu.double_differential_branching_ratio(0.9, 1.2), 1e-14);
                TEST_CHECK_NEARLY_EQUAL(8.170238677828805e-09 , three_lnu.double_differential_branching_ratio(0.6, 1.5), 1e-13);
                TEST_CHECK_NEARLY_EQUAL(1.7053307009828605e-10, three_lnu.double_differential_branching_ratio(0.9, 1.5), 1e-14);
                /*
                 * differential decay width of 5 kinematic variables
                 * quintuple_differential_decay_width(q2, k2, z_gamma, z_w, phi)
                 * q2 is the invariant mass of the off-shell photon in the range 4 m_l'^2 <= q2 <= (m_B - k2)^2
                 * k2 is the invariant mass of the W-meson in the range m_l^2 <= k2 <= (m_B -q2)^2
                 * z_gamma is the angle between the negatively charged lepton l' and the negative z-axis
                 * z_w is the angle between the charged lepton l and the positive z-axis
                 * phi is the angle between the q2 plane and the k2 plane
                 */
                TEST_CHECK_NEARLY_EQUAL(1.690117135105216e-12, three_lnu.quintuple_differential_branching_ratio(1.5, 0.9, 0.0, 0.0, 0.0), 1e-18);
                TEST_CHECK_NEARLY_EQUAL(1.992595975584582e-13, three_lnu.quintuple_differential_branching_ratio(1.5, 0.9, 1.0,-1.0, 3.1), 1e-18);
                TEST_CHECK_NEARLY_EQUAL(1.317143949587225e-14, three_lnu.quintuple_differential_branching_ratio(1.5, 0.9,-1.0, 1.0,-3.1), 1e-18);
                TEST_CHECK_NEARLY_EQUAL(1.992595975584582e-13, three_lnu.quintuple_differential_branching_ratio(1.5, 0.9,-1.0,-1.0,-3.1), 1e-18);
                TEST_CHECK_NEARLY_EQUAL(1.317143949587225e-14, three_lnu.quintuple_differential_branching_ratio(1.5, 0.9, 1.0, 1.0, 3.1), 1e-18);
                TEST_CHECK_NEARLY_EQUAL(8.840457847366234e-13, three_lnu.quintuple_differential_branching_ratio(1.5, 0.9, 0.6, 0.4, 3.1), 1e-18);
                TEST_CHECK_NEARLY_EQUAL(9.06242765770747e-13 , three_lnu.quintuple_differential_branching_ratio(1.5, 0.9, 0.4, 0.6, 0.3), 1e-18);
                TEST_CHECK_NEARLY_EQUAL(9.062269295937419e-13, three_lnu.quintuple_differential_branching_ratio(1.5, 0.9, 0.4, 0.6,-0.3), 1e-18);
                TEST_CHECK_NEARLY_EQUAL(2.050375005253621e-12, three_lnu.quintuple_differential_branching_ratio(1.2, 0.6, 0.4, 0.6,-0.3), 1e-18);
                TEST_CHECK_NEARLY_EQUAL(3.576321002331276e-10, three_lnu.quintuple_differential_branching_ratio(0.6, 1.2, 0.4, 0.6,-0.3), 1e-18);
                TEST_CHECK_NEARLY_EQUAL(2.543292200904078e-14, three_lnu.quintuple_differential_branching_ratio(5.0, 0.1, 0.4, 0.6,-0.3), 1e-18);
                TEST_CHECK_NEARLY_EQUAL(4.357725963221514e-12, three_lnu.quintuple_differential_branching_ratio(0.1, 5.0, 0.4, 0.6,-0.3), 1e-18);
                TEST_CHECK_NEARLY_EQUAL(4.687779112550517e-14, three_lnu.quintuple_differential_branching_ratio(4.0, 1.0, 0.4, 0.6,-0.3), 1e-18);

                TEST_CHECK_NEARLY_EQUAL(-0.03752014033347991, three_lnu.integrated_forward_backward_asymmetry(0.05, 4.2, 0.001, 1.0), 1e-5);
                TEST_CHECK_NEARLY_EQUAL(-0.01182330324675401, three_lnu.integrated_forward_backward_asymmetry(1.0 , 4.7, 0.001, 0.5), 1e-5);
                TEST_CHECK_NEARLY_EQUAL(-0.09288173467575606, three_lnu.integrated_forward_backward_asymmetry(0.1 , 3.7, 1.2  , 1.5), 1e-5);
                TEST_CHECK_NEARLY_EQUAL(-0.04209717482075583, three_lnu.integrated_forward_backward_asymmetry(0.45, 2.0, 0.001, 1.2), 1e-5);
            }
        }
} b_to_3l_nu_test;
